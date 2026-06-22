#!/usr/bin/env python

import os
from past.utils import old_div

import numpy as np
from scipy.stats.stats import pearsonr
from scipy.optimize import curve_fit
import logging
logger = logging.getLogger(__name__)

import opendrift.models.openoil.tamoc_chemical_properties as chem
import adios_db.computation.physical_properties as props  # OD 1.14.9: moved to adios_db
import adios_db.computation.gnome_oil as goil             # OD 1.14.9: moved to adios_db
import opendrift.models.openoil.cicoil_estimations as est  # CICOIL util (kept as standalone)

def _riazi(T, A, B):
    # for non-linear curve fit function
    return 1. - np.exp(-B * ((T - 309.2)/309.2) ** B / A)

class FluidProps(object):
    chemdata_file = 'pseudo_chemdata.csv'
    gas_comp_file = 'gas_composition.csv'
    pseudo_T = np.array([423.15, 453.15, 473.15, 503.15, 553.15, 573.15, 623.15, 653.15])
    sara_T = np.array([423.15, 453.15, 473.15, 503.15, 553.15, 573.15, 623.15, 653.15, 723.0])
    # create the required arrays to store the pseudo-component properties
    component_density = np.array([706.9, 735.3, 740.3, 743.7, 747.5, 783.9, 797.1, 805.2, 855.3,   # Only for dead oil
                                  844.1, 869.5, 829.2, 962.9, 901.3, 989.9, 990.9, 1040.8])
    component_density_Tb = np.array([633.74, 587.84, 568.82, 551.16, 513.58, 484.78, 461.47, 432.49,
                                     747.60, 682.98, 683.05, 609.79, 673.28, 577.38, 565.96, 512.09, 379.42])

    def __init__(self, oiltype, max_cuts):
        __location__ = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

        # Create the full relative path to the default data in ChemData.csv
        chem_path = os.path.join(__location__, self.chemdata_file)
        gas_comp_path = os.path.join(__location__, self.gas_comp_file)
        data, units = chem.load_data(chem_path)
        composition = []
        comp_Mw = []
        comp_Tb = []
        gas_mol_frac = []
        with open(gas_comp_path) as datfile:    # gas composition
            for line in datfile:
                # Get a line of data
                entries = line.strip().split(',')
                composition.append(entries[0])
                comp_Mw.append(data[entries[0]]['M'])
                comp_Tb.append(data[entries[0]]['Tb'])
                gas_mol_frac.append(np.float64(entries[1]))
        self.gas_mol_frac = np.array(gas_mol_frac)
        self.comp_Mw = np.array(comp_Mw) * 1000.
        self.comp_Tb = np.array(comp_Tb)
        self.composition = composition
        self.data = data
        self.units = units
        self.idx_sat = len(self.gas_mol_frac) - 2 * len(self.pseudo_T) - 1
        self.idx_arom = len(self.gas_mol_frac) - len(self.pseudo_T) - 1

        self.component_molar_mass = self.comp_Mw[self.idx_sat:] / 1000
        self.component_boil_temp = self.comp_Tb[self.idx_sat:]
        self.component_no = len(composition)
        self.Cox_A0 = 1.77838e-3 * self.component_boil_temp + 2.186529
        self.Cox_A1 = -1.436577e-7 * self.component_boil_temp - 1.908404e-3
        self.Cox_A2 = -2.734678e-9 * self.component_boil_temp + 2.858042e-6
        self.density = props.Density(oiltype).at_temp(288.15, 'K')
        self.viscosity = props.KinematicViscosity(oiltype).at_temp(288.15, temp_units='K')
        self.dyn_visc = self.viscosity * self.density * 1000.       # dynamic viscosity in mPa*s
        self.record = oiltype
        self._create_TBP(max_cuts)

    def setup_treatment(self, treatment, dispersant, DOR, disp_stock, delay, time_frac_surf, time_frac_ssdi,
                                           effic_surf, effic_ssdi, wind_limit, equipment_rate):
        self.treatment = treatment
        self.dispersant = dispersant
        self.DOR = DOR
        self.disp_stock = disp_stock
        self.delay = delay
        self.time_frac_surf = time_frac_surf
        self.time_frac_ssdi = time_frac_ssdi
        self.effic_surf = effic_surf
        self.effic_ssdi = effic_ssdi
        self.wind_limit = wind_limit
        self.equipment_rate = equipment_rate

    def _create_TBP(self, max_cuts):
        # Calculate pseudo-components properties
        saturates = self.record.sub_samples[0].SARA.saturates
        aromatics = self.record.sub_samples[0].SARA.aromatics
        resins = self.record.sub_samples[0].SARA.resins
        asphaltenes = self.record.sub_samples[0].SARA.asphaltenes

        if saturates is None:
            self.fr_satur = est.saturates_fraction(self.density / 1000, self.dyn_visc)          # correlation uses empirical units density: g/cm3, d.visc: mPa*s
        else:
            self.fr_satur = saturates.converted_to('fraction').value
        if aromatics is None:
            self.fr_arom = est.aromatics_fraction(self.density / 1000, self.dyn_visc, self.fr_satur)
        else:
            self.fr_arom = aromatics.converted_to('fraction').value
        if resins is None:
            self.fr_resin = est.resin_fraction(self.density / 1000, self.dyn_visc)
        else:
            self.fr_resin = resins.converted_to('fraction').value
        if asphaltenes is None:
            self.fr_asphalt = est.asphaltene_fraction(self.density / 1000, self.dyn_visc)
        else:
            self.fr_asphalt = asphaltenes.converted_to('fraction').value

        pseudo_mass = np.zeros(np.shape(self.sara_T))
        if max_cuts > 0:
            cuts = props.get_distillation_cuts(self.record)
            inner_cuts = props.get_distillation_cuts(self.record)
            cum_mass, cut_bp = zip(*[c for c in cuts])
            inner_cuts = [c for c in inner_cuts if c[1] < 775 and c[1] > 375]        # c.vapor_temp_k > 450 and
            logger.debug('Number of inner cuts: %s' % (len(inner_cuts)))
            inner_cum_mass, inner_cut_bp = zip(*[c for c in inner_cuts])
            cut_bp = np.array(cut_bp)
            cum_mass = np.array(cum_mass)
            inner_cut_bp = np.array(inner_cut_bp)
            inner_cum_mass = np.array(inner_cum_mass)

            # Compute the initial boiling temperature
            t_ini = cut_bp[0] - (cut_bp[1] - cut_bp[0]) / (cum_mass[1] - cum_mass[0]) * cum_mass[0]
            # Riazi TBP function
            LnLn_frac = np.log(np.log(1. / (1. - inner_cum_mass)))
            Ln_temp = np.log((inner_cut_bp - t_ini) / t_ini)

            # Solve with least square method to find slope and constant
            X = np.vstack([LnLn_frac, np.ones(len(LnLn_frac))]).T
            slope, bias = np.linalg.lstsq(X, Ln_temp, rcond=None)[0]
            R_coef =   pearsonr(LnLn_frac, Ln_temp)[0]**2
            self.R_slope = slope
            self.R_bias = bias
            self.R_Tini = t_ini
            logger.debug('Riazi slope: %s   Riazi Constant: %s   R^2: %s' % (slope, bias, R_coef))
            logger.debug('LnLn_frac: %s   Ln_time: %s' % (LnLn_frac, Ln_temp))

            Ln_pseudoT = np.log((self.pseudo_T - t_ini) / t_ini)
            LnLn_pseudoM = (Ln_pseudoT - bias) / slope
            pseudo_cum_mass = 1. - 1. / np.exp(np.exp(LnLn_pseudoM))
            pseudo_cum_mass[np.isnan(pseudo_cum_mass)] = 0.
            res_start = pseudo_cum_mass[-1]
            bp_res = np.exp(np.log(np.log(1. / (1. - (res_start +
                                                      (1. - self.fr_asphalt - self.fr_resin - res_start) / 2))))
                            * slope + bias) * t_ini + t_ini
            bp_50 = np.exp(np.log(np.log(1. / (1. - 0.5))) * slope + bias) * t_ini + t_ini
            logger.debug('Riazi Median Bp_50 %s  Riazi before resin-Asph Bp: %s' % (bp_50, bp_res))
            logger.debug('Residue Bp: %s' % (bp_res))
            for i in range(np.min([len(self.pseudo_T), len(inner_cum_mass)])):
                logger.debug('Actual cut cumulative mass:  %s  %s  Riazi pseudo cum mass:  %s  %s' % (inner_cut_bp[i], inner_cum_mass[i], self.pseudo_T[i], pseudo_cum_mass[i]))
        else:
            # cut_bp, cum_mass = self.record.normalized_cut_values(N=15)
            cut_bp, cut_mass = goil.normalized_cut_values_james(self.record, N=15)
            cum_mass = np.cumsum(cut_mass)
            popt, _pcov = curve_fit(goil._linear_curve, cut_bp, cum_mass)
            f_cutoff = goil._linear_curve(732.0, *popt)  # center of asymptote (< 739)
            popt_1 = popt.tolist() + [f_cutoff]
            pseudo_cum_mass = goil._linear_curve(self.pseudo_T, *popt)

            # median 50% mass boiling point
            bp_50 = goil._inverse_linear_curve(0.5, *popt_1)
            bp_res = np.average([self.sara_T[-2], goil._inverse_linear_curve(1 - self.fr_asphalt - self.fr_resin - 0.05, *popt_1)])
            logger.debug('Linear Median Bp_50 %s  Linear before resin-Asph Bp: %s' % (bp_50, bp_res))

            ########################
            inner_bp, inner_cuts = goil.normalized_cut_values_james(self.record, N=15)
            inner_cut_bp = []
            inner_cum_mass = []
            inner_cuts = np.cumsum(inner_cuts)

            for i in range(len(inner_bp)):
                if inner_bp[i] > 450 and inner_bp[i] < 775:
                    inner_cut_bp.append(inner_bp[i])
                    inner_cum_mass.append(inner_cuts[i])
            inner_cut_bp = np.array(inner_cut_bp)
            inner_cum_mass = np.array(inner_cum_mass)

            # calculate the Riazi function with curve fit
            Rpopt, _Rpcov = curve_fit(_riazi, inner_cut_bp, inner_cum_mass, bounds=([0., 1.], [10., 2.]))
            fit_A = Rpopt[0]
            fit_B = Rpopt[1]
            fit_To = 309.2
            fit_cum_mass = 1. - np.exp(-fit_B * ((self.pseudo_T - fit_To) / fit_To) ** fit_B / fit_A)
            logger.debug('non-linear fit results: %s Riazi A: %s Riazi B: %s' % (Rpopt, Rpopt[0], Rpopt[1]))

            self.R_bias = np.log(Rpopt[0] / Rpopt[1]) / Rpopt[1]
            self.R_slope = 1 / Rpopt[1]
            self.R_Tini = fit_To
            logger.debug('Slope: %s Bias: %s' % (self.R_slope, self.R_bias))

            #######################

            for i in range(len(self.pseudo_T)):
                logger.debug('Estimated cuts from API cut cumulative mass: %s  %s  Linear pseudo cumulative mass:  %s %s  Riazi fit: %s' % (cut_bp[i], cum_mass[i], self.pseudo_T[i], pseudo_cum_mass[i], fit_cum_mass[i]))

        # Saturates Aromatics per cut
        # cum_mass[1:] = np.diff(cum_mass)
        self.bp_res = bp_res
        self.sara_T[-1] = bp_res
        pseudo_mass[0] = pseudo_cum_mass[0]
        pseudo_mass[1:-1] = np.diff(pseudo_cum_mass)
        pseudo_mass[-1] = 1 - pseudo_cum_mass[-1] - self.fr_resin - self.fr_asphalt
        pseudo_sat = pseudo_mass / 2.
        pseudo_aro = pseudo_mass / 2.
        logger.debug('Pseudo masses: ' + str(pseudo_mass))
        non_zero_indices = np.where(pseudo_mass != 0.)[0]
        for _i in range(20):
            pseudo_sat[non_zero_indices], pseudo_aro[non_zero_indices] = \
                self.verify_cut_fractional_masses(pseudo_mass[non_zero_indices], self.sara_T[non_zero_indices],
                                                         pseudo_sat[non_zero_indices], pseudo_aro[non_zero_indices])
        # Multiply the cuts saturates and cuts aromatics with a corrections factor to match the SARA fractions
        pseudo_sat = pseudo_sat * self.fr_satur / np.sum(pseudo_sat)
        pseudo_aro = pseudo_aro * self.fr_arom / np.sum(pseudo_aro)
        pseudo_mass_corr = np.ones(len(pseudo_mass))

        pseudo_mass_corr[non_zero_indices] = pseudo_mass[non_zero_indices] / (pseudo_sat[non_zero_indices] + pseudo_aro[non_zero_indices])
        self.pseudo_sat = pseudo_sat * pseudo_mass_corr
        self.pseudo_aro = pseudo_aro * pseudo_mass_corr
        self.pseudo_mass = pseudo_mass

    def get_Riazi_model(self):
        return self.R_bias, self.R_slope, self.R_Tini

    def get_live_composition(self, GOR=0):
        # Calculate detailed pseudo components Mw and moles
        idx_arom = len(self.gas_mol_frac) - len(self.pseudo_aro)
        idx_sat = idx_arom - len(self.pseudo_sat) + 1

        arom_Mw = np.zeros(np.shape(self.sara_T))
        arom_Mw[0:-1] = self.comp_Mw[idx_arom:-1]
        arom_Mw[-1] = est.aromatic_mol_wt(self.sara_T[-1])
        sat_Mw = np.zeros(np.shape(self.sara_T))
        sat_Mw[0:-1] = self.comp_Mw[idx_sat:idx_arom]
        sat_Mw[-1] = est.saturate_mol_wt(self.sara_T[-1])
        arom_mol = self.pseudo_aro / arom_Mw
        sat_mol = self.pseudo_sat / sat_Mw
        resin_mol = self.fr_resin / est.resin_mol_wt(800)
        asphalt_mol = self.fr_asphalt / est.asphaltene_mol_wt(1000)
        total_mol = np.sum(sat_mol) + np.sum(arom_mol) + resin_mol + asphalt_mol
        oil_Mw = (np.sum(self.pseudo_sat) + np.sum(self.pseudo_aro) + self.fr_resin + self.fr_asphalt) / total_mol
        oil_d = props.Density(self.record).at_temp(288.15, 'K') / 1000.

        logger.debug('Saturates pseudo-components Mw: ' + str(sat_Mw))
        logger.debug('Aromatic pseudo-components Mw: ' + str(arom_Mw))
        logger.debug('Oil Mw from detailed pseudo-components: ' + str(oil_Mw))
        Xg = 0
        if GOR != 0: Xg = (1 + oil_d / (7.521 * 10**-6 * GOR * oil_Mw))**-1
        logger.debug('Gas / Oil Xg: ' + str(Xg))

        res_mass = self.pseudo_mass[-1] + self.fr_resin + self.fr_asphalt
        res_mol = arom_mol[-1] + sat_mol[-1] + resin_mol + asphalt_mol
        res_mw = res_mass / res_mol
        self.comp_Mw[-1] = res_mw
        oil_mol_frac = np.zeros(np.shape(self.gas_mol_frac))
        logger.debug('Index arom, index sat: %s %s' % (idx_arom, idx_sat))
        oil_mol_frac[idx_sat:idx_arom] = sat_mol[0:-1]
        oil_mol_frac[idx_arom:-1] = arom_mol[0:-1]
        oil_mol_frac[-1] = res_mol
        oil_mol_frac = oil_mol_frac / (np.sum(sat_mol[0:-1]) + np.sum(arom_mol[0:-1]) + res_mol)
        self.live_oil_mol = Xg * self.gas_mol_frac + (1. - Xg) * oil_mol_frac
        self.oil_mol_frac = oil_mol_frac[idx_sat:]
        self.oil_mass_frac = self.oil_mol_frac * self.comp_Mw[idx_sat:] / \
                             np.sum(self.oil_mol_frac * self.comp_Mw[idx_sat:])
        logger.debug('Live oil mole composition: ' + str(self.live_oil_mol))
        logger.debug('Dead oil mole composition: ' + str(self.oil_mol_frac))
        logger.debug('Dead oil mass composition: ' + str(self.oil_mass_frac))

        density = props.Density(self.record).at_temp(288.15, 'K')
        viscosity = props.KinematicViscosity(self.record).at_temp(288.15, temp_units='K')
        f_sat = est.saturates_fraction(density, viscosity)
        f_aro = 1. - f_sat - self.fr_resin - self.fr_asphalt
        #est_cuts_sara = record.component_mass_fractions_riazi()

        logger.debug('Measured Aromatics fraction: %s  Estimated Aromatic fraction: %s  Estimated cumulative Aromatics from cuts: %s' % (self.fr_arom, f_aro, np.sum(self.pseudo_aro)))
        logger.debug('Measured Saturates fraction: %s  Estimated Saturates fraction: %s  Estimated cumulative Saturates from cuts: %s' % (self.fr_satur, f_sat, np.sum(self.pseudo_sat)))
        logger.debug('All cut aroms: ' + str(self.pseudo_aro))
        logger.debug('All cut sats: ' + str(self.pseudo_sat))
        logger.debug('pseudo mass balance plus resins, asphaltenes' + str(np.sum(self.pseudo_mass) + self.fr_resin + self.fr_asphalt))
        return self.live_oil_mol, self.composition, self.data, self.units

    @classmethod
    def verify_cut_fractional_masses(cls, fmass_i, T_i, f_sat_i, f_arom_i,
                                     prev_f_sat_i=None):
        '''
            Assuming a distillate mass with a boiling point T_i,
            We propose what the component fractional masses might be.

            We calculate what the molecular weights and specific gravities
            likely are for saturates and aromatics at that temperature.

            Then we use these values, in combination with our proposed
            component fractional masses, to produce a proposed average
            molecular weight and specific gravity for the distillate.

            We then use Riazi's formulas (3.77 and 3.78) to obtain the
            saturate and aromatic fractional masses that represent our
            averaged molecular weight and specific gravity.

            If our proposed component mass fractions were correct (or at least
            consistent with Riazi's findings), then our computed component
            mass fractions should match.

            If they don't match, then the computed component fractions should
            at least be a closer approximation to that which is consistent
            with Riazi.

            It is intended that we run this function iteratively to obtain a
            successively approximated value for f_sat_i and f_arom_i.
        '''
        assert np.allclose(fmass_i, f_sat_i + f_arom_i)

        M_w_sat_i = est.saturate_mol_wt(T_i)
        M_w_arom_i = est.aromatic_mol_wt(T_i)

        M_w_avg_i = (old_div(M_w_sat_i * f_sat_i, fmass_i) +
                     old_div(M_w_arom_i * f_arom_i, fmass_i))

        # estimate specific gravity
        rho_sat_i = est.saturate_densities(T_i)
        SG_sat_i = est.specific_gravity(rho_sat_i)

        rho_arom_i = est.aromatic_densities(T_i)
        SG_arom_i = est.specific_gravity(rho_arom_i)

        SG_avg_i = (old_div(SG_sat_i * f_sat_i, fmass_i) +
                    old_div(SG_arom_i * f_arom_i, fmass_i))

        f_sat_i = est.saturate_mass_fraction(fmass_i, M_w_avg_i, SG_avg_i, T_i)
        f_arom_i = fmass_i - f_sat_i

        # Note: Riazi states that eqs. 3.77 and 3.78 only work with
        #       molecular weights less than 200. In those cases,
        #       Chris would like to use the last fraction in which the
        #       molecular weight was less than 200 instead of just guessing
        #       50/50
        # TODO: In the future we might be able to figure out how
        #       to implement CPPF eqs. 3.81 and 3.82, which take
        #       care of cases where molecular weight is greater
        #       than 200.
        above_200 = M_w_avg_i > 200.0
        try:
            if np.any(above_200):
                if np.all(above_200):
                    # once in awhile we get a record where all molecular
                    # weights are over 200, In this case, we have no
                    # choice but to use the 50/50 scale
                    scale_sat_i = 0.5
                else:
                    last_good_sat_i = f_sat_i[above_200 ^ True][-1]
                    last_good_fmass_i = fmass_i[above_200 ^ True][-1]

                    scale_sat_i = old_div(last_good_sat_i, last_good_fmass_i)

                f_sat_i[above_200] = fmass_i[above_200] * scale_sat_i
                f_arom_i[above_200] = fmass_i[above_200] * (1.0 - scale_sat_i)
        except TypeError:
            # numpy array assignment failed, try a scalar assignment
            if above_200:
                # for a scalar, the only way to determine the last
                # successfully computed f_sat_i is to pass it in
                if prev_f_sat_i is None:
                    scale_sat_i = 0.5
                else:
                    scale_sat_i = old_div(prev_f_sat_i, fmass_i)

                f_sat_i = fmass_i * scale_sat_i
                f_arom_i = fmass_i * (1.0 - scale_sat_i)

        return f_sat_i, f_arom_i