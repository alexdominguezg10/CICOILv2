from tamoc import ambient, seawater
from tamoc import dbm
from tamoc import dispersed_phases
from tamoc import bent_plume_model
from tamoc import sintef
from tamoc import model_share
import tamoc.chemical_properties as chem
import numpy as np
import logging
logger = logging.getLogger(__name__)


def run_tamoc(live_comp, composition, chemdata, chemunits, output_name, profile, total_flow, z0, D, Tj, phi_0,
              theta_0, bins, *args, **kwargs):

    if 'tracking' in kwargs:
        tracking = kwargs['tracking']
    else:
        tracking = False
    # load chem data and composition
    data = chemdata
    units = chemunits
    composition = composition
    fluid_mol_frac = live_comp
    fluid = dbm.FluidMixture(composition, user_data=data)

    # Insert a constant crossflow velocity
    # Open an ambient profile object from the netCDF dataset
    v_profile = ambient.Profile(profile, chem_names='all')

    # Jet initial conditions
    U0 = 0.
    Sj = 0.
    cj = 1.

    chem_name = 'tracer'
    masses = fluid.mass_frac(fluid_mol_frac) * total_flow

    # Create the gas phase particles
    droplet = dbm.FluidParticle(composition, fp_type=1, user_data=data)
    bubble = dbm.FluidParticle(composition, fp_type=0, user_data=data)
    disp_phases = []

    # formation of gas and liquid
    T, S, P = v_profile.get_values(z0, ['temperature', 'salinity', 'pressure'])
    logger.debug('temperature=%s salinity=%s pressure=%s' % (T, S, P))
    masses0, xi, K = fluid.equilibrium(masses, T, P)
    rho_g = fluid.density(masses0[0, :], T, P)[0, 0]
    rho_l = fluid.density(masses0[1, :], T, P)[1, 0]
    rho_sw = seawater.density(T, S, P)
    logger.debug('seawater density check: ' + str(seawater.density(273.15+4, 35, 25000000)))
    gas_mol_frac = fluid.mol_frac(masses0[0, :])
    liq_mol_frac = fluid.mol_frac(masses0[1, :])

    # droplet-bubble size distribution
    mu_gas = fluid.viscosity(masses0[0, :], Tj, P)[0, 0]
    mu_liq_correction = 1.
    sigma_liq_correction = 0.5
    mu_liq = fluid.viscosity(masses0[1, :], Tj, P)[1, 0] * mu_liq_correction
    sigma_gas = fluid.interface_tension(masses0[0, :], Tj, S, P)[0, 0]
    sigma_liq = fluid.interface_tension(masses0[1, :], Tj, S, P)[1, 0] * sigma_liq_correction

    if 'dispersant' in kwargs:
        dispersant = kwargs['dispersant']
        applied_DOR = kwargs['applied_DOR']
        sigma_liq = sigma_liq * (dispersant['A'] * 0.5 ** (applied_DOR * 100. / dispersant['t12']) + dispersant['B']) / (
                           dispersant['A'] + dispersant['B'])

    logger.debug('mu gas: %s     mu liq: %s' % (mu_gas, mu_liq))
    logger.debug('sig_gas: %s    sig liq: %s' % (sigma_gas, sigma_liq))

    # Johansen 2013 model
    d50_gas, d50_liq = sintef.modified_We_model(D, rho_g, masses0[0, :], mu_gas,
                                                sigma_gas, rho_l, masses0[1, :], mu_liq, sigma_liq, rho_sw)
    nbins = bins
    bubble_de, bubble_md0 = sintef.rosin_rammler(nbins, d50_gas, np.sum(masses0[0, :]),
                                                 sigma_gas, rho_g, rho_sw)
    droplet_de, droplet_md0 = sintef.rosin_rammler(nbins, d50_liq, np.sum(masses0[1, :]),
                                                   sigma_liq, rho_l, rho_sw)

    bubble_lambda_1 = np.linspace(0.93, 0.85, nbins)
    droplet_lambda_1 = np.linspace(0.95, 0.89, nbins)
    logger.debug('live bubble sizes: ' + str(bubble_de))
    logger.debug('live bubble masses: ' + str(bubble_md0))
    fraction_sum = 0.
    for i in range(len(droplet_de)):
        fraction_sum += droplet_md0[i]/np.sum(droplet_md0)
        logger.debug('live droplet size: %s Droplet mass fraction: %s Cum. mass fraction: %s' % (droplet_de[i], droplet_md0[i]/np.sum(droplet_md0), fraction_sum))


    for i in range(len(bubble_de)):
        if bubble_md0[i] ==0.: continue
        (m0, T0, nb0, P0, Sa, Ta) = dispersed_phases.initial_conditions(v_profile, z0, bubble, gas_mol_frac, bubble_md0[i], 2,
                                                                       bubble_de[i], Tj)
        t_hyd = dispersed_phases.hydrate_formation_time(bubble, z0, m0, T0, v_profile)
        bpm_particle = bent_plume_model.Particle(0., 0., z0, bubble, m0, T0, nb0, bubble_lambda_1[i], P0, Sa, Ta,
                                                 K=1., K_T=1., fdis=1.e-6, t_hyd=t_hyd)
        disp_phases.append(bpm_particle)
        logger.debug('initial t_hyd: ' + str(t_hyd))
    for i in range(len(droplet_de)):
        if droplet_md0[i] == 0.: continue
        (m0, T0, nb0, P0, Sa, Ta) = dispersed_phases.initial_conditions(v_profile, z0, droplet, liq_mol_frac, droplet_md0[i], 2,
                                                                       droplet_de[i], Tj)
        t_hyd = dispersed_phases.hydrate_formation_time(droplet, z0, m0, T0, v_profile)
        bpm_particle = bent_plume_model.Particle(0., 0., z0, droplet, m0, T0, nb0, droplet_lambda_1[i], P0, Sa, Ta,
                                                     K=1., K_T=1., fdis=1.e-6, t_hyd=t_hyd)
        disp_phases.append(bpm_particle)

    # Create the bent plume model object
    bpm = bent_plume_model.Model(v_profile)

    # Run the simulation
    bpm.simulate(np.array([0., 0., z0]), D, U0, phi_0, theta_0, Sj, Tj, cj, chem_name, disp_phases, track=tracking,
                 dt_max=60., sd_max=2000.)

    # Save the plume results
    plume_file = ''.join((output_name, '_plume.nc'))
    bpm.save_sim(plume_file, profile, 'HYCOM profile data in file v_profile.nc')

    # Save the particles results
    title = 'TAMOC Bent Plume Model particles results'
    particles_file = ''.join((output_name, '_particles.nc'))
    nc = model_share.tamoc_nc_file(particles_file, title,
                                   'Near-field exiting particles properties to be used in Far-Field model',
                                   'HYCOM profile data in file v_profile.nc',)
    bpm.q_local.update(bpm.t[-1], bpm.q[-1, :], bpm.profile, bpm.p, bpm.particles)
    dispersed_phases.save_particle_to_nc_file(nc, composition, bpm.particles, bpm.K_T0)

    return bpm