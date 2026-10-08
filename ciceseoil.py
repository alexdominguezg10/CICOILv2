# This file is part of OpenDrift.
#
# OpenDrift is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by of the GNU General Public License as published by
# the Free Software Foundation, version 2
#
# OpenDrift is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OpenDrift.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2015, Knut-Frode Dagestad, MET Norway
# This file was modified by Konstantinos Kotzakoulakis, CICESE Physical Oceanography/SINTEF Ocean

from __future__ import division

from io import open
import os
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpt
import cartopy.crs as ccrs
import logging

logger = logging.getLogger(__name__)

from . import noaa_oil_weathering as noaa
from . import adios
import adios_db.computation.physical_properties as props  # OD 1.14.9: moved to adios_db package
from opendrift.models.physics_methods import oil_wave_entrainment_rate_li2017

from .openoil import Oil, OpenOil
from .fluid_properties import FluidProps
# OD 1.14.9: config level constants moved from class attrs to module-level
from opendrift.config import CONFIG_LEVEL_ESSENTIAL, CONFIG_LEVEL_BASIC, CONFIG_LEVEL_ADVANCED

try:
    from itertools import izip as zip
except ImportError:
    pass

# Defining the oil element properties
class CiceseOil(Oil):
    """Extending Oil with variables relevant for multi pseudo-component oil particles."""

    variable_list = [
        ('spillet_area', {'dtype': np.float32,  # 1-Side surface area of the spillet
                          'units': 'm2',
                          'default': 0.114}),
        ('oil_molar_mass', {'dtype': np.float32,
                            'units': 'kg/mol',
                            'default': 0.2}),
        ('spillet_thickness', {'dtype': np.float32,
                               'units': 'm',
                               'default': 0.01}),
        ('oil_density', {'dtype': np.float32,
                         'units': 'kg/m^3',
                         'default': 880}),
        ('DOR', {'dtype': np.float32,           # Dispersant to Oil Ratio
                         'units': 'unitless',
                         'default': 0}),
        ('IFT', {'dtype': np.float32,           # Oil - Water interfacial tension
                 'units': 'Nm',
                 'default': 30.}),
        ('max_water', {'dtype': np.float32,     # Maximum water content
                 'units': 'unitless',
                 'default': 0}),
        ('mass_dissolved_surface', {'dtype': np.float32,
                                    'units': 'kg',
                                    'default': 0}),
        ('mass_dissolved_subsea', {'dtype': np.float32,
                                   'units': 'kg',
                                   'default': 0}),
        ('is_bubble', {'dtype': np.uint8,
                       'units': 'unitless',
                       'default': 0}),
        ('is_dissolved', {'dtype': np.uint8,
                          'units': 'unitless',
                          'default': 0}),
        ('mass_gas', {'dtype': np.float32,
                      'units': 'kg',
                      'default': 0}),
        ('mass_biodegraded_from_oil', {'dtype': np.float32,
                                       'units': 'kg',
                                       'default': 0}),
        ('mass_biodegraded_from_water', {'dtype': np.float32,
                                         'units': 'kg',
                                         'default': 0}),
        ('fraction_dissolved_surface', {'dtype': np.float32,
                                        'units': '%',
                                        'default': 0}),
        ('fraction_biodegraded_from_oil', {'dtype': np.float32,
                                           'units': '%',
                                           'default': 0}),
        ('fraction_biodegraded_from_water', {'dtype': np.float32,
                                             'units': '%',
                                             'default': 0}),
        # Photooxidation tracking (added for CICOIL+OD 1.14.9 integration)
        ('mass_photooxidized', {'dtype': np.float32,
                                'units': 'kg',
                                'default': 0}),
        ('fraction_photooxidized', {'dtype': np.float32,
                                    'units': '%',
                                    'default': 0}),
        # Oxygenated photoproduct (OP) tracking (audit improvement #4):
        # a fraction of newly-formed photooxidized mass partitions into a
        # dissolved OP pool instead of remaining as an oil-phase residue,
        # and decays at its own (faster) removal rate.
        ('mass_op_dissolved', {'dtype': np.float32,
                               'units': 'kg',
                               'default': 0}),
        ('fraction_op_dissolved', {'dtype': np.float32,
                                   'units': '%',
                                   'default': 0}),
        ('mass_op_degraded', {'dtype': np.float32,
                              'units': 'kg',
                              'default': 0}),
        ('fraction_op_degraded', {'dtype': np.float32,
                                  'units': '%',
                                  'default': 0}),
    ]
    components = ['comp_ali_G1', 'comp_ali_G2', 'comp_ali_G3', 'comp_ali_G4', 'comp_ali_G5', 'comp_ali_G6',
                  'comp_ali_G7',
                  'comp_ali_G8', 'comp_arom_G1', 'comp_arom_G2', 'comp_arom_G3', 'comp_arom_G4', 'comp_arom_G5',
                  'comp_arom_G6', 'comp_arom_G7', 'comp_arom_G8', 'comp_residue']

    #masses = ['mass_', 'evaporated_', 'dissolved_', 'biodeg_in-oil_', 'biodeg_in-water_']

    for pseudo in components:
        variable_list.append((pseudo, {'dtype': np.float32, 'units': 'kg', 'default': 0}))
        variable_list.append((pseudo + '_evaporated', {'dtype': np.float32, 'units': 'kg', 'default': 0}))
        variable_list.append((pseudo + '_dissolved', {'dtype': np.float32, 'units': 'kg', 'default': 0}))
    variables = Oil.add_variables(variable_list)


class OpenCiceseOil(OpenOil):
    """Open source oil trajectory model based on the OpenDrift framework.

        Developed at MET Norway based on oil weathering parameterisations
        found in open/published litterature.

        Under construction.
    """

    ElementType = CiceseOil

    required_variables = {
        'x_sea_water_velocity': {
            'fallback': None
        },
        'y_sea_water_velocity': {
            'fallback': None
        },
        'x_wind': {
            'fallback': None
        },
        'y_wind': {
            'fallback': None
        },
        'upward_sea_water_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_significant_height': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_stokes_drift_x_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_stokes_drift_y_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_period_at_variance_spectral_density_maximum': {
            'fallback': 0,
            'important': False
        },
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment':
            {
                'fallback': 0,
                'important': False
            },
        'sea_ice_area_fraction': {
            'fallback': 0,
            'important': False
        },
        'sea_ice_x_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_ice_y_velocity': {
            'fallback': 0,
            'important': False
        },
        'sea_water_temperature': {
            'fallback': 10,
            'profiles': True
        },
        'sea_water_salinity': {
            'fallback': 34,
            'profiles': True
        },
        'sea_floor_depth_below_sea_level': {
            'fallback': 10000
        },
        'ocean_vertical_diffusivity': {
            'fallback': 0.02,
            'important': False,
            'profiles': True
        },
        'land_binary_mask': {
            'fallback': None
        },
        'ocean_mixed_layer_thickness': {
            'fallback': 50,
            'important': False
        },
        'surface_downward_uv_radiation': {
            'fallback': -1.0,
            'important': False
        },
        # OD 1.14.9: oceandrift.vertical_mixing() reads
        # self.environment.sea_surface_height unconditionally; without this
        # entry (present in OpenOil.required_variables but dropped by this
        # override) it raises AttributeError as soon as vertical_mixing=True.
        'sea_surface_height': {
            'fallback': 0
        },
    }

    # OD 1.14.9: required_profiles_z_range is obsolete.
    # Profile depth is now set via config 'drift:profile_depth' (see __init__).

    max_speed = 1.3  # m/s

    # ── CICESE tunable constants registry ──────────────────────────────────
    # Documented values and literature sources for all physics constants used
    # in CICESE weathering mode.  Methods read from self.ATTR_NAME (class
    # attributes defined alongside their docstrings below); this dict is a
    # single-point reference for inspection, logging, and optional overrides.
    CICESE_DEFAULTS = {
        # Emulsification (Fingas & Fieldhouse 2004, NOAA Bullwinkle)
        'EMUL_STABILITY_NU_REF': 1000.0,   # cSt, reference viscosity
        'EMUL_STABILITY_MIN':    0.3,       # min stability clip factor
        'EMUL_STABILITY_MAX':    3.0,       # max stability clip factor
        'EMUL_RATE_SCALE':       0.0142,    # overall water-uptake rate correction

        # Biodegradation (Prince et al. 2017, Bacosa et al. 2015)
        'BIODEG_T_REF':          20.0,      # C, reference temperature for Q10
        'BIODEG_Q10':            3.0,       # Q10 temperature correction factor
        'BIODEG_RESIDUE':        0.00139,   # d^-1, residue rate (~500 d t1/2)
        'BIODEG_SURFACE_FACTOR': 1./30.,    # surface slick suppression (Prince et al.)
        'BIODEG_DROPLET_D_REF':  0.001,     # m, reference droplet diameter
        'BIODEG_DROPLET_EXPONENT': 1.0,     # (d_ref/d)^n exponent
        'BIODEG_DROPLET_DIAMETER_MIN': 1e-5,  # m (10 um)
        'BIODEG_DROPLET_DIAMETER_MAX': 1e-2,  # m (1 cm)
        'K_BIODEG_WATER_ANCHOR': 0.0578,    # d^-1, t1/2~12d @20C (Adcroft+ 2010)
        'BIODEG_WATER_DEPLETION_FRAC': 0.01,

        # Photooxidation (recal 2026-06-13, Bacosa+ 2015)
        'K_D_UV':           0.10,   # m^-1, UV attenuation (clear GoM)
        'UV_DEPTH_CUTOFF':  50.0,   # m, depth below which UV is negligible
        'OP_DISSOLVED_FRACTION': 0.5,   # dissolved vs oil-phase OP partition
        'K_OP_REMOVAL':     0.116,  # d^-1, OP removal rate (~6d t1/2 @20C)
    }

    # Default colors for plotting
    status_colors = {
        'initial': 'green',
        'active': 'blue',
        'missing_data': 'gray',
        'stranded': 'red',
        'evaporated': 'yellow',
        'dispersed': 'magenta',
        'dissolved': 'olive',
        'released_gas': 'orange',
        'dissolved_subsea': 'darkviolet',
        'biodegraded': 'purple'
    }

    duplicate_oils = ['ALVHEIM BLEND, STATOIL', 'DRAUGEN, STATOIL',
                      'EKOFISK BLEND 2000', 'EKOFISK BLEND, STATOIL',
                      'EKOFISK, CITGO', 'EKOFISK, EXXON', 'EKOFISK, PHILLIPS',
                      'EKOFISK, STATOIL', 'ELDFISK', 'ELDFISK B',
                      'GLITNE, STATOIL', 'GOLIAT BLEND, STATOIL',
                      'GRANE BLEND, STATOIL', 'GUDRUN BLEND, STATOIL',
                      'GULLFAKS A, STATOIL', 'GULLFAKS C, STATOIL',
                      'GULLFAKS, SHELL OIL', 'GULLFAKS SOR',
                      'GULLFAKS, STATOIL', 'HEIDRUN, STATOIL',
                      'NJORD, STATOIL', 'NORNE, STATOIL',
                      'OSEBERG BLEND, STATOIL', 'OSEBERG EXXON',
                      'OSEBERG, PHILLIPS', 'OSEBERG, SHELL OIL',
                      'SLEIPNER CONDENSATE, STATOIL',
                      'STATFJORD BLEND, STATOIL', 'VARG, STATOIL']

    # Workaround as ADIOS oil library uses max water fraction of 0.9 for all crude oils.
    # OD 1.14.9: renamed from max_water_fraction to avoid shadowing the instance-level
    # dict that emulsification_noaa() sets via prepare_run() (JSON-based, different structure).
    _cicoil_max_water_fraction_override = {
        'MARINE GAS OIL 500 ppm S 2017': 0.1,
        'FENJA (PIL) 2015': .75
    }

    dispersants = {'C9500': {'A': 17.,
                             'B': 0.24,
                             't12': 0.11, },
                   'OSR-52': {'A': 17.,
                              'B': 0.28,
                              't12': 0.24, },
                   'Dasic-NS': {'A': 17.,
                                'B': 0.28,
                                't12': 0.24, },
                   'Average': {'A': 17.,
                               'B': 0.26,
                               't12': 0.17, },
                   }


    def __init__(self, weathering_model='noaa', *args, **kwargs):


        # Calling general constructor of parent class
        super(OpenOil, self).__init__(*args, **kwargs)

        self.oil_weathering_model = weathering_model

        if self.oil_weathering_model == 'noaa' or self.oil_weathering_model == 'cicese':
            self.oiltypes = adios.get_oil_names(
                location=kwargs.get('location', None))

            # Update config with oiltypes
            # OD 1.14.9: oil_name_alias removed from adios module
            if hasattr(adios, 'oil_name_alias'):
                self.oiltypes.extend(adios.oil_name_alias.keys())

            # Sort alphabetically, but put GENERIC oils first
            generic_oiltypes = [o for o in self.oiltypes if o[0:7] == 'GENERIC']
            other_oiltypes = [o for o in self.oiltypes if o[0:7] != 'GENERIC']
            self.oiltypes = sorted([o for o in generic_oiltypes]) + sorted([o for o in other_oiltypes])
            self.oiltypes = [ot for ot in self.oiltypes if ot not in self.duplicate_oils]
        else:
            raise ValueError('Weathering model unknown: ' + weathering_model)

        kwargs.pop('location', None)

        if 'released_gas' not in self.status_categories:
            self.status_categories.append('released_gas')
        if 'dissolved_subsea' not in self.status_categories:
            self.status_categories.append('dissolved_subsea')

        self._add_config({
            'seed:m3_per_hour': {
                'type': 'float',
                'default': 1,
                'min': 0,
                'max': 1e10,
                'units': 'm3 per hour',
                'description':
                    'The amount (volume) of oil released per hour (or total amount if release is instantaneous)',
                'level': CONFIG_LEVEL_ESSENTIAL
            },
            'seed:droplet_diameter_min_subsea': {
                'type': 'float',
                'default': 0.0005,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description':
                    'The minimum diameter of oil droplet for a subsea release.',
                'level': CONFIG_LEVEL_BASIC
            },
            'seed:droplet_diameter_max_subsea': {
                'type': 'float',
                'default': 0.005,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description':
                    'The maximum diameter of oil droplet for a subsea release.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:dispersion': {
                'type': 'bool',
                'default': False,
                'description':
                    'Oil is removed from simulation (dispersed), if entrained as very small droplets.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:chemical_dispersion': {
                'type':
                    'enum',
                'enum': ['No_treatment', 'Surface', 'SSDI'],
                'default':
                    'No_treatment',
                'level':
                    CONFIG_LEVEL_ADVANCED,
                'description':
                    'Chemical dispersion method to be used.'
            },
            'processes:evaporation': {
                'type': 'bool',
                'default': True,
                'description': 'Surface oil is evaporated.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:emulsification': {
                'type': 'bool',
                'default': True,
                'description':
                    'Surface oil is emulsified, i.e. water droplets are mixed into oil due to wave mixing, with resulting increas of viscosity.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:biodegradation': {
                'type': 'bool',
                'default': False,
                'description': (
                    'Oil-phase pseudo-components are biodegraded by marine '
                    'bacteria. Per-group first-order rate constants for the '
                    '8 aliphatic and 8 aromatic groups plus residue, scaled '
                    'with temperature via an Adcroft et al. (2010) Q10=3 '
                    'law (T_REF=20°C).'),
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:update_oilfilm_thickness': {
                'type': 'bool',
                'default': False,
                'description':
                    'Oil film thickness is calculated at each time step. The alternative is that oil film thickness is kept constant with value provided at seeding.',
                'level': CONFIG_LEVEL_ADVANCED
            },
            'processes:handle_released_gas': {
                'type': 'bool',
                'default': True,
                'description': 'Handles gas bubbles.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:subsea_dissolution': {
                'type': 'bool',
                'default': True,
                'description': 'Dissolved components from submerged bubbles/droplets.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:spreading': {
                'type': 'bool',
                'default': True,
                'description': 'Surface oil is spreading.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:surface_dissolution': {
                'type': 'bool',
                'default': False,
                'description': 'Surface oil is dissolved.',
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:photooxidation': {
                'type': 'bool',
                'default': False,
                'description': (
                    'Aromatic pseudo-components are photooxidized by solar UV. '
                    'Uses synthetic diurnal UV (Bay of Campeche calibrated) and '
                    'wave-submergence correction (Delvigne & Sweeney 1988). '
                    'Rate constants increase from G1 (benzene-like) to G8 (pyrene-like).'),
                'level': CONFIG_LEVEL_BASIC
            },
            'processes:julia_weathering': {
                'type': 'bool',
                'default': False,
                'description': (
                    'Offload evaporation, emulsification, biodegradation and '
                    'photooxidation to CICOILPhysics.jl (Julia). Automatically '
                    'uses GPU via CUDA.jl when available, otherwise CPU. '
                    'Requires juliacall Python package and Julia >= 1.9.'),
                'level': CONFIG_LEVEL_ADVANCED
            },
            'julia:project_path': {
                'type': 'str',
                'default': '',
                'min_length': 0,
                'max_length': 1024,
                'description': (
                    'Path to CICOILPhysics.jl directory containing Project.toml. '
                    'Required when processes:julia_weathering is True.'),
                'level': CONFIG_LEVEL_ADVANCED
            },
            'wave_entrainment:droplet_size_distribution': {
                'type':
                    'enum',
                'enum': ['Johansen et al. (2015)', 'Li et al. (2017)'],
                'default':
                    'Johansen et al. (2015)',
                'level':
                    CONFIG_LEVEL_ADVANCED,
                'description':
                    'Algorithm to be used for calculating oil droplet size spectrum after entrainment by breaking waves.'
            },
            'wave_entrainment:entrainment_rate': {
                'type':
                    'enum',
                'enum': ['Li et al. (2017)'],
                'default':
                    'Li et al. (2017)',
                'level':
                    CONFIG_LEVEL_ADVANCED,
                'description':
                    'Algorithm to be used for calculating the entrainment rate of oil due to wave breaking.'
            },
            'seed:oil_type': {
                'type':
                    'enum',
                'enum':
                    self.oiltypes,
                'default':
                    self.oiltypes[0],
                'level':
                    CONFIG_LEVEL_ESSENTIAL,
                'description':
                    'Oil type to be used for the simulation, from the NOAA ADIOS database.'
            },
        })

        # OD 1.14.9: OpenOil.__init__ sets this to None; we skip OpenOil.__init__,
        # so must initialize here. prepare_run() populates it from JSON for NOAA mode.
        self.max_water_fraction = None

        self._set_config_default('drift:vertical_advection', False)
        self._set_config_default('drift:vertical_mixing', True)
        self._set_config_default('drift:current_uncertainty', 0.05)
        self._set_config_default('drift:wind_uncertainty', 0.5)
        # OD 1.14.9: replaces class-level required_profiles_z_range = [-20, 0]
        self._set_config_default('drift:profiles_depth', 20)  # OD 1.14.9: plural key

        # ── OD 1.14.9 new config keys — HC-2: biodegradation dispatch ─────────
        self._add_config({
            'biodegradation:method': {
                'type': 'enum',
                'enum': ['Adcroft', 'half_time'],
                'default': 'Adcroft',
                'level': CONFIG_LEVEL_ADVANCED,
                'description': (
                    'Biodegradation temperature-scaling law. CICOIL always '
                    'uses "Adcroft" (Q10=3, T_REF=20°C), applied per-component '
                    'in biodegradation_cicese() rather than the upstream bulk '
                    'biodegradation_adcroft().'),
            },
            # HC-3: subsea droplet size distribution (new in OD 1.14.9)
            'seed:droplet_size_distribution': {
                'type': 'enum',
                'enum': ['uniform', 'normal', 'lognormal'],
                'default': 'uniform',
                'level': CONFIG_LEVEL_ADVANCED,
                'description': 'Droplet size distribution for subsea release.'
            },
            'seed:droplet_diameter_mu': {
                'type': 'float',
                'default': 0.001,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description': 'Mean droplet diameter for normal/lognormal distributions.',
                'level': CONFIG_LEVEL_BASIC
            },
            'seed:droplet_diameter_sigma': {
                'type': 'float',
                'default': 0.0005,
                'min': 1e-8,
                'max': 1,
                'units': 'meters',
                'description': 'Std dev of droplet diameter for normal/lognormal distributions.',
                'level': CONFIG_LEVEL_BASIC
            },
        })
        # Force Adcroft biodegradation — HC-2
        self._set_config_default('biodegradation:method', 'Adcroft')


    def run(self, record_components=False, record_environment=True, *args, **kwargs):
        """Run the simulation.

        :param record_components: if True, export the 51 per-pseudo-component
            arrays (comp_*, comp_*_evaporated, comp_*_dissolved) in addition
            to the bulk mass-balance variables. These dominate output file
            size for large ensembles (17 components x 3 variants x N
            particles x N output steps) and default to off.
        :param record_environment: if False (and record_components is also
            False, and no explicit `export_variables` is given), drop the 20
            `required_variables` environment fields (current, wind, waves,
            SST, salinity, etc.) sampled per element from the output. These
            are useful for QC but for large ensembles can rival the
            bulk mass-balance variables in size; the same values can be
            recovered from the forcing data given each element's
            position/time. Has no effect if record_components=True or an
            explicit `export_variables` is passed.
        """
        if 'export_variables' not in kwargs and not record_components:
            export_variables = [key for key in self.ElementType.variables.keys()
                                if key[0:5]!='comp_']
            if record_environment:
                export_variables += self.required_variables.keys()
            kwargs['export_variables'] = export_variables
        super(OpenCiceseOil, self).run(*args, **kwargs)


    def diffusion_coeff_in_air(self, atm_pressure=101325.):
        """ Vapour diffusion coefficients in Air according to Hilal and Karickhoff (2003), Kim and Monroe (2014)

        :param surface_indices: Indices of elementss on the sea surface
        :param atm_pressure: In Pa
        """
        # Air properties
        ek_air = 97.    # e/k number for air
        r_air = 3.617    # collision radius of air
        molar_mass_air = 28.965

        # seawater temperature at the start and at the location of the first element
        # OD 1.14.9: get_environment moved to self.env; returns (recarray, profiles, missing)
        _env, _, _ = self.env.get_environment(
            ['sea_water_temperature'],
            self.start_time,
            np.array([self.elements_scheduled.lon[0]]),
            np.array([self.elements_scheduled.lat[0]]),
            np.array([0.0]),
        )
        seawater_temp = float(_env['sea_water_temperature'][0])

        boiling_temp = self.fluid_properties.component_boil_temp
        molar_mass = self.fluid_properties.component_molar_mass
        molar_volume = molar_mass / self.fluid_properties.component_density

        ek_vapour = 1.15 * boiling_temp
        ek_air_vapour = np.sqrt(ek_air * ek_vapour)
        r_vapour = 1.18 * np.power(molar_volume * 10.**6, 1./3.)    # collision radius of vapour
        r_air_vapour = (r_air + r_vapour) / 2.

        reduced_temp = seawater_temp / ek_air_vapour

        # Reduced Lenard-Jones integrals from Kim and Monroe (2014)
        integrals = {}
        b_factors = np.array([2.6431984, 6.0432255e-3, -1.5158773e-1, 5.4237938e-2, -9.0468682e-3, 6.1742007e-4], np.float32)
        c_factors = np.array([1.6690746, -6.914589e-1, 1.5502132e-1, -2.0642189e-2, 1.5402077e-3, -4.9729535e-5], np.float32)
        for idx in range(6):
            integrals[idx] = b_factors[idx] / np.power(reduced_temp, idx + 1.) + c_factors[idx] * np.power(np.log(reduced_temp), idx + 1.)

        omega11 = integrals[0] + integrals[1] + integrals[2] + integrals[3] + integrals[4] + integrals[5] - 1.1036729
        Mab = 2. / (1./ molar_mass_air + 1./ (molar_mass * 1000))
        diffusion_coeff = (3.03 - 0.98 / np.sqrt(Mab)) * 10.**-3 * np.power(seawater_temp, 3. / 2.) / \
                          (atm_pressure * 10.**-5 * np.sqrt(Mab) * np.power(r_air_vapour, 2) * omega11) * 10.**-4

        #######################

        air_kin_viscosity = (0.01827 * (291.15 + 120.) / (seawater_temp + 120.) *
                             np.power(seawater_temp / 291.15, 3. / 2.)) / (3.485 * atm_pressure / seawater_temp)

        schmidt_number = air_kin_viscosity / diffusion_coeff
        self.feemnumber_reynolds = np.power(schmidt_number, 2. / 3.) / np.sqrt(1. + np.power(0.0468 / schmidt_number, 2. / 3.))

        self.air_kin_viscosity = air_kin_viscosity
        #######################

        return diffusion_coeff


    def transfer_coeff_air(self, surface_indices):

        surf_no = len(surface_indices)
        components_no = len(self.fluid_properties.oil_mass_frac)

        wind_velocity = np.sqrt(self.environment.x_wind[surface_indices]**2. +
                                self.environment.y_wind[surface_indices]**2.)
        spillet_area = self.elements.spillet_area[surface_indices]
        spillet_diameter = np.sqrt(spillet_area / np.pi) * 2.

        # update diffusion coefficients
        diffusion_coeff = self.fluid_properties.component_diffusion * (np.ones(surf_no).reshape((surf_no, -1)))
        feemnumber_reynolds = self.feemnumber_reynolds * (np.ones(surf_no).reshape((surf_no, -1)))

        reynolds_number = wind_velocity * spillet_diameter / self.air_kin_viscosity
        reynolds_number_matrix = np.ones(components_no) * (reynolds_number.reshape((surf_no, -1)))

        feem_number = reynolds_number_matrix * feemnumber_reynolds

        sherwood_number = 0.6774 * np.sqrt(feem_number) * \
                          np.sqrt(1. + np.power(feem_number / 12500., 3./5.) / \
                                  np.power(1. + np.power(300000. / feem_number, 7./2.), 2./5.))

        spillet_diameter_matrix = np.ones(components_no) * (spillet_diameter.reshape((surf_no, -1)))
        xfer_coeff = sherwood_number * diffusion_coeff / spillet_diameter_matrix
        return xfer_coeff


    def spreading_cicese(self, gravity=9.81, drag_coef=0.003, film_limit=0.000005):
        """ Calculates the gravity-viscous spread according to Nihoul (1984)
            :param gravity: gravity acceleration (m / s2)
            :param drag_coef: drag coefficient (kg / m2 s)
        """

        logging.debug('   Calculating: spreading')

        # elements at surface
        spread = np.where(self.elements.z == 0)[0]
        if len(spread) == 0: return

        # Previous step spillet radius
        spillet_thickness_previous = self.elements.spillet_thickness[spread]

        # T0 = self.environment.sea_water_temperature[spread]
        # S0 = self.environment.sea_water_salinity[spread]
        rho_oil = self.oil_density(spread)
        rho_water = 1025.     # self.sea_water_density(T=T0, S=S0)

        # Current spillet volume
        spillet_volume = 1 / rho_oil    # pseudo volume intedependant of spillet mass selection
        spillet_area_previous = spillet_volume / spillet_thickness_previous
        radius_pr = np.sqrt(spillet_area_previous) / np.sqrt(np.pi)
        element_volume = self.elements.mass_oil[spread] / rho_oil

        g = gravity  # gravity acceleration (m / s2)
        g_r = (rho_water - rho_oil) * g / rho_water  # Reduced gravity

        k = drag_coef  # drag coefficient (kg / m2 s)
        alpha = g_r * rho_oil / k

        # Calculate pseudo-time that corresponds to the radius of the previous
        # step, since spillet volume has changed due to weathering
        pseudo_time = 2. * np.pi * np.power(radius_pr, 6.) / (27. * np.power(spillet_volume, 2.) * alpha)
        spread_time = pseudo_time + self.time_step.total_seconds()
        radius = np.power((27. * np.power(spillet_volume, 2.) * alpha * spread_time / (2. * np.pi)), (1. / 6.))
        spillet_area_new = np.pi * np.power(radius, 2.)
        spillet_thickness_new = spillet_volume / spillet_area_new
        spillet_thickness_new[spillet_thickness_new < film_limit] = film_limit
        self.elements.spillet_thickness[spread] = spillet_thickness_new
        self.elements.spillet_area[spread] = element_volume / spillet_thickness_new


    def oil_density(self, indices):

        surfaceID = indices  # BC-2: OD 1.14.9 uses 0-based element IDs (was -1)
        mass_component = self.cicese_mass_balance['mass_components'][surfaceID, :]
        surf_no = len(indices)
        density_component = self.fluid_properties.component_density * (np.ones(surf_no).reshape((surf_no, -1)))
        volume_component = mass_component / density_component
        return np.sum(mass_component, 1) / np.sum(volume_component, 1)


    def component_molfrac(self, surface):

        surf_no = len(surface)
        components_no = len(self.fluid_properties.oil_mass_frac)
        Mw_component = self.fluid_properties.component_molar_mass * (np.ones(surf_no).reshape((surf_no, -1)))

        surfaceID = self.elements.ID[surface]  # BC-2: 0-based in OD 1.14.9
        mass_component = self.cicese_mass_balance['mass_components'][surfaceID, :]

        moles_component = mass_component / Mw_component

        moles_sum = np.ones(components_no) * (np.sum(moles_component, 1).reshape((surf_no, -1)))
        moles_frac = moles_component / moles_sum

        #TODO Move component Mw update to an appropriate location
        self.elements.oil_molar_mass[surface] = np.sum(moles_frac * Mw_component, 1)

        return moles_frac


    def evaporation_cicese(self, atm_pressure=101325.):
        #############################################
        # Evaporation, for elements at surface only
        #############################################
        logging.debug('    Calculating evaporation - CICESE')
        surface = np.where(self.elements.z == 0)[0]  # of active elements

        #     update_prop = np.where(self.elements.z >= 0 and self.elements.update_properties)[0]

        if len(surface) == 0:
            logging.debug('All elements submerged, no evaporation')
            return
        if self.elements.age_seconds[surface].min() > 3600 * 24 * 120:
            logging.debug('All surface oil elements older than 240 hours, ' +
                          'skipping further evaporation.')
            return

        surfaceID = self.elements.ID[surface]  # BC-2: 0-based in OD 1.14.9
        surf_no = len(surface)
        components_no = len(self.fluid_properties.oil_mass_frac)
        moles_frac = self.component_molfrac(surface)

        spillet_Mw = np.atleast_1d(self.elements.oil_molar_mass[surface])
        spillet_density = np.atleast_1d(self.elements.oil_density[surface])
        seawater_temp = np.atleast_1d(self.environment.sea_water_temperature[surface])  # 299.85 K GM average
        slick_thickness = np.atleast_1d(self.elements.spillet_thickness[surface])

        molar_mass = np.ones(components_no) * (spillet_Mw.reshape((surf_no, -1)))
        density = np.ones(components_no) * (spillet_density.reshape((surf_no, -1)))
        pseudo_temp = np.ones(components_no) * (seawater_temp.reshape((surf_no, -1)))
        pseudo_thickness = np.ones(components_no) * (slick_thickness.reshape((surf_no, -1)))

        transfer_coef = self.transfer_coeff_air(surface)
        pseudo_Tb = self.fluid_properties.component_boil_temp * (np.ones(surf_no).reshape((surf_no, -1)))

        Cox_A0 = self.fluid_properties.Cox_A0 * (np.ones(surf_no).reshape((surf_no, -1)))
        Cox_A1 = self.fluid_properties.Cox_A1 * (np.ones(surf_no).reshape((surf_no, -1)))
        Cox_A2 = self.fluid_properties.Cox_A2 * (np.ones(surf_no).reshape((surf_no, -1)))

        Rgas_constant = 8.31446
        pseudo_Pv = atm_pressure * np.exp((1. - pseudo_Tb / pseudo_temp) * \
                                          np.exp(Cox_A0 + Cox_A1 * pseudo_temp +
                                                 Cox_A2 * np.power(pseudo_temp, 2.)))

        # Store evaporated mass at beginning of timestep
        mass_evaporated_previous = self.elements.mass_evaporated[surface]

        # Calculate the fraction of each pseudo-component that has evaporated
        pseudo_fraction = transfer_coef / pseudo_thickness * molar_mass / density * \
                          pseudo_Pv * moles_frac / (Rgas_constant * pseudo_temp) * self.time_step.total_seconds()
        pseudo_fraction = np.where(pseudo_fraction > 1, 1, pseudo_fraction)
        
        mass_evaporated = self.cicese_mass_balance['mass_components'][surfaceID, :] * pseudo_fraction

        self.cicese_mass_balance['mass_evaporated'][surfaceID, :] += mass_evaporated

        mass_remain = self.cicese_mass_balance['mass_components'][surfaceID, :] - mass_evaporated
        mass_remain[mass_remain < 0] = 0
        self.cicese_mass_balance['mass_components'][surfaceID, :] = mass_remain

        element_mass = np.sum(mass_remain, 1)
        element_evaporated = mass_evaporated_previous + np.sum(mass_evaporated, 1)

        self.elements.mass_evaporated[surface] = element_evaporated
        self.elements.mass_oil[surface] = element_mass

        for i, cm in enumerate(self.fluid_properties.composition[9:]):  # update component masses with remaining masses
            array = getattr(self.elements, 'comp_'+cm)
            array[surface] = mass_remain[:, i]
            setattr(self.elements, 'comp_'+cm, self.elements.variables['comp_'+cm]['dtype'](array))

        # update component evaporated masses with new values
        for i, ce in enumerate(self.fluid_properties.composition[9:]):
            array = getattr(self.elements, 'comp_' + ce + '_evaporated')
            array[surface] = mass_remain[:, i]
            setattr(self.elements, 'comp_' + ce + '_evaporated',
                    self.elements.variables['comp_' + ce + '_evaporated']['dtype'](array))

        # TODO: FIX deactivation
        # evaporated_indices = np.where(element_mass <=0.001)[0]
        # self.deactivate_elements(evaporated_indices, reason='evaporated')

        #self.elements.mass_oil[evaporated_indices] = 0
        self.elements.fraction_evaporated[surface] = element_evaporated / \
                                                    (element_evaporated + element_mass +
                                                     self.elements.mass_dispersed[surface] +
                                                     self.elements.mass_dissolved_surface[surface] +
                                                     self.elements.mass_biodegraded_from_oil[surface] +
                                                     self.elements.mass_biodegraded_from_water[surface])


    def oil_weathering(self):
        if self.time_step.days < 0:
            logging.debug('Skipping oil weathering for backwards run')
            return
        self.timer_start('main loop:updating elements:oil weathering')
        if self.oil_weathering_model == 'cicese':
            self.oil_weathering_cicese()
        else:
            super(OpenCiceseOil, self).oil_weathering()
        self.timer_end('main loop:updating elements:oil weathering')


    def prepare_run(self):

        if self.oil_weathering_model == 'cicese':
            self.fluid_properties.get_live_composition(GOR=0)
            self.fluid_properties.component_diffusion = self.diffusion_coeff_in_air()

            # Populate with seeded mass spread on oiltype.mass_fraction
            self.cicese_mass_balance = {}

            mass_oil = np.atleast_1d(self.elements_scheduled.mass_oil)
            if len(mass_oil) == 1:
                mass_oil = mass_oil*np.ones(self.num_elements_total())
            self.cicese_mass_balance['mass_components'] = self.fluid_properties.oil_mass_frac * \
                                                          (mass_oil.reshape((self.num_elements_total(), -1)))
            self.cicese_mass_balance['mass_evaporated'] = self.cicese_mass_balance['mass_components'] * 0
            self.cicese_mass_balance['mass_dissolved'] = self.cicese_mass_balance['mass_components'] * 0
            self.cicese_mass_balance['mass_biodegraded_from_oil'] = self.cicese_mass_balance['mass_components'] * 0
            # mass_dissolved is populated by handle_subsea_dissolution() for
            # is_dissolved==1 particles, and decayed via K_BIODEG_WATER into
            # mass_biodegraded_from_water by biodegradation_cicese().
            self.cicese_mass_balance['mass_biodegraded_from_water'] = self.cicese_mass_balance['mass_components'] * 0
            # Photooxidation per-component tracking (aromatic groups only, indices 8-15)
            self.cicese_mass_balance['mass_photooxidized'] = self.cicese_mass_balance['mass_components'] * 0
            # Oxygenated photoproduct (OP) pools (audit #4), aromatic groups
            # only (indices 8-15): mass_op_dissolved is the dissolved OP
            # pool fed by photooxidation_cicese(); mass_op_degraded is its
            # decay product (K_OP_REMOVAL).
            self.cicese_mass_balance['mass_op_dissolved'] = self.cicese_mass_balance['mass_components'] * 0
            self.cicese_mass_balance['mass_op_degraded'] = self.cicese_mass_balance['mass_components'] * 0

            self.elements_scheduled.spillet_area = mass_oil / self.elements_scheduled.oil_density / 0.01
            self.elements_scheduled.water_fraction = 0.

            n_elem = self.num_elements_total()
            n_mb = self.cicese_mass_balance['mass_components'].shape[0]
            assert n_mb == n_elem, (
                f'cicese_mass_balance shape[0]={n_mb} != num_elements={n_elem}')

            dt_s = self.time_step.total_seconds() if self.time_step else 3600.0
            v_max = 1.5  # m/s, typical GoM current upper bound
            dx = 1000.0  # m, conservative grid scale estimate
            for rdr in getattr(self.env, 'readers', getattr(self, 'readers', {})).values():
                if hasattr(rdr, 'delta_x') and rdr.delta_x is not None:
                    dx = min(dx, rdr.delta_x)
                    break
            courant = dt_s * v_max / dx
            if courant > 1.0:
                logger.warning('Courant number estimate %.2f > 1 '
                               '(dt=%.0fs, v_max=%.1fm/s, dx=%.0fm)',
                               courant, dt_s, v_max, dx)
            else:
                logger.info('Courant number estimate %.2f '
                            '(dt=%.0fs, v_max=%.1fm/s, dx=%.0fm)',
                            courant, dt_s, v_max, dx)
        else:
            super(OpenCiceseOil, self).prepare_run()

    # ── Julia/CICOILPhysics backend ─────────────────────────────────────────
    _jl = None          # juliacall Main — initialised once
    _jl_state = None    # SimState opaque handle

    def _init_julia_backend(self):
        """Lazy-init Julia + CICOILPhysics.jl + SimState on first call."""
        if OpenCiceseOil._jl is not None:
            return
        julia_project = self.get_config('julia:project_path')
        if not julia_project:
            raise ValueError(
                'processes:julia_weathering is True but julia:project_path '
                'is not set. Point it to the CICOILPhysics.jl directory.')
        import os
        julia_project = os.path.abspath(julia_project)

        # Use JULIA_PROJECT env var if set (e.g. envs/gpu or envs/cpu),
        # otherwise auto-detect: prefer envs/gpu if it exists, then envs/cpu,
        # then fall back to the bare project (requires Pkg.instantiate).
        julia_env = os.environ.get('JULIA_PROJECT', '')
        if not julia_env:
            for env in ['envs/gpu', 'envs/cpu']:
                candidate = os.path.join(julia_project, env)
                if os.path.isfile(os.path.join(candidate, 'Manifest.toml')):
                    julia_env = candidate
                    break
            if not julia_env:
                julia_env = julia_project
        julia_env = os.path.abspath(julia_env)
        os.environ["JULIA_PROJECT"] = julia_env

        from juliacall import Main as jl
        jl.seval(f'import Pkg; Pkg.activate("{julia_env}", io=devnull)')
        jl.seval('using CICOILPhysics')
        OpenCiceseOil._jl = jl
        logger.info('CICOILPhysics.jl loaded (env: %s) — backend: %s',
                     julia_env,
                     str(jl.CICOILPhysics.backend_name(jl.CICOILPhysics.backend())))

    def _init_julia_state(self):
        """Allocate Julia SimState from the current oil type."""
        if OpenCiceseOil._jl_state is not None:
            return
        jl = OpenCiceseOil._jl
        fp = self.fluid_properties
        OpenCiceseOil._jl_state = jl.CICOILPhysics.py_init_state(
            N_particles=int(self.num_elements_total()),
            oil_Mw=np.array(fp.comp_Mw[-17:], dtype=np.float32),
            oil_Tb=np.array(fp.component_boil_temp[-17:], dtype=np.float32),
            oil_mass_frac=np.array(fp.oil_mass_frac[-17:], dtype=np.float32),
            initial_mass_kg=float(np.sum(
                self.cicese_mass_balance['mass_components'][0])),
        )
        logger.info('Julia SimState allocated for %d particles (%s)',
                     self.num_elements_total(),
                     str(jl.CICOILPhysics.backend_name(jl.CICOILPhysics.backend())))

    def _julia_weathering_step(self):
        """Run the 4 weathering kernels via CICOILPhysics.jl."""
        self._init_julia_backend()
        self._init_julia_state()
        jl = OpenCiceseOil._jl
        state = OpenCiceseOil._jl_state

        mb = self.cicese_mass_balance
        el = self.elements
        N = len(el.z)
        dt = float(self.time_step.total_seconds())
        uv_irradiance = float(self._compute_uv_irradiance())

        # Mutable float32 buffers — Julia writes results into these
        buf_mass_oil       = np.array(el.mass_oil, dtype=np.float32)
        buf_mass_evap      = np.array(el.mass_evaporated, dtype=np.float32)
        buf_mass_biodeg    = np.array(el.mass_biodegraded, dtype=np.float32)
        buf_mass_biodeg_oil= np.array(el.mass_biodegraded_from_oil, dtype=np.float32)
        buf_mass_biodeg_w  = np.array(el.mass_biodegraded_from_water, dtype=np.float32)
        buf_mass_photoox   = np.array(el.mass_photooxidized, dtype=np.float32)
        buf_mass_op_diss   = np.array(el.mass_op_dissolved, dtype=np.float32)
        buf_mass_op_deg    = np.array(el.mass_op_degraded, dtype=np.float32)
        buf_water_frac     = np.array(el.water_fraction, dtype=np.float32)
        buf_kvisc          = np.zeros(N, dtype=np.float32)
        buf_mc             = np.array(mb['mass_components'], dtype=np.float32)
        buf_me             = np.array(mb['mass_evaporated'], dtype=np.float32)
        buf_mb             = np.array(mb['mass_biodegraded_from_oil'], dtype=np.float32)
        buf_mp             = np.array(mb['mass_photooxidized'], dtype=np.float32)
        buf_mod            = np.array(mb['mass_op_dissolved'], dtype=np.float32)
        buf_mog            = np.array(mb['mass_op_degraded'], dtype=np.float32)

        jl.CICOILPhysics.py_update_weathering(
            state,
            np.array(el.z, dtype=np.float32),
            np.array(el.diameter, dtype=np.float32),
            np.array(el.age_seconds, dtype=np.float32),
            np.array(el.is_dissolved, dtype=np.uint8),
            np.array(el.oil_molar_mass, dtype=np.float32),
            np.array(el.oil_density, dtype=np.float32),
            np.array(el.spillet_thickness, dtype=np.float32),
            buf_water_frac,
            buf_kvisc,
            np.array(el.max_water, dtype=np.float32),
            buf_mass_oil,
            buf_mass_evap,
            np.array(el.mass_dispersed, dtype=np.float32),
            buf_mass_biodeg,
            buf_mass_biodeg_oil,
            buf_mass_biodeg_w,
            buf_mass_photoox,
            buf_mass_op_diss,
            buf_mass_op_deg,
            buf_mc, buf_me, buf_mb, buf_mp, buf_mod, buf_mog,
            np.array(mb['mass_dissolved'], dtype=np.float32),
            np.array(self.environment.sea_water_temperature, dtype=np.float32),
            np.array(self.environment.x_wind, dtype=np.float32),
            np.array(self.environment.y_wind, dtype=np.float32),
            np.array(self.environment.x_sea_water_velocity, dtype=np.float32),
            np.array(self.environment.y_sea_water_velocity, dtype=np.float32),
            np.array(getattr(self.environment,
                    'sea_surface_wave_significant_height',
                    np.zeros(N, dtype=np.float32)), dtype=np.float32),
            np.array(getattr(self.environment,
                    'sea_surface_wave_period_at_variance_spectral_density_maximum',
                    np.zeros(N, dtype=np.float32)), dtype=np.float32),
            uv_irradiance,
            bool(self.get_config('processes:evaporation')),
            bool(self.get_config('processes:emulsification')),
            bool(self.get_config('processes:biodegradation')),
            bool(self.get_config('processes:photooxidation')),
            dt,
        )

        # Write results back into OpenDrift element arrays
        el.mass_oil[:]                    = buf_mass_oil
        el.mass_evaporated[:]             = buf_mass_evap
        el.mass_biodegraded[:]            = buf_mass_biodeg
        el.mass_biodegraded_from_oil[:]   = buf_mass_biodeg_oil
        el.mass_biodegraded_from_water[:] = buf_mass_biodeg_w
        el.mass_photooxidized[:]          = buf_mass_photoox
        el.mass_op_dissolved[:]           = buf_mass_op_diss
        el.mass_op_degraded[:]            = buf_mass_op_deg
        el.water_fraction[:]              = buf_water_frac
        el.viscosity[:]                   = buf_kvisc * 1e-6  # cSt → m²/s
        mb['mass_components'][:]          = buf_mc
        mb['mass_evaporated'][:]          = buf_me
        mb['mass_biodegraded_from_oil'][:]= buf_mb
        mb['mass_photooxidized'][:]       = buf_mp
        mb['mass_op_dissolved'][:]        = buf_mod
        mb['mass_op_degraded'][:]         = buf_mog

    def oil_weathering_cicese(self):
        '''Oil weathering scheme adopted from NOAA PyGNOME model:
        https://github.com/NOAA-ORR-ERD/PyGnome
        '''
        logging.debug('CICESE oil weathering')
        # C to K
        self.environment.sea_water_temperature[
            self.environment.sea_water_temperature < 100] += 273.15

        #########################################################
        # Update density and viscosity according to temperature
        #########################################################
        self.timer_start('main loop:updating elements:oil weathering:updating viscosities')
        oil_viscosity = self.oiltype.kvis_at_temp(
            self.environment.sea_water_temperature)
        self.timer_end('main loop:updating elements:oil weathering:updating viscosities')
        self.timer_start('main loop:updating elements:oil weathering:updating densities')
        active_elements = np.where(self.elements.z<=0)[0]
        oil_density = self.oil_density(active_elements)    # New oil density set by Cicese calculations
        self.elements.oil_density[active_elements] = oil_density
        self.timer_end('main loop:updating elements:oil weathering:updating densities')

        # Calculate emulsion density
        self.elements.density = (
            self.elements.water_fraction*self.sea_water_density() +
           (1 - self.elements.water_fraction) * self.elements.oil_density)

        # Calculate emulsion viscosity
        visc_f_ref = 0.84  # From PyGNOME
        visc_curvfit_param = 1.5e3 # units are sec^0.5 / m
        fw_d_fref = self.elements.water_fraction/visc_f_ref
        kv1 = np.sqrt(oil_viscosity)*visc_curvfit_param
        kv1[kv1<1] = 1
        kv1[kv1>10] = 10
        self.elements.viscosity = (
            oil_viscosity*np.exp(kv1*self.elements.fraction_evaporated)*   # New fraction evaporated method calculation
                (1 + (fw_d_fref / (1.187 - fw_d_fref))) ** 2.49)

        # Dry-oil (no-emulsion), evaporation-corrected kinematic viscosity [m^2/s],
        # used by emulsification_cicese() for the viscosity-stability scaling factor.
        self._oil_viscosity_dry = oil_viscosity * np.exp(kv1 * self.elements.fraction_evaporated)

        # ── Processes that always run in Python ────────────────────────────────
        if self.get_config('processes:spreading') is True:
            self.timer_start('main loop:updating elements:oil weathering:spreading')
            self.spreading_cicese()
            self.timer_end('main loop:updating elements:oil weathering:spreading')

        if self.get_config('processes:handle_released_gas') is True:
            self.timer_start('main loop:updating elements:oil weathering:gas_particles')
            self.handle_gas_particles()
            self.timer_end('main loop:updating elements:oil weathering:gas_particles')

        if self.get_config('processes:subsea_dissolution') is True:
            self.timer_start('main loop:updating elements:oil weathering:subsea_dissolution')
            self.handle_subsea_dissolution()
            self.timer_end('main loop:updating elements:oil weathering:subsea_dissolution')

        if self.get_config('processes:dispersion') is True:
            self.timer_start('main loop:updating elements:oil weathering:dispersion')
            self.dispersion_cicese()
            self.timer_end('main loop:updating elements:oil weathering:dispersion')

        if self.get_config('processes:surface_dissolution') is True:
            self.timer_start('main loop:updating elements:oil weathering:surface_dissolution')
            self.surface_dissolution()
            self.timer_end('main loop:updating elements:oil weathering:surface_dissolution')

        # ── Evaporation / emulsification / biodegradation / photooxidation ───
        # When julia_weathering is enabled, all 4 run in a single Julia call.
        # Otherwise, the original Python methods are called individually.
        if self.get_config('processes:julia_weathering') is True:
            self.timer_start('main loop:updating elements:oil weathering:julia_kernels')
            self._julia_weathering_step()
            self.timer_end('main loop:updating elements:oil weathering:julia_kernels')
        else:
            if self.get_config('processes:evaporation') is True:
                self.timer_start('main loop:updating elements:oil weathering:evaporation')
                self.evaporation_cicese()
                self.timer_end('main loop:updating elements:oil weathering:evaporation')

            if self.get_config('processes:emulsification') is True:
                self.timer_start('main loop:updating elements:oil weathering:emulsification')
                self.emulsification_cicese()
                self.timer_end('main loop:updating elements:oil weathering:emulsification')

            if self.get_config('processes:biodegradation') is True:
                self.timer_start('main loop:updating elements:oil weathering:biodegradation')
                self.biodegradation_cicese()
                self.timer_end('main loop:updating elements:oil weathering:biodegradation')

            if self.get_config('processes:photooxidation') is True:
                self.timer_start('main loop:updating elements:oil weathering:photooxidation')
                self.photooxidation_cicese()
                self.timer_end('main loop:updating elements:oil weathering:photooxidation')


    def handle_gas_particles(self):
        """for now gas bubble mass is recorded and gas bubbles are deactivated """
        bubbles = np.where(self.elements.is_bubble == 1)[0]
        if len(bubbles) > 0:
            self.deactivate_elements(bubbles, reason='released_gas')


    def handle_subsea_dissolution(self):
        """Dissolved-phase mass bookkeeping for is_dissolved==1 particles.

        On the seeding step, each element's per-component mass is
        transferred once from `cicese_mass_balance['mass_components']` into
        `cicese_mass_balance['mass_dissolved']` (zeroing mass_components),
        instead of deactivating the element immediately. This gives
        biodegradation_cicese() a persistent dissolved-mass pool to decay via
        K_BIODEG_WATER (mass_biodegraded_from_water). Once that pool has
        decayed below BIODEG_WATER_DEPLETION_FRAC of its initial mass
        (elements.mass_dissolved_subsea), the element is deactivated.
        """
        dissolved = np.where(self.elements.is_dissolved == 1)[0]
        if len(dissolved) == 0:
            return

        dissolvedID = self.elements.ID[dissolved]
        mass_dis = self.cicese_mass_balance['mass_dissolved'][dissolvedID, :]
        mass_dis_total = np.sum(mass_dis, axis=1)

        # One-time transfer of component mass into the dissolved pool, for
        # elements that haven't been transferred yet (mass_dissolved == 0).
        # elements.mass_dissolved_subsea is overwritten with the transferred
        # total, so it is a self-consistent "initial dissolved mass"
        # reference for the depletion check and fraction_biodegraded_from_water,
        # independent of whatever value it held at seeding.
        not_yet_transferred = mass_dis_total == 0
        if np.any(not_yet_transferred):
            idx = dissolved[not_yet_transferred]
            idxID = self.elements.ID[idx]
            transferred = self.cicese_mass_balance['mass_components'][idxID, :]
            self.cicese_mass_balance['mass_dissolved'][idxID, :] = transferred
            self.cicese_mass_balance['mass_components'][idxID, :] = 0.0
            self.elements.mass_dissolved_subsea[idx] = np.sum(transferred, axis=1)

        # Deactivate elements whose dissolved pool has mostly decayed away
        # (skip elements transferred this step, so decay gets at least one step).
        initial_dissolved = self.elements.mass_dissolved_subsea[dissolved]
        mask_init = initial_dissolved > 0
        depleted = np.zeros(len(dissolved), dtype=bool)
        depleted[mask_init] = (mass_dis_total[mask_init] <=
                                self.BIODEG_WATER_DEPLETION_FRAC * initial_dissolved[mask_init])
        to_deactivate = dissolved[depleted & ~not_yet_transferred]
        if len(to_deactivate) > 0:
            self.deactivate_elements(to_deactivate, reason='dissolved_biodegraded')


    def dispersion_cicese(self):
        logging.debug('    Calculating: dispersion - NOAA')
        # From NOAA PyGnome model:
        # https://github.com/NOAA-ORR-ERD/PyGnome/
        c_disp = np.power(self.wave_energy_dissipation(), 0.57) * \
            self.sea_surface_wave_breaking_fraction()
        # Roy's constant
        C_Roy = 2400.0 * np.exp(-73.682*np.sqrt(
            self.elements.viscosity/self.elements.density))
        v_entrain = 3.9E-8
        q_disp = C_Roy * c_disp * v_entrain / self.elements.density
        fraction_dispersed = (q_disp * self.time_step.total_seconds() *
                         self.elements.density)
        oil_mass_loss = fraction_dispersed*self.elements.mass_oil

        self.cicese_mass_balance['mass_components'][self.elements.ID, :] = \
            self.cicese_mass_balance['mass_components'][self.elements.ID, :]*(1-fraction_dispersed[:, np.newaxis])  # BC-2

        self.elements.mass_oil -= oil_mass_loss
        self.elements.mass_dispersed += oil_mass_loss

    def oil_wave_entrainment_rate(self):
        er = self.get_config('wave_entrainment:entrainment_rate')
        if er == 'Li et al. (2017)':
            if self.oil_weathering_model == 'cicese':
                entrainment_rate = oil_wave_entrainment_rate_li2017(
                    dynamic_viscosity=self.elements.viscosity *
                    self.elements.density,
                    oil_density=self.elements.density,
                    interfacial_tension=self.elements.IFT,
                    significant_wave_height=self.significant_wave_height(),
                    wave_breaking_fraction=self.sea_surface_wave_breaking_fraction(
                    ),
                    sea_water_density=self.sea_water_density())
            else:
                entrainment_rate = oil_wave_entrainment_rate_li2017(
                    dynamic_viscosity=self.elements.viscosity *
                                      self.elements.density,
                    oil_density=self.elements.density,
                    interfacial_tension=self.oil_water_interfacial_tension,
                    significant_wave_height=self.significant_wave_height(),
                    wave_breaking_fraction=self.sea_surface_wave_breaking_fraction(
                    ),
                    sea_water_density=self.sea_water_density())
            return entrainment_rate
        else:
            logger.error("no entrainment rate mechanism configured")
            raise Exception("no entrainment rate mechanism configured")

    def get_wave_breaking_droplet_diameter_liz2017(self):
        # Li,Zhengkai, M. Spaulding, D. French-McCay, D. Crowley, J.R. Payne: "Development of a unified oil droplet size distribution model
        # with application to surface breaking waves and subsea blowout releases considering dispersant effects" Mar. Pol. Bul.
        # DOI: 10.1016/j.marpolbul.2016.09.008
        # Should be prefered when the oil film thickness is unknown.
        if not hasattr(self, 'droplet_spectrum_pdf'):
            # Generate droplet spectrum as in Li (Zhengkai) et al. (2017)
            # Bounds are hardcoded to 1 micron and 3mm
            logger.debug('Generating wave breaking droplet size spectrum')
            self.droplet_spectrum_diameter = np.linspace(1e-6, 3e-3, 1000000)
            g = 9.81
            if self.oil_weathering_model == 'cicese':
                interfacial_tension = self.elements.IFT
            else:
                interfacial_tension = self.oil_water_interfacial_tension
            delta_rho = self.sea_water_density() - self.elements.density
            d_o = 4 * (interfacial_tension / (delta_rho * g))**0.5
            we = (self.sea_water_density() * g *
                  self.significant_wave_height() * d_o) / interfacial_tension
            oh = self.elements.viscosity * self.elements.density * (
                self.elements.density * interfacial_tension *
                d_o)**-0.5  # From kin. to dyn. viscosity by * density
            r = 1.791
            p = 0.460
            q = -0.518
            dV_50 = d_o * r * (
                1 + 10 * oh
            )**p * we**q  # median droplet diameter in volume distribution
            sd = 0.4  # log standard deviation in log10 units
            Sd = np.log(10) * sd  # log standard deviation in natural log units
            # TODO: calculation below with scalars, but we have arrays, with varying oil properties
            # treat all particle in one go:
            dV_50 = np.mean(dV_50)  # mean log diameter
            dN_50 = np.exp(
                np.log(dV_50) - 3 *
                Sd**2)  # convert number distribution to volume distribution
            logger.debug(
                'Droplet distribution median diameter dV_50: %f, dN_50: %f ' %
                (dV_50, np.mean(dN_50)))
            spectrum = (np.exp(
                -(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 /
                (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd *
                                 np.sqrt(2 * np.pi))
            self.droplet_spectrum_pdf = spectrum / np.sum(spectrum)
        if ~np.isfinite(np.sum(self.droplet_spectrum_pdf)) or \
                np.abs(np.sum(self.droplet_spectrum_pdf) - 1) > 1e-6:
            logger.warning('Could not update droplet diameters.')
            return self.elements.diameter
        else:
            return np.random.choice(self.droplet_spectrum_diameter,
                                    size=self.num_elements_active(),
                                    p=self.droplet_spectrum_pdf)

    def get_wave_breaking_droplet_diameter_johansen2015(self):
        # Johansen O, Reed M, Bodsberg NR, Natural dispersion revisited
        # DOI: 10.1016/j.marpolbul.2015.02.026
        # requires oil film thickness
        if not hasattr(self, 'droplet_spectrum_pdf') or self.get_config(
                'processes:update_oilfilm_thickness') is True:
            # Generate droplet spectrum as in Johansen et al. (2015)
            # Bounds are hardcoded to 1micron and 3mm
            logger.debug('Generating wave breaking droplet size spectrum')
            self.droplet_spectrum_diameter = np.linspace(1e-6, 3e-3, 1000000)
            g = 9.81
            if self.oil_weathering_model == 'cicese':
                interfacial_tension = self.elements.IFT
            else:
                interfacial_tension = self.oil_water_interfacial_tension
            H = self.significant_wave_height(
            )  # fall height = 2 * wave amplitude
            # Reyolds number (Eq. 7a from Johansen et al. 2015)
            re = (self.elements.density * self.elements.oil_film_thickness *
                  (g * H)**0.5) / (self.elements.viscosity *
                                   self.elements.density)
            # Weber number (Eq. 7b from Johansen et al.2015)
            we = (self.elements.density * self.elements.oil_film_thickness *
                  g * H) / interfacial_tension  # Weber number
            A = 2.251  # parameters from Johansen et al. 2015
            Bp = 0.027
            B = A * Bp
            dN_50 = (A * self.elements.oil_film_thickness * we**-0.6) + (
                B * self.elements.oil_film_thickness * re**-0.6)
            # median droplet diameter in number distribution
            sd = 0.4  # log standard deviation in log10 units
            Sd = np.log(10) * sd  # log standard deviation in natural log units
            # Convert number distribution to volume distribution
            dV_50 = np.exp(np.log(dN_50) + 3 * Sd**2)
            # TODO: calculation below with scalars, but we have
            # arrays, with varying oil properties
            # treat all particle in one go:
            dV_50 = np.mean(dV_50)  # mean log diameter
            logger.debug(
                'Droplet distribution median diameter dV_50: %f, dN_50: %f ' %
                (dV_50, np.mean(dN_50)))
            spectrum = (np.exp(
                -(np.log(self.droplet_spectrum_diameter) - np.log(dV_50))**2 /
                (2 * Sd**2))) / (self.droplet_spectrum_diameter * Sd *
                                 np.sqrt(2 * np.pi))
            self.droplet_spectrum_pdf = spectrum / np.sum(spectrum)
        if ~np.isfinite(np.sum(self.droplet_spectrum_pdf)) or \
                np.abs(np.sum(self.droplet_spectrum_pdf) - 1) > 1e-6:
            logger.warning('Could not update droplet diameters.')
            return self.elements.diameter
        else:
            return np.random.choice(self.droplet_spectrum_diameter,
                                    size=self.num_elements_active(),
                                    p=self.droplet_spectrum_pdf)

    #: Viscosity-stability scaling for emulsification, ported from
    #: abkatun_2026/oil_weathering.py emulsification_rate() (Fingas &
    #: Fieldhouse 2004): emulsification is favoured for moderately viscous
    #: oil (stability rises through 1x-3x as nu_dry goes from
    #: EMUL_STABILITY_NU_REF to 3*EMUL_STABILITY_NU_REF) and suppressed for
    #: very fluid, fresh oil (stability floors at EMUL_STABILITY_MIN).
    #:   k_emul_eff = k_emul * clip(nu_dry/EMUL_STABILITY_NU_REF,
    #:                               EMUL_STABILITY_MIN, EMUL_STABILITY_MAX)
    #:                       * EMUL_RATE_SCALE
    #: nu_dry = dry-oil (no-emulsion), evaporation-corrected kinematic
    #: viscosity [cSt] (self._oil_viscosity_dry, set in weather_cicese()).
    #: Set EMUL_STABILITY_MIN = EMUL_STABILITY_MAX = EMUL_RATE_SCALE = 1.0
    #: to recover the original (viscosity-independent) NOAA Bullwinkle
    #: behaviour.
    EMUL_STABILITY_NU_REF = 1000.0  # cSt
    EMUL_STABILITY_MIN = 0.3
    EMUL_STABILITY_MAX = 3.0

    #: Overall multiplicative correction to the stock ADIOS2/NOAA Bullwinkle
    #: water-uptake rate constant k_emul = water_uptake_coefficient(...)
    #: (noaa_oil_weathering.py: k_emul = 6*K0Y*U^2/drop_max, K0Y=2.024e-6,
    #: drop_max=1e-5).
    #:
    #: Why this is needed (Abkatun 2026 calibration, June 2026):
    #: The interfacial-area update in emulsification_cicese() drives
    #:   water_fraction = IA*drop_max / (6 + IA*drop_max)
    #: which HALF-SATURATES already at IA = 6/drop_max = 6e5. The
    #: exp(-k_emul/S_max * age) term that is meant to slow IA's growth as
    #: it approaches S_max = (6/drop_min)*(Y_max/(1-Y_max)) ~ 5.4e7 (at
    #: Y_max=0.9, drop_min=1e-6) therefore never engages -- water_fraction
    #: is already saturated at ~1% of S_max. The result is a SINGLE
    #: characteristic timescale
    #:   t_char = (6/drop_max) / k_emul_eff
    #: which, for unscaled k_emul (~121 at 10 m/s wind) with
    #: EMUL_STABILITY in [0.3, 1.0], gives t_char ~ 1.4-4.6 HOURS --
    #: i.e. mousse forms within hours regardless of the oil's weathering
    #: state.
    #:
    #: For MAYA crude, Fingas & Fieldhouse (2004)-consistent emulsification
    #: is days-to-weeks: abkatun_2026/oil_weathering.py's independent ODE
    #: (K_EM = 2.0e-9 s^-1 (m/s)^-2, stability capped at 3.0) gives
    #: t_char ~ 17 days at 10 m/s -- matching the day-17.2 tar-ball
    #: saturation target from TarballFormation_MPB.tex. A first analytical
    #: estimate (t_char(10 m/s) ~ 17 days => k_emul_needed/k_emul_base(10 m/s)
    #: ~ 0.41/121.4 ~ 3e-3) was in the right ballpark and consistent with the
    #: K_SCALE = 2e-3 empirically calibrated for CICOILv1.0 against the same
    #: day-17.2 target, but underestimated the required rate by ~7x once
    #: tested in the full ciceseoil/PyGNOME run (the analytical estimate
    #: ignores the time evolution of wind, evaporation and EMUL_STABILITY
    #: over the run, and the non-linear approach of water_fraction to its
    #: cap).
    #:
    #: EMPIRICAL CALIBRATION (Abkatun 2026, TAMOC+GLORYS realistic
    #: SST/salinity/droplet-size run, calibrate_emulsification_viscosity_tamoc.py,
    #: June 2026), with default EMUL_STABILITY_NU_REF/MIN/MAX above:
    #:   EMUL_RATE_SCALE   tar-ball crossing day (NU_TARBALL=10,000 cSt)
    #:   1.0 (no scaling)  0.76
    #:   0.024             10.21
    #:   0.0142            16.92   <- matches target (17.2) to within 1.6%
    #:   0.0072            none within 25 days (7,388 cSt @ day 25)
    #:   0.002             none within 25 days (2,419 cSt @ day 25)
    #: These four points are consistent with crossing_day ~ const/EMUL_RATE_SCALE
    #: (log-log slope ~ -1.0 between the 0.0072 and 0.024 points).
    #:
    #: EMUL_STABILITY_MIN/MAX alone (a <=10x range, 0.3-3.0) cannot supply
    #: this ~2-orders-of-magnitude rate correction -- hence this separate,
    #: unbounded multiplicative factor. Set EMUL_RATE_SCALE = 1.0 (with
    #: EMUL_STABILITY_MIN = EMUL_STABILITY_MAX = 1.0) to recover the stock
    #: NOAA Bullwinkle rate.
    EMUL_RATE_SCALE = 0.0142

    def emulsification_cicese(self):
        #############################################
        # Emulsification (surface only?)
        #############################################
        logging.debug('    Calculating emulsification - NOAA')
        emul_time = 7200   # self.oiltype.bulltime
        emul_constant = self.oiltype.bullwinkle
        # max water content fraction per spillet. Depending on dispersant treatment
        Y_max = self.elements.max_water
        # emulsion
        emulsifiable = np.where(self.elements.max_water > 0.)
        if len(emulsifiable) == 0:
            logging.debug('Oil does not emulsify, returning.')
            return
        # Constants for droplets
        drop_min = 1.0e-6
        drop_max = 1.0e-5
        S_max = (6. / drop_min) * (Y_max / (1.0 - Y_max))
        S_min = (6. / drop_max) * (Y_max / (1.0 - Y_max))

        # Emulsify...
        # Modified to take cicese fraction calculation into account
        fraction_evaporated = self.elements.fraction_evaporated

        start_emulsion = np.where(
            ((self.elements.age_seconds >= emul_time) & (emul_time >= 0)) |
            ((fraction_evaporated >= emul_constant) & (emul_constant > 0))
            )[0]
        if len(start_emulsion) == 0:
            logging.debug('        Emulsification not yet started')
            return

        if self.oiltype.bulltime > 0:  # User has set value
            start_time = self.oiltype.bulltime*np.ones(len(start_emulsion))
        else:
            start_time = self.elements.age_seconds[start_emulsion]
            start_time[self.elements.age_seconds[start_emulsion]
                       >= 0] = self.elements.bulltime[start_emulsion]
        # Update droplet interfacial area
        k_emul = noaa.water_uptake_coefficient(
                    self.oiltype, self.wind_speed()[start_emulsion])

        # Viscosity-stability scaling + overall rate correction (see
        # EMUL_STABILITY_* / EMUL_RATE_SCALE docstrings above)
        nu_dry_cst = self._oil_viscosity_dry[start_emulsion] * 1.0e6  # m^2/s -> cSt
        stability = np.clip(nu_dry_cst / self.EMUL_STABILITY_NU_REF,
                             self.EMUL_STABILITY_MIN, self.EMUL_STABILITY_MAX)
        k_emul = k_emul * stability * self.EMUL_RATE_SCALE

        self.elements.interfacial_area[start_emulsion] = \
            self.elements.interfacial_area[start_emulsion] + \
            (k_emul*self.time_step.total_seconds()*
             np.exp((-k_emul/S_max[start_emulsion])*(
                self.elements.age_seconds[start_emulsion] - start_time)))
        self.elements.interfacial_area[start_emulsion] = np.minimum(self.elements.interfacial_area[start_emulsion],
                                                                    S_max[start_emulsion])

        # Update water fraction
        self.elements.water_fraction[start_emulsion] = (
            self.elements.interfacial_area[start_emulsion]*drop_max/
            (6.0 + (self.elements.interfacial_area[start_emulsion]
             *drop_max)))
        over_max = np.where(self.elements.interfacial_area >= ((6.0 / drop_max)*(Y_max/(1.0 - Y_max))))
        self.elements.water_fraction[over_max] = Y_max[over_max]
        # self.elements.water_fraction[self.elements.interfacial_area >=
        #     ((6.0 / drop_max)*(Y_max/(1.0 - Y_max)))] = Y_max


    # ── Biodegradation ─────────────────────────────────────────────────────────

    #: Reference temperature (°C) and Q10 factor for the temperature-scaling
    #: law applied to all per-component biodegradation rate constants below:
    #:   rate(T) = rate(T_REF) * BIODEG_Q10 ** ((T - T_REF) / 10)
    #: Q10=3, T_REF=20°C reproduces the Adcroft et al. (2010) GoM
    #: dissolved-oil decay timescale tau(T) = 12 * 3**((20-T)/10) days
    #: (62 d @ 5°C, 12 d @ 20°C, 7 d @ 25°C).
    BIODEG_T_REF = 20.0
    BIODEG_Q10 = 3.0

    #: Per-group base biodegradation rate constants (d⁻¹) at T_REF=20°C for
    #: the 8 aliphatic pseudo-groups (G1=light/<C10 .. G8=heavy waxes/C36+),
    #: for DISPERSED/ENTRAINED oil (droplets in the water column, z<0; high
    #: surface-area-to-volume for microbial colonization). Half-lives 5-40 d,
    #: consistent with Prince et al. (2003, 2017 review "The Rate of Crude
    #: Oil Biodegradation in the Sea"): dispersed oil at environmentally
    #: relevant (sub-ppm) concentrations has an apparent half-life of
    #: 1-3 weeks, with mesocosm aliphatic (C8-C40) bulk half-lives ~8-19 d;
    #: C10-C22 most bioavailable/fastest, <C14 fastest, heavier waxy
    #: components slower. Applied to surface-slick particles (z=0) after
    #: BIODEG_SURFACE_FACTOR suppression -- see below.
    K_BIODEG_ALI = np.array([
        0.1386,  # G1  t½ ~5 d
        0.1155,  # G2  t½ ~6 d
        0.0990,  # G3  t½ ~7 d   (C10-C22 range, most bioavailable)
        0.0866,  # G4  t½ ~8 d
        0.0693,  # G5  t½ ~10 d
        0.0462,  # G6  t½ ~15 d
        0.0277,  # G7  t½ ~25 d  (waxes, slower)
        0.0173,  # G8  t½ ~40 d  (heavy waxes)
    ], dtype=np.float32)

    #: Per-group base biodegradation rate constants (d⁻¹) at T_REF=20°C for
    #: the 8 aromatic pseudo-groups (G1=benzene-like .. G8=pyrene+), for
    #: DISPERSED/ENTRAINED oil (z<0) -- same basis as K_BIODEG_ALI
    #: (Prince et al. dispersed-oil 1-3 week half-life range). Half-lives
    #: 8-60 d. Ordering follows Bacosa et al. (2015): biodegradation is
    #: dominant/comparable to photooxidation for heavy (4-5 ring) PAHs
    #: (G7/G8 half-lives here are shorter than the corresponding
    #: K_PHOTO_AROM half-lives), with light aromatics (G1-G3, BTEX-like)
    #: intermediate.
    K_BIODEG_AROM = np.array([
        0.0866,  # G1  t½ ~8 d   benzene-like
        0.0693,  # G2  t½ ~10 d  toluene-like
        0.0578,  # G3  t½ ~12 d  xylenes
        0.0462,  # G4  t½ ~15 d  naphthalene (2-ring)
        0.0347,  # G5  t½ ~20 d  acenaphthene (2-3 ring)
        0.0277,  # G6  t½ ~25 d  phenanthrene (3 ring)
        0.0173,  # G7  t½ ~40 d  chrysene (4 ring)
        0.0116,  # G8  t½ ~60 d  pyrene+ (4-5 ring)
    ], dtype=np.float32)

    #: Residue (asphaltenes+resins, comp_residue) -- essentially
    #: non-biodegradable on simulation timescales (t½ ~500 d, dispersed-oil
    #: basis; even longer at the surface after BIODEG_SURFACE_FACTOR).
    K_BIODEG_RESIDUE = 0.00139

    #: Surface-slick biodegradation suppression factor (dimensionless,
    #: applied multiplicatively to K_BIODEG_* for particles with z==0).
    #:
    #: Prince et al. (2003, 2017 review "The Rate of Crude Oil
    #: Biodegradation in the Sea") report that floating surface slicks are
    #: "almost immune to detectable biodegradation" because of their low
    #: surface-area-to-volume ratio for microbial colonization: apparent
    #: half-lives are "many months to many years", vs. 1-3 weeks for
    #: dispersed oil at the same bulk concentration -- roughly a 1-2 order
    #: of magnitude difference. BIODEG_SURFACE_FACTOR=1/30 maps the
    #: dispersed-oil K_BIODEG_* half-lives (5-60 d) onto ~5 months -- ~5
    #: years at the surface, spanning the reported "months to years" range.
    BIODEG_SURFACE_FACTOR = 1.0 / 30.0

    #: Reference droplet diameter (m) at which K_BIODEG_ALI/AROM/RESIDUE are
    #: calibrated -- matches `seed:droplet_diameter_mu` default (1 mm),
    #: representative of the "dispersed/entrained oil" droplets referenced
    #: by Adcroft et al. (2010) / Prince et al. (2003, 2017).
    BIODEG_DROPLET_D_REF = 0.001  # m

    #: Exponent for the droplet-size surface-area-to-volume scaling factor
    #: K_i_eff = K_i * (BIODEG_DROPLET_D_REF / d_i) ** BIODEG_DROPLET_EXPONENT.
    #: Exponent=1 follows from SA/V ~ 1/d for a sphere: smaller droplets
    #: present proportionally more oil-water interface per unit oil volume
    #: to colonizing bacteria (e.g. King 1992 microbial degradation
    #: kinetics; Brakstad et al. 2015 droplet-size biodegradation
    #: experiments).
    BIODEG_DROPLET_EXPONENT = 1.0

    #: `elements.diameter` is clamped to this range (m) before computing the
    #: droplet-size scaling factor, bounding the factor to [0.1, 100] for
    #: BIODEG_DROPLET_D_REF=1 mm. Avoids unphysically extreme rates for
    #: sub-resolution droplets from the wave-entrainment spectrum (down to
    #: ~1 um) or unusually large subsea droplets.
    BIODEG_DROPLET_DIAMETER_MIN = 1e-5  # m, 10 um
    BIODEG_DROPLET_DIAMETER_MAX = 1e-2  # m, 1 cm

    #: Dissolved-phase biodegradation rate constants (d⁻¹) at T_REF=20°C,
    #: applied per pseudo-component to `cicese_mass_balance['mass_dissolved']`
    #: for is_dissolved==1 particles (e.g. TAMOC far-field dissolved plume
    #: mass, seeded via seed_plume_elements). Dissolved hydrocarbons are the
    #: most bioavailable fraction (fully solubilized, maximal
    #: surface-area-to-volume for microbial uptake), so the relative
    #: ordering across the 17 pseudo-components follows the same
    #: lighter/more-soluble-degrades-faster pattern as K_BIODEG_ALI/AROM/
    #: RESIDUE for dispersed oil, rescaled by K_BIODEG_WATER_ANCHOR so that
    #: the mean rate across the 16 hydrocarbon groups matches Adcroft et al.
    #: (2010)'s 12-day apparent half-life (k=0.0578 d⁻¹) for dissolved GoM
    #: oil-plume hydrocarbons at T_REF=20°C. Residue is rescaled by the same
    #: factor and remains effectively non-biodegradable in the dissolved
    #: phase too.
    K_BIODEG_WATER_ANCHOR = 0.0578  # d^-1, t1/2 ~ 12 d @ 20C (Adcroft et al. 2010)
    K_BIODEG_WATER = (
        K_BIODEG_WATER_ANCHOR
        / np.mean(np.concatenate([K_BIODEG_ALI, K_BIODEG_AROM]))
        * np.concatenate([K_BIODEG_ALI, K_BIODEG_AROM, [K_BIODEG_RESIDUE]])
    ).astype(np.float32)  # (17,)

    #: Fraction of the initial dissolved-pool mass (elements.mass_dissolved_subsea)
    #: below which a dissolved (is_dissolved==1) element is deactivated by
    #: handle_subsea_dissolution(), once mass_dissolved has decayed via
    #: K_BIODEG_WATER. Avoids tracking near-zero-mass dissolved particles
    #: indefinitely.
    BIODEG_WATER_DEPLETION_FRAC = 0.01

    def biodegradation_cicese(self):
        """Biodegradation of oil-phase pseudo-components by marine bacteria.

        Process:
            dM_i = M_i * (1 - exp(-k_i(T) * dt))
            k_i(T) = k_i(T_REF) * BIODEG_Q10 ** ((T - T_REF) / 10)

        where k_i (d⁻¹) is the per-group base rate constant (K_BIODEG_ALI /
        K_BIODEG_AROM / K_BIODEG_RESIDUE), T is the local sea water
        temperature (°C), and dt is the timestep in days. The exponential
        form (rather than linear k*dt) matches the upstream Adcroft/
        half_time formulas and remains stable for large dt.

        All 17 pseudo-components (8 aliphatic + 8 aromatic + residue) of
        `cicese_mass_balance['mass_components']` are affected, for all
        active particles (surface slick and droplets) -- unlike
        photooxidation, biodegradation is not light- or depth-limited.

        K_BIODEG_ALI/AROM/RESIDUE are calibrated for DISPERSED/ENTRAINED oil
        (z<0, droplets with high surface-area-to-volume for microbial
        colonization). For surface-slick particles (z==0), k_i is further
        multiplied by BIODEG_SURFACE_FACTOR (1/30), following Prince et al.
        (2003, 2017): floating slicks are "almost immune to detectable
        biodegradation" (apparent half-life "many months to many years"),
        vs. 1-3 weeks for dispersed oil at comparable bulk concentration.

        Droplet-size scaling: k_i is additionally multiplied by
        (BIODEG_DROPLET_D_REF / d) ** BIODEG_DROPLET_EXPONENT, where d is
        `elements.diameter` clamped to
        [BIODEG_DROPLET_DIAMETER_MIN, BIODEG_DROPLET_DIAMETER_MAX].
        K_BIODEG_ALI/AROM/RESIDUE are calibrated at BIODEG_DROPLET_D_REF
        (1 mm); smaller droplets (higher surface-area-to-volume) degrade
        faster, larger droplets slower. Particles with diameter==0 (the
        default for most surface "spillet" particles, for which diameter is
        not a meaningful droplet size) are treated as if at
        BIODEG_DROPLET_D_REF, i.e. factor=1 -- surface suppression is
        already handled separately via BIODEG_SURFACE_FACTOR above.

        Dissolved-phase pathway (`mass_biodegraded_from_water`):
        is_dissolved==1 particles (e.g. TAMOC far-field dissolved plume mass)
        have their initial mass transferred into
        `cicese_mass_balance['mass_dissolved']` by handle_subsea_dissolution()
        (called earlier in oil_weathering_cicese, on the seeding step) instead
        of being deactivated immediately. Here, that dissolved pool decays
        per pseudo-component at rate K_BIODEG_WATER (Q10-scaled, no surface
        suppression -- dissolved mass is never at the surface), moving mass into
        `mass_biodegraded_from_water`. Once the dissolved pool falls below
        BIODEG_WATER_DEPLETION_FRAC of its initial mass, the element is
        deactivated by handle_subsea_dissolution() on a subsequent step.

        References
        ----------
        Adcroft et al. (2010) Geophys. Res. Lett. 37:L18605
            (Q10=3, T_REF=20°C temperature scaling)
        Bacosa et al. (2015) Mar. Pollut. Bull. 95:265
            (biodegradation vs. photooxidation magnitudes by compound class)
        Prince et al. (2003, 2017) Environ. Sci. Technol.
            "The Rate of Crude Oil Biodegradation in the Sea"
            (dispersed-oil vs. surface-slick biodegradation rate contrast)
        """
        logging.debug('    Calculating: biodegradation - CICESE')

        active = np.arange(len(self.elements.lon))
        if len(active) == 0:
            return

        dt_days = self.time_step.total_seconds() / (3600 * 24)

        # Temperature scaling (Adcroft Q10=3, T_REF=20°C); env temp is in K
        # at this point in oil_weathering_cicese (converted C->K at its top).
        swt_C = self.environment.sea_water_temperature[active] - 273.15
        q10_factor = self.BIODEG_Q10 ** ((swt_C - self.BIODEG_T_REF) / 10.0)  # (n_active,)

        # Per-group rate array: [8 aliphatic, 8 aromatic, 1 residue] = 17
        K_base = np.concatenate([self.K_BIODEG_ALI, self.K_BIODEG_AROM,
                                  [self.K_BIODEG_RESIDUE]])  # (17,)
        K = K_base[np.newaxis, :] * q10_factor[:, np.newaxis]  # (n_active, 17)

        # Surface-slick suppression (Prince et al.): K_BIODEG_* are
        # calibrated for dispersed/entrained oil (z<0); particles at the
        # surface (z==0) biodegrade BIODEG_SURFACE_FACTOR times slower.
        is_surface = self.elements.z[active] == 0  # (n_active,)
        surface_mult = np.where(is_surface, self.BIODEG_SURFACE_FACTOR, 1.0)
        K = K * surface_mult[:, np.newaxis]

        # Droplet-size scaling (audit #2): K_i_eff = K_i * (d_ref/d)^n.
        # diameter==0 (unset, e.g. most surface "spillet" particles, for
        # which diameter is not a meaningful droplet size) is treated as
        # d_ref -> factor=1, leaving K unchanged.
        diameter = self.elements.diameter[active].copy()
        n_unset = int((diameter <= 0).sum())
        if n_unset > 0:
            logger.debug('biodegradation_cicese: %d particles with diameter<=0 '
                         'replaced with d_ref=%.4f m', n_unset, self.BIODEG_DROPLET_D_REF)
        diameter[diameter <= 0] = self.BIODEG_DROPLET_D_REF
        diameter = np.clip(diameter, self.BIODEG_DROPLET_DIAMETER_MIN,
                            self.BIODEG_DROPLET_DIAMETER_MAX)
        size_factor = (self.BIODEG_DROPLET_D_REF / diameter) ** self.BIODEG_DROPLET_EXPONENT
        K = K * size_factor[:, np.newaxis]

        activeID = self.elements.ID[active]  # 0-based (BC-2 fix)
        assert activeID.min() >= 0, (
            f'Element IDs must be >= 0 for cicese_mass_balance indexing, '
            f'got min={activeID.min()}')
        mass_comp = self.cicese_mass_balance['mass_components'][activeID, :]  # (n_active, 17)

        dM_biodeg = mass_comp * (1.0 - np.exp(-K * dt_days))
        dM_biodeg = np.clip(dM_biodeg, 0.0, mass_comp)

        # Update component mass balance
        self.cicese_mass_balance['mass_components'][activeID, :] -= dM_biodeg
        self.cicese_mass_balance['mass_biodegraded_from_oil'][activeID, :] += dM_biodeg

        # Update per-particle totals
        total_biodeg = np.sum(dM_biodeg, axis=1)  # (n_active,)
        self.elements.mass_biodegraded_from_oil[active] += total_biodeg
        self.elements.mass_biodegraded[active]          += total_biodeg  # generic var, used by reporting/plots
        self.elements.mass_oil[active]                  -= total_biodeg

        # Update fraction_biodegraded_from_oil
        initial_mass = (self.elements.mass_oil[active] +
                        self.elements.mass_evaporated[active] +
                        self.elements.mass_dispersed[active] +
                        self.elements.mass_photooxidized[active] +
                        self.elements.mass_biodegraded[active] +
                        self.elements.mass_op_dissolved[active] +
                        self.elements.mass_op_degraded[active])
        mask = initial_mass > 0
        self.elements.fraction_biodegraded_from_oil[active[mask]] = (
            self.elements.mass_biodegraded_from_oil[active[mask]] /
            initial_mass[mask])

        logging.debug(
            f'    Biodeg: {total_biodeg.sum():.3f} kg biodegraded this step '
            f'from {len(active)} elements')

        # Dissolved-phase biodegradation (mass_biodegraded_from_water): decays
        # cicese_mass_balance['mass_dissolved'], populated by
        # handle_subsea_dissolution() for is_dissolved==1 elements.
        dissolved = np.where(self.elements.is_dissolved[active] == 1)[0]
        if len(dissolved) > 0:
            dis_idx = active[dissolved]
            disID = self.elements.ID[dis_idx]

            K_water = (self.K_BIODEG_WATER[np.newaxis, :]
                       * q10_factor[dissolved][:, np.newaxis])  # (n_dissolved, 17)
            mass_dis = self.cicese_mass_balance['mass_dissolved'][disID, :]  # (n_dissolved, 17)

            dM_water = mass_dis * (1.0 - np.exp(-K_water * dt_days))
            dM_water = np.clip(dM_water, 0.0, mass_dis)

            self.cicese_mass_balance['mass_dissolved'][disID, :] -= dM_water
            self.cicese_mass_balance['mass_biodegraded_from_water'][disID, :] += dM_water

            total_water = np.sum(dM_water, axis=1)  # (n_dissolved,)
            self.elements.mass_biodegraded_from_water[dis_idx] += total_water
            self.elements.mass_biodegraded[dis_idx]            += total_water

            initial_dissolved = self.elements.mass_dissolved_subsea[dis_idx]
            mask_w = initial_dissolved > 0
            self.elements.fraction_biodegraded_from_water[dis_idx[mask_w]] = (
                self.elements.mass_biodegraded_from_water[dis_idx[mask_w]] /
                initial_dissolved[mask_w])

            logging.debug(
                f'    Biodeg (dissolved): {total_water.sum():.3f} kg biodegraded '
                f'this step from {len(dissolved)} dissolved elements')

    # ── Photooxidation ─────────────────────────────────────────────────────────

    #: Per-group photoox rate constants (m² J⁻¹) for aromatic groups G1–G8.
    #: G1/G2/G3 from surface_fate_coupled.py K_PHOTO dict (unchanged, not
    #: part of the rescaling below; G3=xylenes already calibrated there for
    #: real ERA5 UV — no compound-specific literature alternative was found,
    #: see TODO #3 note below). G4-G8 anchors (naphthalene/phenanthrene/
    #: pyrene+) calibrated against full 30-day surface_fate_coupled.py time
    #: series (data_local/surface_fate/surface_fate_base.csv), 2026-06-13:
    #: least-squares rescale f*=0.3834 of the prior synthetic-UV-calibrated
    #: values (previously f=0.4046 from a single day-5 snapshot match).
    #: Full-series RMSE improves 0.822 -> 0.782 percentage points; a residual
    #: shape mismatch remains (depth-resolved scheme saturates by ~day 15,
    #: target rises gradually to day 30) -- documented as a known limitation.
    #: G5 (acenaphthene) and G7 (chrysene) have no direct surface_fate
    #: equivalent and remain log-linearly interpolated vs. boiling point
    #: between their bracketing recalibrated anchors (G4-G6 and G6-G8
    #: respectively). A WebSearch (2026-06-13) for compound-specific
    #: oil-slick photooxidation rate constants for xylenes/acenaphthene/
    #: chrysene found none; the resulting half-lives (G5 ~7.8 d, G7 ~25.3 d)
    #: fall within the Bacosa et al. (2015) 2-3-ring (0.7-8.7 d) and 4-5-ring
    #: (10-69 d) half-life ranges respectively, so interpolation is retained.
    #: Groups are ordered from lightest (G1, benzene-like) to heaviest
    #: (G8, pyrene-like) boiling point.
    K_PHOTO_AROM = np.array([
        2.9e-8,    # G1  benzene       bp ~80°C   t½ ~14 d
        2.0e-8,    # G2  toluene       bp ~111°C  t½ ~20 d
        5.0e-8,    # G3  xylenes       bp ~138°C  t½ ~8 d    (surface_fate K_PHOTO, unchanged)
        4.260e-8,  # G4  naphthalene   bp ~218°C  t½ ~9.4 d  (full-series recal 2026-06-13, f*=0.3834)
        5.145e-8,  # G5  acenaphth.    bp ~280°C  t½ ~7.8 d  (interpolated G4-G6)
        6.176e-8,  # G6  phenanthrene  bp ~340°C  t½ ~6.5 d  (full-series recal 2026-06-13, f*=0.3834)
        1.581e-8,  # G7  chrysene      bp ~375°C  t½ ~25.3 d (interpolated G6-G8)
        5.110e-9,  # G8  pyrene+       bp ~404°C  t½ ~79.8 d (full-series recal 2026-06-13, f*=0.3834)
    ], dtype=np.float64)

    #: UV diffuse attenuation coefficient (m⁻¹) for the water column, used to
    #: scale the surface UV irradiance to each particle's depth via
    #: Beer-Lambert: I_UV(z) = I_UV(0) * exp(-K_D_UV * |z|).
    #: Value from surface_fate_coupled.py ("clear GoM" condition).
    K_D_UV = 0.10  # m⁻¹

    #: Depth (m) below which attenuated UV is considered negligible
    #: (exp(-K_D_UV * UV_DEPTH_CUTOFF) ≈ 0.007, i.e. <1% of surface value)
    #: -- elements deeper than this are skipped for photooxidation.
    UV_DEPTH_CUTOFF = 50.0  # m

    #: Fraction of newly-formed photooxidized aromatic mass that partitions
    #: into the water column as dissolved oxygenated photoproducts (OPs:
    #: oxy-PAHs, quinones, carboxylic acids, etc.) rather than remaining as
    #: an oil-phase residue (`mass_photooxidized`, terminal sink). OPs are
    #: typically markedly more water-soluble than their parent PAHs (Aeppli
    #: et al. 2012), so an even 50/50 split between "stays in/on oil" and
    #: "partitions to water column" is taken as a first-pass estimate
    #: pending compound-class-resolved partitioning data.
    OP_DISSOLVED_FRACTION = 0.5

    #: First-order removal rate (d⁻¹) at T_REF=20°C (Q10-scaled like
    #: K_BIODEG_*) for the dissolved oxygenated photoproduct pool
    #: (`mass_op_dissolved`), representing further
    #: biodegradation/photodegradation/dilution of these compounds in the
    #: water column. OPs are smaller and more polar than their parent PAHs
    #: (Aeppli et al. 2012) and are generally more bioavailable/labile, so
    #: K_OP_REMOVAL is set to 2x K_BIODEG_WATER_ANCHOR (t½~12d@20°C) ->
    #: t½~6d@20°C, a first-pass estimate pending OP-specific
    #: biodegradation-rate literature.
    K_OP_REMOVAL = 0.116  # d^-1, t1/2 ~ 6 d @ 20C

    def _compute_uv_irradiance(self):
        """Solar UV irradiance just below the sea surface (W m⁻²).

        Prefers real UV forcing from `surface_downward_uv_radiation`
        (e.g. ERA5 'uvb', converted from J m⁻² hourly-accumulated to W m⁻²
        and CF-tagged so `reader_netCDF_CF_generic` can supply it), averaged
        over all active particles (this is a 2D surface field and does not
        depend on particle depth). Negative/missing values (the -1.0
        fallback, meaning no reader provided this variable) fall back to a
        synthetic diurnal cosine solar-zenith model calibrated to Bay of
        Campeche July conditions (peak ~50 W m⁻², UV-A/B effective for
        photooxidation; calibrated against naphthalene t½ ~2 d).

        This is the boundary value at z=0; `photooxidation_cicese()` applies
        Beer-Lambert depth attenuation (`K_D_UV`) to obtain the local
        irradiance at each particle's depth.

        Returns
        -------
        float
            UV irradiance in W m⁻² at the sea surface. Zero at nighttime.
        """
        if len(self.elements.lon) == 0:
            return 0.0

        uv_env = getattr(self.environment, 'surface_downward_uv_radiation', None)
        if uv_env is not None:
            uv_all = np.asarray(uv_env)
            valid = uv_all[uv_all >= 0.0]
            if len(valid) > 0:
                return max(0.0, float(np.mean(valid)))

        t = self.time
        lon_mean = float(np.mean(self.elements.lon))
        # Solar time offset from UTC (hours) — approx from longitude
        utc_offset = lon_mean / 15.0
        local_hour = (t.hour + t.minute / 60.0 + utc_offset) % 24.0
        # Cosine solar-zenith model; peak at local noon
        cos_z = np.cos(np.radians((local_hour - 12.0) * 15.0))
        I_UV = max(0.0, float(cos_z)) * 50.0  # W m⁻², peak ~50 at noon
        return I_UV


    def photooxidation_cicese(self):
        """Photooxidation of aromatic pseudo-components by solar UV.

        Process:
            dM_i = k_i · E_UV_eff(z) · M_i
            E_UV_eff(z) = I_UV(0) · exp(-K_D_UV · |z|) · Δt

        where k_i (m² J⁻¹) is the per-group rate constant, I_UV(0) is the
        solar UV irradiance just below the sea surface (real ERA5 or
        synthetic diurnal), K_D_UV (m⁻¹) is the water-column UV diffuse
        attenuation coefficient, z is the particle depth (≤0, m), and Δt is
        the timestep. Elements deeper than `UV_DEPTH_CUTOFF` are skipped
        (negligible attenuated dose).

        Only the 8 aromatic pseudo-groups (indices 8–15 in mass_components)
        are affected. Aliphatic groups and the residue are not photooxidized.

        Oxygenated photoproduct (OP) tracking
        --------------------------------------
        Of the photooxidized mass dM_phox removed from the oil-phase
        aromatics each step, a fraction OP_DISSOLVED_FRACTION partitions
        into the water column as dissolved oxygenated photoproducts
        (mass_op_dissolved); the remainder (1 - OP_DISSOLVED_FRACTION)
        stays in mass_photooxidized as before (an oil-phase terminal sink).
        This reflects the markedly higher water solubility of OPs (oxy-PAHs,
        quinones, carboxylic acids) relative to their parent PAHs (Aeppli
        et al. 2012).

        The dissolved OP pool decays at its own first-order rate
        K_OP_REMOVAL (Q10-scaled like K_BIODEG_*, t1/2~6 d @ 20C),
        representing further biodegradation/photodegradation/dilution in the
        water column; removed mass accumulates in mass_op_degraded. This
        removal step is unconditional (computed before the nighttime
        I_UV0<=0 check), since it is a water-column process independent of
        solar UV.

        fraction_op_dissolved and fraction_op_degraded are reported relative
        to the same (conserved) initial_mass denominator used for
        fraction_photooxidized, which now also includes mass_op_dissolved
        and mass_op_degraded.

        Note on f_surface
        ------------------
        No separate `_compute_f_surface_cicese()` discount is applied (cf.
        surface_fate_coupled.py's fixed-w_rise f_surface). In CICOIL,
        OpenDrift's vertical_mixing + per-particle terminal_velocity (from
        the wave-entrainment droplet spectrum) already places particles at
        z<0 when submerged, and the Beer-Lambert attenuation above is the
        continuous, depth-resolved realization of that same wave-driven
        submergence effect (replacing the previous binary z==0 cutoff).
        Re-applying a 0-D f_surface on top would double-count submergence
        with an inconsistent (fixed w_rise=0.27 m/s) rise velocity.

        References
        ----------
        Aeppli et al. (2012) Env. Sci. Technol. 46:13093
        Delvigne & Sweeney (1988) Oil & Chem. Pollut. 4:281
        surface_fate_coupled.py rate-constant calibration (CICESE 2026)
        """
        logging.debug('    Calculating: photooxidation - CICESE')

        active = np.where(self.elements.z >= -self.UV_DEPTH_CUTOFF)[0]
        if len(active) == 0:
            logging.debug('All elements below UV depth cutoff, no photooxidation')
            return

        activeID = self.elements.ID[active]  # 0-based (BC-2 fix)
        assert activeID.min() >= 0, (
            f'Element IDs must be >= 0 for cicese_mass_balance indexing, '
            f'got min={activeID.min()}')
        dt = self.time_step.total_seconds()
        dt_days = dt / (3600.0 * 24.0)

        # --- OP removal (K_OP_REMOVAL): decays the dissolved oxygenated
        # photoproduct pool, Q10-scaled, independent of solar UV/lighting so
        # it also proceeds at night. ---
        swt_C = self.environment.sea_water_temperature[active] - 273.15
        q10_factor = self.BIODEG_Q10 ** ((swt_C - self.BIODEG_T_REF) / 10.0)  # (n_active,)
        K_op = self.K_OP_REMOVAL * q10_factor  # (n_active,)

        mass_op = self.cicese_mass_balance['mass_op_dissolved'][activeID, :]  # (n_active, 17)
        dM_op_removed = mass_op * (1.0 - np.exp(-K_op[:, np.newaxis] * dt_days))
        dM_op_removed = np.clip(dM_op_removed, 0.0, mass_op)

        self.cicese_mass_balance['mass_op_dissolved'][activeID, :] -= dM_op_removed
        self.cicese_mass_balance['mass_op_degraded'][activeID, :] += dM_op_removed

        total_op_removed = np.sum(dM_op_removed, axis=1)  # (n_active,)
        self.elements.mass_op_dissolved[active] -= total_op_removed
        self.elements.mass_op_degraded[active]  += total_op_removed

        # --- Solar UV irradiance at the sea surface ---
        I_UV0 = self._compute_uv_irradiance()  # W m⁻², at z=0
        if I_UV0 <= 0.0:
            logging.debug('Nighttime — no photooxidation (OP removal still applied)')
        else:
            # --- Depth-resolved UV dose (Beer-Lambert attenuation) ---
            depth = np.abs(self.elements.z[active])             # (n_active,) m
            atten = np.exp(-self.K_D_UV * depth)                # (n_active,)
            E_UV_eff = I_UV0 * atten * dt                       # (n_active,) J m⁻²

            logging.debug(
                f'    Photoox: I_UV(0)={I_UV0:.1f} W/m², '
                f'atten range=[{atten.min():.3f}, {atten.max():.3f}], '
                f'E_UV_eff range=[{E_UV_eff.min():.1f}, {E_UV_eff.max():.1f}] J/m²')

            # --- Apply to aromatic groups (indices 8–15 in the 17-component array) ---
            mass_arom = self.cicese_mass_balance['mass_components'][np.ix_(activeID, range(8, 16))]

            # Rate: dM_i = k_i * E_UV_eff(z) * M_i  (broadcast over particles and groups)
            dM_phox = self.K_PHOTO_AROM[np.newaxis, :] * E_UV_eff[:, np.newaxis] * mass_arom  # (n_active, 8)
            dM_phox = np.minimum(dM_phox, mass_arom)            # cap at available mass
            dM_phox = np.maximum(dM_phox, 0.0)                  # non-negative

            # Split into a dissolved-OP fraction and an oil-phase-residue
            # fraction (audit #4): OP_DISSOLVED_FRACTION partitions to the
            # water column (mass_op_dissolved), the remainder stays as the
            # terminal-sink mass_photooxidized as before.
            dM_op = dM_phox * self.OP_DISSOLVED_FRACTION
            dM_residue = dM_phox - dM_op

            # Update component mass balance
            self.cicese_mass_balance['mass_components'][
                np.ix_(activeID, range(8, 16))] -= dM_phox
            self.cicese_mass_balance['mass_photooxidized'][
                np.ix_(activeID, range(8, 16))] += dM_residue
            self.cicese_mass_balance['mass_op_dissolved'][
                np.ix_(activeID, range(8, 16))] += dM_op

            # Update per-particle totals
            total_phox = np.sum(dM_phox, axis=1)      # (n_active,)
            total_residue = np.sum(dM_residue, axis=1)
            total_op = np.sum(dM_op, axis=1)
            self.elements.mass_photooxidized[active] += total_residue
            self.elements.mass_op_dissolved[active]  += total_op
            self.elements.mass_oil[active]           -= total_phox

            logging.debug(
                f'    Photoox: {total_phox.sum():.3f} kg photooxidized this step '
                f'from {len(active)} elements within UV depth cutoff '
                f'({total_op.sum():.3f} kg routed to dissolved OP pool)')

        # --- Update fraction_photooxidized / fraction_op_* ---
        # initial_mass is the (conserved) total mass of each particle, now
        # split across mass_oil, mass_evaporated, mass_dispersed,
        # mass_photooxidized, mass_biodegraded, mass_op_dissolved and
        # mass_op_degraded.
        initial_mass = (self.elements.mass_oil[active] +
                        self.elements.mass_evaporated[active] +
                        self.elements.mass_dispersed[active] +
                        self.elements.mass_photooxidized[active] +
                        self.elements.mass_biodegraded[active] +
                        self.elements.mass_op_dissolved[active] +
                        self.elements.mass_op_degraded[active])
        mask = initial_mass > 0
        self.elements.fraction_photooxidized[active[mask]] = (
            self.elements.mass_photooxidized[active[mask]] /
            initial_mass[mask])
        self.elements.fraction_op_dissolved[active[mask]] = (
            self.elements.mass_op_dissolved[active[mask]] /
            initial_mass[mask])
        self.elements.fraction_op_degraded[active[mask]] = (
            self.elements.mass_op_degraded[active[mask]] /
            initial_mass[mask])

    # ── Oil type / dispersant ──────────────────────────────────────────────────

    def set_oiltype(self, oiltype):
        # OD 1.14.9: oil_name_alias removed — use identity if not present
        _alias = getattr(adios, 'oil_name_alias', {})
        oiltype = _alias.get(oiltype, oiltype)
        print('setting oil_type to: ', oiltype)

        self.oil_name = oiltype
        # OD 1.14.9: set_config() requires Mode.Config only; __set_seed_config__
        # accepts [Mode.Config, Mode.Ready] — correct for seed-related keys.
        # This allows set_oiltype() to be called after the first seed_elements()
        # (which advances mode Config→Ready), enabling seed_plume_elements() to
        # call seed_elements() multiple times for different particle classes.
        self.__set_seed_config__('seed:oil_type', oiltype)
        if self.oil_weathering_model == 'cicese':
            oils = adios.oils(query=oiltype)
            max_cuts = 0
            best_oil = None
            best_id = None
            if len(oils) == 0:
                raise ValueError('Oil type "%s" not found in NOAA database' % oiltype)
            elif len(oils) == 1:
                self.oiltype = oils[0]
                best_id = oils[0].id
                cuts = props.get_distillation_cuts(oils[0].oil)
                max_cuts = len(cuts)
                self.fluid_properties = FluidProps(self.oiltype.oil, max_cuts)
            else:
                for oil in oils:
                    cuts = props.get_distillation_cuts(oil.oil)
                    no_cuts = len(cuts)
                    id = oil.id
                    logger.info(f'Checked Oil {oiltype} with Adios ID {id} and found {no_cuts} cuts.')

                    kvis = round(oil.kvis_at_temp(285) * 10 ** 6, 2)
                    print('_____ Name= ' + oil.name + ';   Adios ID= ' + oil.id + ';   API= ' +
                          str(round(oil.oil.metadata.API, 1)) + ';   K.Visc (cSt)= ' + str(
                        kvis) + ';   No of cuts= ' + str(no_cuts))

                    if no_cuts > max_cuts:
                        max_cuts = no_cuts
                        best_oil = oil
                        best_id = oil.id
                self.oiltype = best_oil
                self.fluid_properties = FluidProps(self.oiltype.oil, max_cuts)

            print('Selected Oil ', oiltype, ' with Adios ID ', best_id, ' and ', max_cuts, ' cuts.')
            if max_cuts == 0:
                logger.warning(f'#### Selected Oil {oiltype} with Adios ID {best_id} has no distillation cuts! ###')
        else:
            super(OpenCiceseOil, self).set_oiltype(oiltype)

    def dispersant_efficiency_sampled(self, target_DOR, efficiency, samples):
        # Uses standard normal distribution to calculate actual DOR of spillets so that the ratio of integrated
        # spillets DOR to the target is equal to the application efficiency

        # Approximation of efficiency to number of sigmas
        times_sigma = 29.52748 * efficiency ** 4 - 95.01174 * efficiency ** 3 + 109.75335 * efficiency ** 2 - 58.21448 * \
                      efficiency + 14.12874
        sigmas = np.random.random_sample((samples,)) * times_sigma
        applied_DORs = target_DOR * np.exp(-1/2 * sigmas ** 2)

        return applied_DORs

    def treated_IFT_and_WC(self, dispersant, applied_DORs):
        # Calculation of the oil-water interfacial tension based on the applied DOR

        dispersant = self.dispersants[dispersant]
        IFT_untreated = self.oiltype.oil_water_surface_tension()
        treated_IFTs = IFT_untreated * (dispersant['A'] * 0.5 ** (applied_DORs * 100. / dispersant['t12']) + dispersant['B']) / (
                    dispersant['A'] + dispersant['B'])
        relative_IFTs = treated_IFTs / IFT_untreated

        # max water content calculation after dispersant treatment
        alpha = 0.928
        beta = 0.082
        max_treated_WCs = relative_IFTs * alpha - beta
        if np.isscalar(max_treated_WCs):
            max_treated_WCs = max(max_treated_WCs, beta)
        else:
            max_treated_WCs[max_treated_WCs < beta] = beta

        return treated_IFTs, max_treated_WCs

    def get_max_water_content(self):
        # Calculation of the maximum water content that can be reached based on the type of oil

        max_WC = None
        oil = self.oiltype
        emuls = oil.oil.sub_samples[0].environmental_behavior.emulsions
        if len(emuls) != 0:
            max_WC = emuls[0].water_content.value
        else:
            oiltype = self.oil_name
            oils = adios.oils(query=oiltype)
            for oil in oils:
                emuls = oil.oil.sub_samples[0].environmental_behavior.emulsions
                if len(emuls) != 0:
                    max_WC = emuls[0].water_content.value
            if max_WC == None:
                max_WC = oil.emulsion_water_fraction_max
                if self.oil_name in self._cicoil_max_water_fraction_override:
                    max_water_fraction = self._cicoil_max_water_fraction_override[self.oil_name]
                    logging.debug('Overriding max water fraction with value %f instead of default %f'
                                  % (max_water_fraction, max_WC))
                    max_WC = max_water_fraction
        return max_WC

    def weather_permitting_time(self, lon, lat, treatment_time, wind_limit, time_step):
        # Calculation of the time intervals where wind velocity is below the operational limit

        duration = (treatment_time[1] - treatment_time[0]).total_seconds()
        steps_no = round(duration / time_step)
        start_time = datetime.timestamp(treatment_time[0])
        end_time = datetime.timestamp(treatment_time[1])
        time_stamps = np.linspace(start_time, end_time, num=steps_no + 1)
        dates = np.array([datetime.fromtimestamp(stamp) for stamp in time_stamps])
        wind_u = []
        wind_v = []
        for date in dates:
            # OD 1.14.9: get_environment → self.env.get_environment returning (recarray, profiles, missing)
            _eu, _, _ = self.env.get_environment(['x_wind'], date, np.array([lon]), np.array([lat]), np.array([0.0]))
            _ev, _, _ = self.env.get_environment(['y_wind'], date, np.array([lon]), np.array([lat]), np.array([0.0]))
            wind_u.append(float(_eu['x_wind'][0]))
            wind_v.append(float(_ev['y_wind'][0]))
        wind_u = np.array(wind_u)
        wind_v = np.array(wind_v)
        wind_velocity = np.sqrt(wind_u**2 + wind_v**2)

        allowed = np.where(wind_velocity <= wind_limit)[0]
        forbitten = np.where(wind_velocity > wind_limit)[0]
        allowed_end = np.where(np.diff(allowed) > 1)[0]
        up_times = np.split(allowed, allowed_end + 1)
        forbitten_end = np.where(np.diff(forbitten) > 1)[0]
        down_times = np.split(forbitten, forbitten_end + 1)
        allowed_periods = []
        allowed_winds = []
        forbitten_periods = []
        forbitten_winds = []
        for i in up_times:
            allowed_periods.append([dates[i[0]], dates[i[-1]] + timedelta(seconds=time_step)])
            allowed_winds.append([wind_velocity[i[0]:i[-1]+1]])
            print('Dispersant operations are allowed during: ', allowed_periods[-1])
            print('Allowed winds in the period: ', allowed_winds[-1])
        for i in down_times:
            forbitten_periods.append([dates[i[0]], dates[i[-1]] + timedelta(seconds=time_step)])
            forbitten_winds.append([wind_velocity[i[0]:i[-1]+1]])
            print('Dispersant operations are forbitten during: ', forbitten_periods[-1])
            print('Allowed winds in the period: ', forbitten_winds[-1])
        print(dates[allowed])
        return allowed_periods, forbitten_periods, dates[allowed]

    def seed_with_dispersant(self, lon, lat, seed_time, appl_method, DOR, dispersant='C9500',
                             time_frac_surf=0.5, time_frac_ssdi=1.0, effic_surf=0.6, effic_ssdi=0.95, wind_limit=None,
                             delay=None, stock_kg=None, time_step=None, appl_rate=None, number=5000, *args, **kwargs):
        # Main dispersant treatment controlling method

        spill_rate = kwargs['m3_per_hour']
        if time_step is None:
            time_step = self.get_config('general:time_step_minutes') * 60.
        else:
            time_step = time_step * 3600
        if delay is None:
            delay = 0
        delay = timedelta(hours=delay)
        if appl_rate is None:
            appl_rate = spill_rate * DOR
        if stock_kg is not None:
            stock_steps = round(stock_kg / (appl_rate * time_step))
        if DOR <= 0.0:
            application_rate_fraction = 0.0
        else:
            application_rate_fraction = min(appl_rate / (spill_rate * DOR), 1.0)

        treatment_time = [seed_time[0] + delay, seed_time[1]]
        if treatment_time[0] > treatment_time[1]:
            treatment_time[0] = treatment_time[1]
        treatment_duration = (treatment_time[1] - treatment_time[0]).total_seconds()
        seed_duration = (seed_time[1] - seed_time[0]).total_seconds()
        wind_limit = None
        if wind_limit is None:
            permitted_periods = [treatment_time,]
            application_steps = round(treatment_duration / time_step)
        else:
            permitted_periods, forbitten_periods, permitted_dates = \
                self.weather_permitting_time(lon, lat, treatment_time, wind_limit, time_step)
            application_steps = permitted_dates.size
        if application_steps > stock_steps:
            application_steps = stock_steps
            if wind_limit is None:
                treatment_time[1] = treatment_time[0] + timedelta(seconds=stock_steps * time_step)
                treatment_duration = (treatment_time[1] - treatment_time[0]).total_seconds()
                permitted_periods = [treatment_time,]
            else:
                new_periods = []
                runout_date = permitted_dates[stock_steps - 1] + timedelta(seconds=time_step)
                for period in permitted_periods:
                    if period[1] > runout_date:
                        new_periods.append(period)
                    else:
                        period[1] = runout_date
                        new_periods.append(period)
                        permitted_periods = new_periods
                        break
        if appl_method == 'Surface':
            time_used_fraction = min(time_frac_surf * application_steps * time_step / seed_duration, 1.0)
        else:
            time_used_fraction = min(time_frac_ssdi * application_steps * time_step / seed_duration, 1.0)

        IFT_untreated = self.oiltype.oil_water_surface_tension() * 1000.
        max_WC_untreated = self.get_max_water_content()
        number_treated = round(number * time_used_fraction * application_rate_fraction)
        number_untreated = number - number_treated
        if appl_method == 'SSDI':
            plume = kwargs['plume']
            del kwargs['plume']
            if number_treated > 0:
                IFT_treated, max_WC_treated = self.treated_IFT_and_WC(dispersant, DOR*effic_ssdi)
                kwargs['IFT'] = IFT_treated
                kwargs['max_water'] = max_WC_treated
                kwargs['DOR'] = DOR*effic_ssdi
                self.seed_plume_elements(lon, lat, plume, seed_time, number=number_treated, *args, **kwargs)
            if number_untreated > 0:
                kwargs['IFT'] = IFT_untreated
                kwargs['max_water'] = max_WC_untreated
                kwargs['DOR'] = 0.
                self.seed_plume_elements(lon, lat, plume, seed_time, number=number_untreated, *args, **kwargs)
        elif appl_method == 'Surface':
            for period in permitted_periods:
                period_steps = (period[1] - period[0]).total_seconds() / time_step
                time_fraction = period_steps / application_steps
                treated_per_period = round(number_treated * time_fraction)
                untreated_per_period = round(number_untreated * time_fraction)
                number_per_period = treated_per_period + untreated_per_period
                applied_DORs = self.dispersant_efficiency_sampled(target_DOR=DOR, efficiency=effic_surf,
                                                              samples=treated_per_period)
                treated_IFTs, max_treated_WCs = self.treated_IFT_and_WC(dispersant, applied_DORs)
                if untreated_per_period > 0:
                    max_untreated_WCs = np.ones(untreated_per_period) * max_WC_untreated
                    period_max_WCs = np.concatenate((max_treated_WCs, max_untreated_WCs))
                    untreated_IFTs = np.ones(untreated_per_period) * IFT_untreated
                    period_IFTs = np.concatenate((treated_IFTs, untreated_IFTs))
                    untreated_DORs = np.zeros(untreated_per_period)
                    period_DORs = np.concatenate((applied_DORs, untreated_DORs))
                else:
                    period_max_WCs = max_treated_WCs
                    period_IFTs = treated_IFTs
                    period_DORs = applied_DORs
                    number_per_period = treated_per_period
                self.seed_elements(lon, lat, time=period, IFT=period_IFTs, DOR=period_DORs, max_water=period_max_WCs,
                                     number=number_per_period, *args, **kwargs)


    def seed_elements(self, *args, **kwargs):

        if len(args) == 2:
            kwargs['lon'] = args[0]
            kwargs['lat'] = args[1]
            args = {}

        self.store_oil_seed_metadata(**kwargs)

        if 'number' not in kwargs:
            number = self.get_config('seed:number')
        else:
            number = kwargs['number']
        if 'diameter' in kwargs:
            logger.info('Droplet diameter is provided, and will '
                        'be kept constant during simulation')
            self.keep_droplet_diameter = True
        else:
            self.keep_droplet_diameter = False
        if 'z' not in kwargs or kwargs['z'] is None:
            if self.get_config('seed:seafloor') is True:
                kwargs['z'] = 'seafloor'
            else:
                kwargs['z'] = self.get_config('seed:z')
        if isinstance(kwargs['z'], str) and \
                kwargs['z'][0:8] == 'seafloor':
            z = -np.ones(number)
        else:
            z = np.atleast_1d(kwargs['z'])
        if len(z) == 1:
            z = z * np.ones(number)  # Convert scalar z to array
        subsea = z < 0
        if np.sum(subsea) > 0 and 'diameter' not in kwargs:
            # Droplet min and max for particles seeded below sea surface
            sub_dmin = self.get_config('seed:droplet_diameter_min_subsea')
            sub_dmax = self.get_config('seed:droplet_diameter_max_subsea')
            logger.info('Using particle diameters between %s and %s m for '
                        'elements seeded below sea surface.' %
                        (sub_dmin, sub_dmax))
            kwargs['diameter'] = np.random.uniform(sub_dmin, sub_dmax, number)

        if 'oiltype' in kwargs:
            logger.warning(
                'Seed argument *oiltype* is deprecated, use *oil_type* instead'
            )
            kwargs['oil_type'] = kwargs['oiltype']
            del kwargs['oiltype']

        if 'oil_type' in kwargs:
            # OD 1.14.9: use __set_seed_config__ (accepts Config+Ready) not set_config (Config only)
            self.__set_seed_config__('seed:oil_type', kwargs['oil_type'])
            del kwargs['oil_type']
        else:
            logger.info('Oil type not specified, using default: ' +
                        self.get_config('seed:oil_type'))
        self.set_oiltype(self.get_config('seed:oil_type'))

        if self.oil_weathering_model == 'cicese':
            oil_viscosity = self.oiltype.kvis_at_temp(285)
            if 'oil_density' in kwargs:
                oil_density = kwargs['oil_density']
                logger.info('Using viscosity %s of oiltype %s' %
                            (oil_viscosity, self.get_config('seed:oil_type')))
            else:
                oil_density = self.oiltype.density_at_temp(285)
                logger.info('Using density %s and viscosity %s of oiltype %s' %
                            (oil_density, oil_viscosity, self.get_config('seed:oil_type')))
            kwargs['density'] = oil_density
            kwargs['oil_density'] = oil_density
            kwargs['viscosity'] = oil_viscosity

        elif self.oil_weathering_model == 'noaa':
            oil_density = self.oiltype.density_at_temp(285)
            oil_viscosity = self.oiltype.kvis_at_temp(285)
            # OD 1.14.9: oil_weathering_noaa() uses self.Density / self.KinematicViscosity
            # objects (adios_db style). Set them here since CICOIL skips OpenOil.seed_elements.
            from adios_db.computation.physical_properties import KinematicViscosity, Density
            self.Density = Density(self.oiltype.oil)
            self.KinematicViscosity = KinematicViscosity(self.oiltype.oil)
            logger.info(
                'Using density %s and viscosity %s of oiltype %s' %
                (oil_density, oil_viscosity, self.get_config('seed:oil_type')))
            kwargs['density'] = oil_density
            kwargs['viscosity'] = oil_viscosity

        if 'm3_per_hour' in kwargs:
            m3_per_hour = kwargs['m3_per_hour']
            del kwargs['m3_per_hour']
        else:
            m3_per_hour = self.get_config('seed:m3_per_hour')

        if 'number' in kwargs:
            num_elements = kwargs['number']
        else:
            num_elements = self.get_config('seed:number')
        time = kwargs['time']
        if type(time) is list:
            duration_hours = ((time[1] - time[0]).total_seconds()) / 3600
            if duration_hours == 0:
                duration_hours = 1.
        else:
            duration_hours = 1.  # For instantaneous spill, we use 1h
        kwargs['mass_oil'] = (m3_per_hour * duration_hours / num_elements *
                              kwargs['density'])

        if 'is_bubble' in kwargs:
            if kwargs['is_bubble'] == 1:
                kwargs['mass_gas'] = (m3_per_hour * duration_hours / num_elements * kwargs['density'])
                kwargs['mass_oil'] = 0.000001
        if 'is_dissolved' in kwargs:
            if kwargs['is_dissolved'] == 1:
                kwargs['mass_dissolved_subsea'] = (m3_per_hour * duration_hours / num_elements * kwargs['density'])
                kwargs['mass_oil'] = 0.000001

        super(OpenOil, self).seed_elements(*args, **kwargs)


    def seed_plume_elements(self, lon, lat, plume, seed_time, *args, **kwargs):
        """Seed a given number of particles around given position(s).

        Arguments:
            lon, lat: longitude, latitude of the blowout. They must be float numbers.
            plume: A Tamoc_Plume class object containing all the required
                particle data, including location, depth, exit-from-plume time
                and particle properties.
            seed_time: datetime list, the start and end time of the blowout release.
            z_uncertainty: each size class particles are released around the
                tamoc particle depth z from a normal distribution with 2s=z_uncertainty (m).
            number: integer, total number of particles to be seeded
                Elements are spread around the tamoc-particle location
                with a radius equal to the final tamoc plume half width.
            kwargs: keyword arguments containing properties/attributes and
                values corresponding to the actual particle type (ElementType).
                These are forwarded to the ElementType class. All properties
                for which there are no default value must be specified.
        """
        from opendrift.models.openoil.tamoc_plume import Plume
        import pyproj  # direct import — basereader no longer re-exports pyproj in OD 1.14.9

        if 'number' in kwargs:
            number = kwargs['number']
            del kwargs['number']
        else:
            number = 5000

        if 'z_uncertainty' in kwargs:
            z_uncertainty = kwargs['z_uncertainty']
            del kwargs['z_uncertainty']
        else:
            z_uncertainty = 0.

        if not isinstance(plume, Plume):
            raise ValueError('A tamoc_plume object must be specified')
        if not isinstance(seed_time[0], datetime):
            raise ValueError('start_time must be a datetime object')
        if not isinstance(seed_time[1], datetime):
            raise ValueError('end_time must be a datetime object')
        if not isinstance(lon, float):
            raise ValueError('longitute must be a single float')
        if not isinstance(lat, float):
            raise ValueError('latitude must be a single float')

        duration = seed_time[1] - seed_time[0]

        mass_total_droplets = 0.
        mass_total_bubbles = 0.
        if plume.tracked:
            mass_total_dissolved = np.sum(plume.kg_s_diss_plume + plume.kg_s_diss_farfield)
        else:
            mass_total_dissolved = np.sum(plume.cps)
        for i in range(plume.nparticles):
            if plume.fp_type[i] == 0:
                mass_total_bubbles += np.sum(plume.M_p[i])
            else:
                mass_total_droplets += np.sum(plume.M_p[i])
        total_elements = 0
        # OD 1.14.9: get_environment → self.env.get_environment.
        # env may not be finalized at seed time → use fallback temperature if so.
        try:
            _et, _, _ = self.env.get_environment(['sea_water_temperature'], seed_time[1],
                                                 np.array([lon]), np.array([lat]), np.array([0.0]))
            seawater_temp = np.float64(_et['sea_water_temperature'][0]) + 273.15
        except (AssertionError, AttributeError):
            # Fallback: use configured or canonical Gulf of Mexico surface temperature
            _T_fallback = self.get_config('environment:fallback:sea_water_temperature')
            seawater_temp = np.float64(_T_fallback if _T_fallback else 28.0) + 273.15
            logger.debug(f'seed_plume_elements: env not finalized, using T={seawater_temp:.1f} K')
        rho_p = props.Density(self.fluid_properties.record).at_temp(seawater_temp, 'K')
        oil_viscosity = props.KinematicViscosity(self.fluid_properties.record).at_temp(seawater_temp, temp_units='K')

        # residue_mid_Tb = self.riazi_Tb_ini * (1. + np.exp(self.riazi_intercept) * \
        #                                       np.power(np.log(2.), self.riazi_slope))
        seeded_mass = 0.
        if mass_total_dissolved > 0.:
            dissolved_particle_no = int(np.ceil(duration.total_seconds() / 3600))  # One dissolved particle per hour
            dissolved_particle_mass = mass_total_dissolved * duration.total_seconds() / dissolved_particle_no
            z = -np.random.normal(plume.z[-1], z_uncertainty, dissolved_particle_no)
            z[np.where(z >= 0)] = 0.
            kwargs['z'] = z
            kwargs['is_dissolved'] = 1
            # kwargs['moving'] = 0
            kwargs['m3_per_hour'] = mass_total_dissolved * 3600. / rho_p
            longit, latit = lon, lat
            geod = pyproj.Geod(ellps='WGS84')

            x = np.random.normal(plume.x[-1], plume.b, dissolved_particle_no)
            y = np.random.normal(plume.y[-1], plume.b, dissolved_particle_no)

            az = np.degrees(np.arctan2(x, y))
            ones = np.ones(dissolved_particle_no)
            longit = longit * ones
            latit = latit * ones
            dist = [np.sqrt(x[i] * x[i] + y[i] * y[i]) for i in range(dissolved_particle_no)]
            longit, latit, az = geod.fwd(longit, latit, az, dist, radians=False)
            if dissolved_particle_no > 1:
                self.seed_elements(longit, latit, time=seed_time, number=dissolved_particle_no, *args, **kwargs)
            else:
                self.seed_elements(longit[0], latit[0], time=seed_time[0], number=1, *args, **kwargs)

        for idx in range(plume.nparticles):
            if np.sum(plume.mp[idx, :]) == 0: continue
            if plume.fp_type[idx] == 0:
                element_no = int(np.ceil(duration.total_seconds() / 3600))  # One bubble per hour
                element_mass = np.sum(plume.M_p[idx]) * duration.total_seconds() / element_no
                kwargs['is_bubble'] = 1
                kwargs['is_dissolved'] = 0
                # kwargs['moving'] = 0
                z = np.zeros(element_no)

            else:
                element_no = int(np.ceil(np.sum(plume.M_p[idx]) * number / mass_total_droplets))
                element_mass = np.sum(plume.M_p[idx]) * duration.total_seconds() / element_no
                kwargs['is_bubble'] = 0
                kwargs['is_dissolved'] = 0
                # kwargs['moving'] = 1
                z = -np.random.normal(plume.zp[idx], z_uncertainty, element_no)
                z[np.where(z >= 0)] = 0.

            # Calculate the initial spillet area A0 for the gravity-viscous phase
            kwargs['m3_per_hour'] = np.sum(plume.M_p[idx]) * 3600. / rho_p
            v_w = 9.8 * 10. ** -7.  # kinematic viscosity of seawater @ 23oC, 35o Salinity
            rho_water = 1024.  # Surface seawater density @ 23oC
            g = 9.81  # gravity
            g_r = (rho_water - rho_p) * g / rho_water  # Reduced gravity
            V0 = element_mass / rho_p  # Initial element volume
            A0 = 3.4 * np.pi * np.power(g_r * V0 ** 5. / v_w ** 2., (1. / 6.))  # Initial element area
            droplet_diameter = 2 * (np.sum(plume.mp[idx, :]) / rho_p * 3. / 4. / np.pi) ** (1. / 3.)
            kwargs['diameter'] = droplet_diameter
            kwargs['spillet_area'] = A0
            kwargs['oil_film_thickness'] = V0 / A0
            kwargs['oil_density'] = rho_p
            kwargs['density'] = rho_p
            kwargs['viscosity'] = oil_viscosity

            longit, latit = lon, lat
            total_elements += element_no

            # seeded_mass += elemenmass_total_bubbles = 0.t_no[idx] * element_mass

            geod = pyproj.Geod(ellps='WGS84')
            seed_time = [seed_time[0] + timedelta(seconds=plume.tp[idx]),
                         seed_time[1] + timedelta(seconds=plume.tp[idx])]

            x = np.random.normal(plume.xp[idx], plume.b, element_no)
            y = np.random.normal(plume.yp[idx], plume.b, element_no)

            az = np.degrees(np.arctan2(x, y))
            ones = np.ones(element_no)
            longit = longit * ones
            latit = latit * ones
            dist = [np.sqrt(x[i] * x[i] + y[i] * y[i]) for i in range(element_no)]
            longit, latit, az = geod.fwd(longit, latit, az, dist, radians=False)

            kwargs['z'] = z
            self.keep_droplet_diameter = False
            if element_no > 1:
                self.seed_elements(longit, latit, time=seed_time, number=element_no, *args, **kwargs)
            else:
                self.seed_elements(longit[0], latit[0], time=seed_time[0], number=1, *args, **kwargs)


    def get_oil_budget(self):
        """Get oil budget for the current simulation

        The oil budget consists of the following categories:

        * surface: the sum of variable mass_oil for all active elements where z = 0
        * submerged: the sum of variable mass_oil for all active elements where z < 0
        * stranded: the sum of variable mass_oil for all elements which are stranded
        * evaporated: the sum of variable mass_evaporated for all elements
        * dispersed: the sum of variable mass_dispersed for all elements

        The sum (total mass) shall equal the mass released. Note that the mass of oil
        is conserved, whereas the volume may change with changes in density and
        water uptake (emulsification). Therefore mass should be used for budgets,
        eventually converted to volume (by dividing on density) in the final step
        before presentation.

        Note that mass_oil is the mass of pure oil. The mass of oil emulsion
        (oil containing entrained water droplets) can be calculated as:

        .. code::

            mass_emulsion = mass_oil / (1 - water_fraction)

        I.e. water_fraction = 0 means pure oil, water_fraction = 0.5 means mixture of
        50% oil and 50% water, and water_fraction = 0.9 (which is maximum)
        means 10% oil and 90% water.
        """

        if self.time_step.days < 0:  # Backwards simulation
            return None

        # OD 1.14.9: get_property() returns xarray DataArrays, not masked arrays.
        # Use .values to convert and np.nan_to_num for safety.
        def _prop(name):
            """Return (values_array, status_array) as numpy, summed over trajectories."""
            val, st = self.get_property(name)
            if hasattr(val, 'values'):
                val = val.values
            if hasattr(st, 'values'):
                st = st.values
            return np.nan_to_num(val), st

        z, dummy   = _prop('z')
        mass_oil, status = _prop('mass_oil')
        mass_dissolved_subsea, _ = _prop('mass_dissolved_subsea')
        mass_gas, _  = _prop('mass_gas')
        density = self.get_property('density')[0]
        density = density.values[0, 0] if hasattr(density, 'values') else float(density[0, 0])

        if 'stranded' not in self.status_categories:
            self.status_categories.append('stranded')
        strand_idx = self.status_categories.index('stranded')

        mass_submerged = np.where(
            (status == strand_idx) | (z == 0.0), 0.0, mass_oil)
        mass_submerged = np.nansum(mass_submerged, axis=1)
        mass_dissolved_subsea = np.nansum(mass_dissolved_subsea, axis=1)
        mass_gas  = np.nansum(mass_gas, axis=1)
        mass_surface = np.where(
            (status == strand_idx) | (z < 0.0), 0.0, mass_oil)
        mass_surface  = np.nansum(mass_surface, axis=1)
        mass_stranded = np.nansum(np.where(status != strand_idx, 0.0, mass_oil), axis=1)

        mass_evaporated, _ = _prop('mass_evaporated')
        mass_evaporated = np.nansum(mass_evaporated, axis=1)
        mass_dispersed, _  = _prop('mass_dispersed')
        mass_dispersed  = np.nansum(mass_dispersed, axis=1)
        mass_biodegraded, _ = _prop('mass_biodegraded')
        mass_biodegraded = np.nansum(mass_biodegraded, axis=1)
        mass_photooxidized, _ = _prop('mass_photooxidized')
        mass_photooxidized = np.nansum(mass_photooxidized, axis=1)

        oil_budget = {
            'oil_density':
            density,
            'mass_dispersed':
            mass_dispersed,
            'mass_submerged':
            mass_submerged,
            'mass_surface':
            mass_surface,
            'mass_stranded':
            mass_stranded,
            'mass_dissolved_subsea':
            mass_dissolved_subsea,
            'mass_gas':
            mass_gas,
            'mass_evaporated':
            mass_evaporated,
            'mass_biodegraded':
            mass_biodegraded,
            'mass_photooxidized':
            mass_photooxidized,
            'mass_total': (mass_dispersed + mass_submerged + mass_surface + mass_stranded +
                           mass_dissolved_subsea + mass_gas + mass_evaporated +
                           mass_biodegraded + mass_photooxidized)
        }

        return oil_budget

    def plot_oil_budget(self,
                        filename=None,
                        ax=None,
                        show_density_viscosity=True,
                        show_wind_and_current=True):

        if self.time_step.days < 0:  # Backwards simulation
            fig = plt.figure(figsize=(10, 6.))
            plt.text(0.1, 0.5, 'Oil weathering deactivated for '
                     'backwards simulations')
            plt.axis('off')
            if filename is not None:
                plt.savefig(filename)
                plt.close()
            else:
                plt.show()
            return

        b = self.get_oil_budget()

        oil_budget = np.row_stack(
            (b['mass_dispersed'], b['mass_submerged'], b['mass_dissolved_subsea'], b['mass_surface'],
             b['mass_stranded'], b['mass_evaporated'], b['mass_gas'], b['mass_biodegraded']))
        oil_density = b['oil_density']

        budget = np.cumsum(oil_budget, axis=0)

        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds() / 3600. for t in time_relative])

        if ax is None:
            # Left axis showing oil mass
            nrows = 1
            if show_density_viscosity is True:
                nrows = nrows + 1
            if show_wind_and_current is True:
                nrows = nrows + 1
            fig, axs = plt.subplots(
                nrows=nrows, ncols=1,
                figsize=(10, 6. + (nrows - 1) * 3))  # Suitable aspect ratio
            #ax1 = fig.add_subplot(nrows=nrows, 1, 1)
            if nrows == 1:
                ax1 = axs
            elif nrows >= 2:
                ax1 = axs[0]
                if show_density_viscosity is True:
                    self.plot_oil_watercontent_and_viscosity(ax=axs[1], show=False)
                if show_wind_and_current is True:
                    self.plot_environment(ax=axs[nrows - 1], show=False)
        else:
            ax1 = ax

        # Hack: make some emply plots since fill_between does not support label
        if np.sum(b['mass_dispersed']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='darkslategrey',
                              label='dispersed'))
            ax1.fill_between(time, 0, budget[0, :], facecolor='darkslategrey')
        if np.sum(b['mass_submerged']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='darkblue',
                              label='submerged'))
            ax1.fill_between(time,
                             budget[0, :],
                             budget[1, :],
                             facecolor='darkblue')
        if np.sum(b['mass_dissolved_subsea']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='darkviolet',
                              label='dissolved subsea'))
            ax1.fill_between(time,
                             budget[1, :],
                             budget[2, :],
                             facecolor='darkviolet')
        if np.sum(b['mass_surface']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0), 0, 0, color='royalblue',
                              label='surface'))
            ax1.fill_between(time,
                             budget[2, :],
                             budget[3, :],
                             facecolor='royalblue')
        if np.sum(b['mass_stranded']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0), 0, 0, color='black', label='stranded'))
            ax1.fill_between(time,
                             budget[3, :],
                             budget[4, :],
                             facecolor='black')
        if np.sum(b['mass_evaporated']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='skyblue',
                              label='evaporated'))
            ax1.fill_between(time,
                             budget[4, :],
                             budget[5, :],
                             facecolor='skyblue')

        if np.sum(b['mass_gas']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='orange',
                              label='released gas'))
            ax1.fill_between(time,
                             budget[5, :],
                             budget[6, :],
                             facecolor='orange')
        if np.sum(b['mass_biodegraded']) > 0:
            ax1.add_patch(
                plt.Rectangle((0, 0),
                              0,
                              0,
                              color='indigo',
                              label='biodegraded'))
            ax1.fill_between(time,
                             budget[6, :],
                             budget[7, :],
                             facecolor='indigo')

        ax1.set_ylim([0, budget.max()])
        ax1.set_xlim([0, time.max()])
        ax1.set_ylabel('Mass oil  [%s]' %
                       self.elements.variables['mass_oil']['units'])
        ax1.set_xlabel('Time  [hours]')
        # Right axis showing volume
        ax2 = ax1.twinx()
        mass_total = b['mass_total'][-1]
        ax2.set_ylim([0, mass_total / oil_density])
        ax2.set_ylabel('Volume oil [m3]')
        plt.title('%s (%.1f kg/m3) - %s to %s' %
                  (self.get_oil_name(), oil_density,
                   self.start_time.strftime('%Y-%m-%d %H:%M'),
                   self.time.strftime('%Y-%m-%d %H:%M')))
        # Shrink current axis's height by 10% on the bottom
        box = ax1.get_position()
        ax1.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax2.set_position(
            [box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax1.legend(bbox_to_anchor=(0., -0.10, 1., -0.03),
                   loc=1,
                   ncol=6,
                   mode="expand",
                   borderaxespad=0.,
                   fontsize=10)
        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()

    def plot_density(self, value_array, lon_array, lat_array, title=None, label=None, min_value=None, max_value=None,
                     cmapname=None):
        lat_array, lon_array = np.meshgrid(lat_array, lon_array)
        # map, plt, x, y, index_of_first, index_of_last = self.set_up_map(buffer=.2)
        # ax = plt.gcf().gca()
        fig, ax, crs, x, y, index_of_first, index_of_last = self.set_up_map(buffer=.2)

        ax.set_title(title)
        if cmapname == None:
            cmap = colors.LinearSegmentedColormap.from_list('Oilspill', ['#fcedca', '#cd731c', '#41260c'],
                                                        N=255)
        else:
            cmap = mpt.colormaps[cmapname]
        cmap.set_under('w')
        if min_value == None: min_value = 0.1 * np.sum(value_array) / (len(lon_array) * len(lat_array))
        if max_value == None: max_value = 20. * np.sum(value_array) / (len(lon_array) * len(lat_array))
        gcrs = ccrs.PlateCarree()
        pm = ax.pcolormesh(lon_array, lat_array, value_array, transform=gcrs, cmap=cmap,
                            vmin=min_value, vmax=max_value)
        fig.canvas.draw()
        fig.set_tight_layout(True)
        cb = fig.colorbar(pm, orientation='horizontal', pad=.05, aspect=30, shrink=.8)
                          #location='bottom',
                          #size='3%', pad='5%')
        cb.set_label(label)
        cb.set_alpha(1)
        cb.draw_all()
        plt.show()


    def plot_stats_property(self, value_array, lon_array, lat_array, title=None, label=None, min_value=None, max_value=None,
                     cmapname=None):
        lat_array, lon_array = np.meshgrid(lat_array, lon_array)
        corners = [np.nanmin(lon_array), np.nanmax(lon_array), np.nanmin(lat_array), np.nanmax(lat_array)]
        fig, ax, crs, index_of_first, index_of_last = self.create_new_map(corners=corners, buffer=.2)
        ax.set_title(title)
        if cmapname == None:
            cmap = colors.LinearSegmentedColormap.from_list('Oilspill', ['#fcedca', '#cd731c', '#41260c'],
                                                        N=255)
        else:
            cmap = mpt.colormaps[cmapname]
        cmap.set_under('w')
        if min_value == None: min_value = 0.1 * np.sum(value_array) / (len(lon_array) * len(lat_array))
        if max_value == None: max_value = 20. * np.sum(value_array) / (len(lon_array) * len(lat_array))

        #cb = map.colorbar(pm, label=label, location='bottom',
        #                  size='3%', pad='5%')
        gcrs = ccrs.PlateCarree()
        pm = ax.pcolormesh(lon_array, lat_array, value_array, transform=gcrs, cmap=cmap,
                           vmin=min_value, vmax=max_value)
        fig.canvas.draw()
        fig.set_tight_layout(True)
        cb = fig.colorbar(pm, orientation='horizontal', pad=.05, aspect=30, shrink=.8)
        cb.set_label(label)
        cb.set_alpha(1)
        cb.draw_all()
        plt.show()


    def create_new_map(self,
                   corners=None,
                   buffer=.1,
                   lscale=None,
                   fast=False):
        """
        Generate Figure instance on which trajectories are plotted.

        :param hide_landmask: do not plot landmask (default False)
        :type hide_landmask: bool

        provide corners=[lonmin, lonmax, latmin, latmax] for specific map selection
        """
        from opendrift.readers import reader_global_landmask

        # Initialise map
        if corners is not None:  # User provided map corners
            lonmin = corners[0] - buffer * 2
            lonmax = corners[1] + buffer * 2
            latmin = corners[2] - buffer
            latmax = corners[3] + buffer

        crs = ccrs.Mercator()
        if lscale is None:
            lscale = 'auto'

        globe = crs.globe

        meanlat = (latmin + latmax) / 2
        aspect_ratio = float(latmax - latmin) / (float(lonmax - lonmin))
        aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))
        if aspect_ratio > 1:
            fig = plt.figure(figsize=(11. / aspect_ratio, 11.))
        else:
            fig = plt.figure(figsize=(11., 11. * aspect_ratio))

        ax = fig.add_subplot(111, projection=crs)
        ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree(globe=globe))

        ocean_color = 'white'
        land_color = 'gray'

        reader_global_landmask.plot_land(ax, lonmin, latmin, lonmax,
                                         latmax, fast, ocean_color,
                                         land_color, lscale, globe)

        gl = ax.gridlines(ccrs.PlateCarree(globe=globe), draw_labels=True)
        gl.top_labels = None

        fig.canvas.draw()
        fig.set_tight_layout(True)

        index_of_first = None
        index_of_last = None

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except Exception as e:
            logger.debug('Figure zoom not available: %s', e)

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except Exception as e:
            logger.debug('Figure maximize not available: %s', e)

        return fig, ax, crs, index_of_first, index_of_last


    def get_bins(self, bin_size_deg, lat_min, lat_max, lon_min, lon_max, depths=None):
        lon = self.get_property('lon')[0]
        lat = self.get_property('lat')[0]

        lat_array = np.arange(lat_min, lat_max + bin_size_deg, bin_size_deg)
        lon_array = np.arange(lon_min, lon_max + bin_size_deg, bin_size_deg)
        bins = (lon_array, lat_array)

        z = self.get_property('z')[0]
        z_layer = {}
        status = self.get_property('status')[0]
        lon_stranded = lon.copy()
        lat_stranded = lat.copy()
        # try:
        #     strandnum = self.status_categories.index('stranded')
        #     lon_stranded[status != strandnum] = 1000
        #     lat_stranded[status != strandnum] = 1000
        #     contains_stranded = True
        # except ValueError:
        #     lon_stranded[:] = 1000
        #     lat_stranded[:] = 1000
        #     contains_stranded = False
        #
        # try:
        #     strandnum = self.status_categories.index('stranded')
        #     lon[status == strandnum] = 1000
        #     lat[status == strandnum] = 1000
        #     contains_stranded = True
        # except ValueError:
        #     # lon[:] = 1000
        #     # lat[:] = 1000
        #     contains_stranded = False

        de_bin = {}
        p_groups = {}
        used_de_sizes = [1,]

        # used_de_sizes = np.zeros(6)
        # if depths is None: depths = [10, -0.01, -10000]
        # diameter = self.get_property('diameter')[0]
        # de_index = np.array([0, 2, 4, 6, 8, 10], dtype=np.int32)
        # de_sizes = np.unique(diameter[diameter < 1.]).data
        #
        # ## indexes 0.14, 0.24, 0.33, 0.43, 0.52 corresponds to 2.6%, 6.77%, 15.25%, 31.22%, 56.49% of cumulative mass
        # # de_index_position = np.array([0.14, 0.24, 0.33, 0.43, 0.52])
        # # de_index = np.array(np.unique(np.rint(de_index_position * len(de_sizes.data) - 1)), dtype=np.int32)
        # # de_index = de_index[de_index >=0]
        #
        # for i in range(len(depths)):
        #     if i==0:
        #         z_layer[i] = [lon_stranded, lat_stranded]
        #     else:
        #         z_layer[i] = [lon.copy(), lat.copy()]
        #         z_layer[i][0][z >= depths[i-1]] = 1000
        #         z_layer[i][1][z >= depths[i-1]] = 1000
        #         z_layer[i][0][z < depths[i]] = 1000
        #         z_layer[i][1][z < depths[i]] = 1000
        #
        #     for s in range(len(de_index)):
        #         used_de_sizes[s] = de_sizes[de_index[s]]
        #         de_bin[s] = [z_layer[i][0].copy(), z_layer[i][1].copy()]
        #         if s == 0:
        #             de_bin[s][0][diameter > de_sizes[de_index[s + 1]]] = 1000
        #             de_bin[s][1][diameter > de_sizes[de_index[s + 1]]] = 1000
        #         elif s == len(de_index) - 1:
        #             de_bin[s][0][diameter <= de_sizes[de_index[s]]] = 1000
        #             de_bin[s][1][diameter <= de_sizes[de_index[s]]] = 1000
        #         else:
        #             de_bin[s][0][diameter <= de_sizes[de_index[s]]] = 1000
        #             de_bin[s][1][diameter <= de_sizes[de_index[s]]] = 1000
        #             de_bin[s][0][diameter > de_sizes[de_index[s + 1]]] = 1000
        #             de_bin[s][1][diameter > de_sizes[de_index[s + 1]]] = 1000
        #         p_groups[i * len(de_index) + s] = de_bin[s]

        # fake_de_sizes = np.ones(self.get_property('ID')[0].shape)
        if depths is None: depths = [0, 1]
        for p in range(len(used_de_sizes)):
            for z in range(len(depths)):
                if z == 0:
                    z_layer[z] = [lon_stranded, lat_stranded]
                else:
                    z_layer[z] = [lon, lat]
                    # for s in range(1):
                    #     used_de_sizes[s] = 1
                    p_groups[z] = z_layer[z]

        return bins, p_groups, depths, used_de_sizes

    def get_property_binned(self, bins, p_groups, depths, de_sizes, property='particles_no'):

        times = self.get_time_array()[0]
        weight_array = np.ones(self.get_property('ID')[0].shape)
        if property != 'particles_no':
            weight_array = self.get_property(property)[0]

        H = np.zeros((len(depths), len(de_sizes), len(bins[0]) - 1, len(bins[1]) - 1, len(times)))
        for i in range(len(depths)):
            for s in range(len(de_sizes)):
                for t in range(len(times)):
                    weights = weight_array[t, :]
                    H[i, s, :, :, t], dummy, dummy = \
                        np.histogram2d(p_groups[1][0][t, :], p_groups[1][1][t, :],
                        # np.histogram2d(p_groups[i * len(de_sizes) + s][0][t, :], p_groups[i * len(de_sizes) + s][1][t, :],
                                       weights=weights, bins=bins)
        return H

    def get_unique_visits(self, bins, p_groups, depths, de_sizes):

        times = self.get_time_array()[0]
        weight_visit = np.ones(self.get_property('status')[0].shape)
        weight_mass = self.get_property('mass_oil')[0]
        traj_no = len(weight_visit[0, :])
        Hvisit = np.zeros((len(depths), len(de_sizes), len(bins[0]) - 1, len(bins[1]) - 1, traj_no))
        Hmass = np.zeros((len(depths), len(de_sizes), len(bins[0]) - 1, len(bins[1]) - 1, traj_no))

        for i in range(len(depths)):
            for s in range(len(de_sizes)):
                for j in range(traj_no):
                    Hvisit[i, s, :, :, j], dummy, dummy = np.histogram2d(p_groups[1][0][:, j], p_groups[1][1][:, j],
                                                                         weights=weight_visit[:, j], bins=bins)
                        # np.histogram2d(p_groups[i * len(de_sizes) + s][0][t, :], p_groups[i * len(de_sizes) + s][1][t, :],
                    Hmass[i, s, :, :, j], dummy, dummy = np.histogram2d(p_groups[1][0][:, j], p_groups[1][1][:, j],
                                                                         weights=weight_mass[:, j], bins=bins)

                    Hmass[i, s, :, :, j] = np.where(Hvisit[i, s, :, :, j] > 1., Hmass[i, s, :, :, j] /
                                                    Hvisit[i, s, :, :, j], Hmass[i, s, :, :, j])
                    # Hvisit[i, s, :, :, j] = np.where(Hvisit[i, s, :, :, j] > 1., Hvisit[i, s, :, :, j] /
                    #                                  Hvisit[i, s, :, :, j], Hvisit[i, s, :, :, j])
                    Hvisit[i, s, :, :, j][Hvisit[i, s, :, :, j]!=0.]=1.
        visits = np.sum(np.sum(Hvisit, axis=1), axis=3)
        mass = np.sum(np.sum(Hmass, axis=1), axis=3)
        return visits[0], mass[0], traj_no

    def get_density_array_framed(self, bin_size_deg, lat_min, lat_max, lon_min, lon_max, weight=None):
        lon = self.get_property('lon')[0]
        lat = self.get_property('lat')[0]
        times = self.get_time_array()[0]

        lat_array = np.arange(lat_min - bin_size_deg,
                              lat_max + bin_size_deg, bin_size_deg)
        lon_array = np.arange(lon_min - bin_size_deg,
                              lon_max + bin_size_deg, bin_size_deg)
        bins = (lon_array, lat_array)
        z = self.get_property('z')[0]
        if weight is not None:
            weight_array = self.get_property(weight)[0]

        status = self.get_property('status')[0]
        lon_submerged = lon.copy()
        lat_submerged = lat.copy()
        lon_stranded = lon.copy()
        lat_stranded = lat.copy()
        lon_submerged[z >= 0] = 1000
        lat_submerged[z >= 0] = 1000
        lon[z < 0] = 1000
        lat[z < 0] = 1000
        H = np.zeros((len(times), len(lon_array) - 1,
                      len(lat_array) - 1))  # .astype(int)
        H_submerged = H.copy()
        H_stranded = H.copy()
        try:
            strandnum = self.status_categories.index('stranded')
            lon_stranded[status != strandnum] = 1000
            lat_stranded[status != strandnum] = 1000
            contains_stranded = True
        except ValueError:
            contains_stranded = False

        for i in range(len(times)):
            if weight is not None:
                weights = weight_array[i, :]
            else:
                weights = None
            H[i, :, :], dummy, dummy = \
                np.histogram2d(lon[i, :], lat[i, :],
                               weights=weights, bins=bins)
            H_submerged[i, :, :], dummy, dummy = \
                np.histogram2d(lon_submerged[i, :], lat_submerged[i, :],
                               weights=weights, bins=bins)
            if contains_stranded is True:
                H_stranded[i, :, :], dummy, dummy = \
                    np.histogram2d(lon_stranded[i, :], lat_stranded[i, :],
                                   weights=weights, bins=bins)

        return H, H_submerged, H_stranded, lon_array, lat_array