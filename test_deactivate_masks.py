#!/usr/bin/env python3
"""deactivate_elements() masks in OpenCiceseOil (OpenDrift 1.14.9): gas bubbles and depleted dissolved elements.

OpenDrift 1.14.9 deactivate_elements() computes `indices & (status == 0)`, so it needs a boolean mask over the active elements.
b032140 passed integer indices (np.where) in handle_gas_particles() and handle_subsea_dissolution(): with 4 bubbles among 19
elements the first time step stopped with "operands could not be broadcast together with shapes (4,) (19,)".
Checks (constant forcing, a few elements, 2 h):
  1. oil elements plus is_bubble elements, processes:handle_released_gas True: the run completes and the bubbles end 'released_gas'
  2. oil elements plus is_dissolved elements with biodegradation (depletion threshold forced high so they deplete): the run
     completes and the dissolved elements end 'dissolved_biodegraded'
Run:  python test_deactivate_masks.py   (CICOIL_FILE=/path/ciceseoil.py to test another copy without installing it)
"""
import sys, os, logging, warnings
from datetime import datetime, timedelta
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
if os.environ.get('CICOIL_FILE'):
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_constant import Reader as ConstReader

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)

def model(**cfg):
    o = OpenCiceseOil(loglevel=50, weathering_model='cicese', seed=1)
    o.add_reader(ConstReader({'x_sea_water_velocity': 0.1, 'y_sea_water_velocity': 0.0, 'x_wind': 5.0, 'y_wind': 0.0,
                              'sea_water_temperature': 28.0, 'sea_water_salinity': 36.0, 'sea_floor_depth_below_sea_level': 50.0,
                              'sea_surface_height': 0.0, 'ocean_vertical_diffusivity': 0.01,
                              'sea_surface_wave_significant_height': 1.0, 'sea_surface_wave_stokes_drift_x_velocity': 0.02,
                              'sea_surface_wave_stokes_drift_y_velocity': 0.0,
                              'sea_surface_wave_period_at_variance_spectral_density_maximum': 6.0,
                              'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 5.0}))
    for k, v in cfg.items():
        o.set_config(k, v)
    o.set_oiltype('MAYA')
    return o

def run(o):
    o.run(duration=timedelta(hours=2), time_step=timedelta(minutes=30), time_step_output=timedelta(hours=1))
    return o.result['status'].values, list(o.status_categories)      # (trajectory, time); a deactivated element keeps its code at the
                                                                    # step it was deactivated and is NaN afterwards

t0 = datetime(2026, 2, 6)
print('1. gas bubbles among oil elements')
o = model(**{'processes:handle_released_gas': True})
o.seed_elements(-92.2, 19.27, z=-5.0, number=15, time=t0, m3_per_hour=10., oil_type='MAYA', max_water=0.9, diameter=3e-3)
o.seed_elements(-92.2, 19.27, z=-5.0, number=4, time=t0, m3_per_hour=1., oil_type='MAYA', max_water=0.9, is_bubble=1, diameter=3e-3)
try:
    st, cats = run(o)
    n = int(np.sum((st == cats.index('released_gas')).any(axis=1)))
    check('run with 4 bubbles among 19 elements completes; bubbles are deactivated as released_gas', n == 4, f'{n} released_gas')
except Exception as e:
    check('run with 4 bubbles among 19 elements completes', False, f'{type(e).__name__}: {str(e)[:90]}')

print('2. dissolved elements with biodegradation')
OpenCiceseOil.BIODEG_WATER_DEPLETION_FRAC = 10.0     # test only: depletion threshold above 1 so the dissolved pool counts as depleted after the transfer step
o = model(**{'processes:biodegradation': True, 'processes:subsea_dissolution': True})
o.seed_elements(-92.2, 19.27, z=-5.0, number=15, time=t0, m3_per_hour=10., oil_type='MAYA', max_water=0.9, diameter=3e-3)
o.seed_elements(-92.2, 19.27, z=-5.0, number=4, time=t0, m3_per_hour=1., oil_type='MAYA', max_water=0.9, is_dissolved=1, diameter=3e-3)
try:
    st, cats = run(o)
    n = int(np.sum((st == cats.index('dissolved_biodegraded')).any(axis=1)))
    check('run with 4 dissolved among 19 elements completes; they are deactivated as dissolved_biodegraded', n == 4, f'{n} dissolved_biodegraded')
except Exception as e:
    check('run with 4 dissolved among 19 elements completes', False, f'{type(e).__name__}: {str(e)[:90]}')

print('\n' + ('ALL PASSED' if not FAILS else 'FAILED: ' + ', '.join(FAILS)))
sys.exit(1 if FAILS else 0)
