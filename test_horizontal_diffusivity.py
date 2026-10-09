"""
Horizontal diffusivity test for OpenCiceseOil (OpenDrift 1.14.9).

Regression test for the bug where drift:horizontal_diffusivity had no effect
(the variable was missing from required_variables, so the random walk never ran).

Checks (constant forcing, no current, no wind, transport-only, 6 h):
  1. D = 0 (default)        -> particles do not move
  2. D = 10 m2/s via drift:horizontal_diffusivity -> spread std ~ sqrt(2 D t)
  3. D = 10 m2/s via environment:fallback:horizontal_diffusivity -> same spread
  4. D = 40 spreads ~2x more than D = 10

Run:  python test_horizontal_diffusivity.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy)
"""
import sys, logging, warnings
from datetime import datetime, timedelta
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
import os
if os.environ.get('CICOIL_FILE'):     # test a patched ciceseoil.py without installing it
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_constant import Reader as ConstReader

HOURS, DT_MIN, N = 6, 15, 2000
FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)

def spread_m(D=0.0, via='drift'):
    o = OpenCiceseOil(weathering_model='noaa', loglevel=50)
    o.set_oiltype('MAYA')
    for k in ('evaporation', 'emulsification', 'dispersion', 'spreading'):
        o.set_config(f'processes:{k}', False)
    o.set_config('drift:vertical_mixing', False)
    o.set_config('drift:current_uncertainty', 0.0); o.set_config('drift:wind_uncertainty', 0.0)
    if via == 'drift':
        o.set_config('drift:horizontal_diffusivity', D)
    else:
        o.set_config('environment:fallback:horizontal_diffusivity', D)
    o.add_reader(ConstReader({'x_sea_water_velocity': 0., 'y_sea_water_velocity': 0., 'x_wind': 0.,
        'y_wind': 0., 'sea_water_temperature': 302., 'sea_water_salinity': 36.,
        'sea_floor_depth_below_sea_level': 50., 'sea_surface_height': 0.,
        'ocean_vertical_diffusivity': 0.01, 'sea_surface_wave_significant_height': 1.,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 6.,
        'sea_surface_wave_stokes_drift_x_velocity': 0., 'sea_surface_wave_stokes_drift_y_velocity': 0.,
        'sea_ice_area_fraction': 0., 'land_binary_mask': 0.}))
    o.seed_elements(lon=-92., lat=19., time=datetime(2026, 10, 7), number=N, radius=0,
                    m3_per_hour=1.0, oil_type='MAYA')
    np.random.seed(1)
    o.run(duration=timedelta(hours=HOURS), time_step=timedelta(minutes=DT_MIN), time_step_output=timedelta(hours=HOURS))
    lon, lat = o.result.lon.values[:, -1], o.result.lat.values[:, -1]
    x = (lon + 92.) * 111190. * np.cos(np.deg2rad(19.)); y = (lat - 19.) * 111190.
    return float(np.sqrt(0.5 * (x.std() ** 2 + y.std() ** 2)))

t = HOURS * 3600.
s0 = spread_m(0.0)
check('D=0: no spreading', s0 < 1.0, f'{s0:.2f} m')
s10 = spread_m(10.0)
exp10 = np.sqrt(2 * 10. * t)
check('D=10 (drift:) spreads as sqrt(2Dt)', abs(s10 / exp10 - 1) < 0.1, f'{s10:.0f} m vs {exp10:.0f} m')
s10f = spread_m(10.0, via='fallback')
check('D=10 (environment:fallback) spreads as sqrt(2Dt)', abs(s10f / exp10 - 1) < 0.1, f'{s10f:.0f} m')
s40 = spread_m(40.0)
check('D=40 spreads ~2x D=10', abs(s40 / s10 - 2) < 0.2, f'ratio {s40 / s10:.2f}')
print('\nRESULT:', 'PASS' if not FAILS else f'FAIL ({len(FAILS)})'); sys.exit(1 if FAILS else 0)
