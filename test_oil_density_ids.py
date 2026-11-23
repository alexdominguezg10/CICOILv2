#!/usr/bin/env python3
"""oil_density() must index cicese_mass_balance by element ID, not by array position.

Staggered seeding (plume classes at different times: gas, dissolved, oil) plus removal of the bubbles in the first step makes
positions differ from IDs. With positional indexing an oil element read the zeroed mass row of a dissolved element (0/0 = NaN):
it stayed status 0 with NaN density and mass ("ghost") and dropped out of the mass budget.
Check: after a 36 h run, no OIL element (not dissolved, not bubble) is active with NaN mass_oil.
Run:  python test_oil_density_ids.py   (CICOIL_FILE=/path/ciceseoil.py to test another copy without installing it)
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
from opendrift.models.openoil.tamoc_plume import Plume
from opendrift.readers.reader_constant import Reader as ConstReader


class P(Plume):
    """Minimal stand-in for a TAMOC plume: 2 gas classes, 3 oil classes, dissolved pool; released at different times."""
    def __init__(self):
        self.fp_type = [0, 0, 1, 1, 1]; self.nparticles = 5; self.tracked = False
        self.xp = np.linspace(0, 2, 5); self.yp = np.linspace(0, 1, 5); self.zp = np.full(5, -0.5); self.b = 2.0
        self.tp = np.array([0., 0., 300., 600., 900.])
        mass = np.array([0.5, 0.4, 5.0, 4.0, 3.0])
        self.M_p = [np.array([m]) for m in mass]
        self.mp = np.array([[m / 50.0] for m in mass])
        self.cps = np.array([0.4]); self.z = np.array([-0.5]); self.x = np.array([0.0]); self.y = np.array([0.0])


o = OpenCiceseOil(loglevel=50, weathering_model='cicese', seed=1)
o.add_reader(ConstReader({'x_sea_water_velocity': 0.1, 'y_sea_water_velocity': 0.0, 'x_wind': 5.0, 'y_wind': 0.0,
                          'sea_water_temperature': 28.0, 'sea_water_salinity': 36.0, 'sea_floor_depth_below_sea_level': 50.0,
                          'sea_surface_height': 0.0, 'ocean_vertical_diffusivity': 0.01,
                          'sea_surface_wave_significant_height': 1.0, 'sea_surface_wave_stokes_drift_x_velocity': 0.02,
                          'sea_surface_wave_stokes_drift_y_velocity': 0.0,
                          'sea_surface_wave_period_at_variance_spectral_density_maximum': 6.0,
                          'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 5.0}))
for k, v in {'drift:vertical_mixing': True, 'processes:handle_released_gas': True, 'processes:subsea_dissolution': True,
             'processes:dispersion': True, 'processes:biodegradation': True}.items():
    o.set_config(k, v)
o.set_oiltype('MAYA'); t0 = datetime(2026, 2, 6)
o.seed_plume_elements(lon=-92.2, lat=19.27, plume=P(), seed_time=[t0, t0 + timedelta(hours=24)], number=60, z_uncertainty=2.0, max_water=0.9)
n = o.num_elements_total()
o.run(duration=timedelta(hours=36), time_step=timedelta(minutes=30), time_step_output=timedelta(hours=1))
r = o.result
m = r['mass_oil'].values; st = r['status'].values
ib = np.nanmax(r['is_bubble'].values, axis=1); idis = np.nanmax(r['is_dissolved'].values, axis=1)
last = m.shape[1] - 1
oil = (ib != 1) & (idis != 1)
ghost_oil = int(((st[:, last] == 0) & ~np.isfinite(m[:, last]) & oil).sum())
print(f'{n} elements, {int(oil.sum())} oil; oil elements active with NaN mass at the end: {ghost_oil}')
ok = ghost_oil == 0 and int(oil.sum()) > 0
print('PASSED' if ok else 'FAILED')
sys.exit(0 if ok else 1)
