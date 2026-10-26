"""
Sea-water density test for OpenCiceseOil (OpenDrift 1.14.9).

Regression test for the bug where sea_water_density() was called without
arguments (emulsion density, Li et al. entrainment, droplet spectrum), which
returns the density for the OpenDrift defaults (T = 10 degC, S = 35) instead of
the water the elements are in.

Checks:
  1. the default call still returns 1026.95 kg/m3 (documents the trap)
  2. local_sea_water_density() at 30.2 degC, S 36.6 returns 1022.86 kg/m3
  3. same result when the environment temperature is already in K
  4. element-wise: two elements with different T, S get different densities
  5. no bare self.sea_water_density() call is left in ciceseoil.py

Run:  python test_local_sea_water_density.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy)
"""
import sys, logging, warnings, os, re
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
if os.environ.get('CICOIL_FILE'):     # test a patched ciceseoil.py without installing it
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.models.physics_methods import PhysicsMethods

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)

class Env:                              # stands in for the environment recarray
    def __init__(self, t, s):
        self.sea_water_temperature = np.asarray(t, dtype='f8')
        self.sea_water_salinity = np.asarray(s, dtype='f8')

o = OpenCiceseOil(weathering_model='noaa', loglevel=50)
ref = PhysicsMethods.sea_water_density(T=30.2, S=36.6)

check('default call is the 10 degC / S 35 value', abs(o.sea_water_density() - 1026.952) < 0.01,
      f'{float(o.sea_water_density()):.3f}')

o.environment = Env([30.2, 30.2], [36.6, 36.6])
d = o.local_sea_water_density()
check('local density at 30.2 degC, S 36.6', np.allclose(d, ref, atol=1e-6) and abs(d[0] - 1022.86) < 0.01,
      f'{d[0]:.3f} vs {ref:.3f}')

o.environment = Env([30.2 + 273.15, 30.2 + 273.15], [36.6, 36.6])
check('same result when the environment temperature is in K', np.allclose(o.local_sea_water_density(), ref, atol=1e-6))

o.environment = Env([30.2, 20.0], [36.6, 35.0])
d = o.local_sea_water_density()
check('element-wise: different T, S give different density', d.shape == (2,) and d[1] > d[0] + 1.0,
      f'{d[0]:.2f}, {d[1]:.2f}')

src = open(os.environ.get('CICOIL_FILE', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ciceseoil.py'))).read()
check('no bare self.sea_water_density() call in ciceseoil.py', not re.search(r'self\.sea_water_density\(\s*\)', src))

print('\nFAILED: ' + ', '.join(FAILS) if FAILS else '\nALL PASSED')
sys.exit(1 if FAILS else 0)
