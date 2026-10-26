"""
Slick thickness model option tests for OpenCiceseOil (OpenDrift 1.14.9): configuration and guards; the field itself is tested in test_slick_field.py.

Checks:
  1. processes:slick_thickness_model defaults to legacy; field is accepted; an invalid value is rejected; slick_field:dx defaults to 500 m
  2. field with processes:spreading True, or with processes:julia_weathering True, is refused (ValueError)
  3. _linear_weights interpolates a linear function exactly (and clamps outside the axis)

Run:  python test_slick_coupling.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy; slick_field.py must be next to it)
"""
import sys, logging, warnings, os
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
if os.environ.get('CICOIL_FILE'):     # test a patched ciceseoil.py without installing it
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil, _linear_weights

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)


def raises(f, exc=Exception):
    try: f()
    except exc: return True
    return False


o = OpenCiceseOil(weathering_model='cicese', loglevel=50)
check('slick_thickness_model defaults to legacy', o.get_config('processes:slick_thickness_model') == 'legacy')
check('slick_field:dx defaults to 500 m', o.get_config('slick_field:dx') == 500.)
o.set_config('processes:slick_thickness_model', 'field')
check('field is accepted', o.get_config('processes:slick_thickness_model') == 'field')
check('invalid model is rejected', raises(lambda: o.set_config('processes:slick_thickness_model', 'nonsense')))

o.set_config('processes:spreading', True)
check('field with processes:spreading True is refused', raises(o._slick_field_init, ValueError))
o.set_config('processes:spreading', False); o.set_config('processes:julia_weathering', True)
check('field with processes:julia_weathering True is refused', raises(o._slick_field_init, ValueError))

axis = np.array([0., 1., 2.5, 4.]); q = np.array([-1., 0., 0.4, 1., 1.75, 3.9, 5.])
M = _linear_weights(axis, q)
f = lambda x: 3. * x - 2.
check('_linear_weights is exact for a linear function inside the axis', np.allclose((M @ f(axis))[1:6], f(q[1:6])))
check('_linear_weights clamps outside the axis', np.allclose((M @ f(axis))[[0, 6]], [f(0.), f(4.)]))
check('_linear_weights rows sum to 1', np.allclose(M.sum(1), 1.))

print('\nFAILED: ' + ', '.join(FAILS) if FAILS else '\nALL PASSED')
sys.exit(1 if FAILS else 0)
