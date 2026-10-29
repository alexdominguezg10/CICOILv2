"""
Emulsification constants of OpenCiceseOil: the Abkatun-2026 tuning (EMUL_STABILITY_MIN 0.3, EMUL_STABILITY_MAX 3.0,
EMUL_RATE_SCALE 0.0142) was removed on 2026-10-29 (it was tuned to a tar-ball day-17.2 target with no data support).
With all three at 1.0 the viscosity stability clip is 1 for any viscosity and the water-uptake rate is the stock NOAA one.

Checks:
  1. source level (ast, no OpenDrift needed): the three class attributes and the three CICESE_DEFAULTS entries are 1.0
  2. source level: emulsification_cicese still multiplies by stability and EMUL_RATE_SCALE (so 1.0 makes the factor neutral)
  3. class level (only if opendrift is importable): the attributes are 1.0 and np.clip(nu/NU_REF, MIN, MAX) == 1 for 1e-3..1e7 cSt

Run:  python test_emul_constants.py   (CICOIL_FILE=/path/ciceseoil.py to test another copy)
"""
import ast, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FILE = os.environ.get('CICOIL_FILE', os.path.join(HERE, 'ciceseoil.py'))
NAMES = ('EMUL_STABILITY_MIN', 'EMUL_STABILITY_MAX', 'EMUL_RATE_SCALE')
FAILS = []


def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond:
        FAILS.append(name)


src = open(FILE).read()
tree = ast.parse(src)
cls = next(n for n in ast.walk(tree) if isinstance(n, ast.ClassDef) and n.name == 'OpenCiceseOil')
attrs = {t.id: n.value.value for n in cls.body if isinstance(n, ast.Assign) and isinstance(n.value, ast.Constant)
         for t in n.targets if isinstance(t, ast.Name) and t.id in NAMES}
dflt = {}
for n in cls.body:
    if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'CICESE_DEFAULTS' for t in n.targets):
        dflt = {k.value: v.value for k, v in zip(n.value.keys, n.value.values) if isinstance(k, ast.Constant) and k.value in NAMES}
print('file:', FILE)
for k in NAMES:
    check('class attribute %s == 1.0' % k, attrs.get(k) == 1.0, str(attrs.get(k)))
    check('CICESE_DEFAULTS[%s] == 1.0' % k, dflt.get(k) == 1.0, str(dflt.get(k)))
fn = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'emulsification_cicese')
seg = ast.get_source_segment(src, fn)
check('emulsification_cicese applies stability and EMUL_RATE_SCALE', 'k_emul * stability * self.EMUL_RATE_SCALE' in seg)

try:
    from opendrift.models.openoil.ciceseoil import OpenCiceseOil
except ImportError:
    print('  (opendrift not importable here: class-level check skipped)')
else:
    for k in NAMES:
        check('imported OpenCiceseOil.%s == 1.0' % k, getattr(OpenCiceseOil, k) == 1.0, str(getattr(OpenCiceseOil, k)))
    nu = np.logspace(-3, 7, 41)
    s = np.clip(nu / OpenCiceseOil.EMUL_STABILITY_NU_REF, OpenCiceseOil.EMUL_STABILITY_MIN, OpenCiceseOil.EMUL_STABILITY_MAX)
    check('stability factor == 1 for 1e-3..1e7 cSt', bool(np.all(s == 1.0)))

print('FAILED: %s' % FAILS if FAILS else 'OK')
sys.exit(1 if FAILS else 0)
