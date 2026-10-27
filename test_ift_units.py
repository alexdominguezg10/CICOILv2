"""
Interfacial-tension units test for OpenCiceseOil (OpenDrift 1.14.9).

Regression test: the element field IFT defaulted to 30 (mN/m) but the Li et al. (2017) entrainment rate and the
Johansen droplet spectrum take N/m, so the cicese branch used sigma 1000 times too large (characteristic diameter
d_o 0.72 m instead of 0.023 m, entrainment rate lower by about 2.5).

Checks:
  1. the IFT default is 0.03 N/m
  2. the Li rate with the default IFT equals an independent hand calculation (rho_o 930, rho_w 1025, Hs 1.5, f 0.03, mu 0.05)
  3. the rate with the old default (30) is lower by a factor 2 to 3 (documents the size of the error)
  4. no `* 1000.` is left on the oil-water surface tension in ciceseoil.py

Run:  python test_ift_units.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy)
"""
import sys, logging, warnings, os, re, math
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
if os.environ.get('CICOIL_FILE'):
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
    SRC = os.environ['CICOIL_FILE']
else:
    import opendrift.models.openoil.ciceseoil as _mod; SRC = _mod.__file__
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.models.physics_methods import oil_wave_entrainment_rate_li2017

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)

def hand(sig, mu, ro=930., rw=1025., Hs=1.5, f=0.03, g=9.81):
    do = 4 * math.sqrt(sig / ((rw - ro) * g)); we = rw * g * Hs * do / sig; oh = mu / math.sqrt(ro * sig * do)
    return 4.604e-10 * we ** 1.805 * oh ** -1.023 * f

dflt = float(OpenCiceseOil.ElementType.variables['IFT']['default'])
check('IFT default is 0.03 N/m', abs(dflt - 0.03) < 1e-9, f'{dflt}')
kw = dict(dynamic_viscosity=0.05, oil_density=930., significant_wave_height=1.5, wave_breaking_fraction=0.03, sea_water_density=1025.)
r_new = float(oil_wave_entrainment_rate_li2017(interfacial_tension=dflt, **kw))
check('Li rate with the default IFT matches the hand calculation', abs(r_new / hand(0.03, 0.05) - 1) < 1e-6, f'{r_new:.4e} vs {hand(0.03, 0.05):.4e}')
r_old = float(oil_wave_entrainment_rate_li2017(interfacial_tension=30., **kw))
check('old default (30) gives a rate 2-3 times lower', 2.0 < r_new / r_old < 3.0, f'ratio {r_new / r_old:.2f}')
src = open(SRC).read()
check('no * 1000. on the surface tension left', re.search(r'oil_water_surface_tension\(\)\s*\*\s*1000', src) is None)
print('FAILED: ' + ', '.join(FAILS) if FAILS else 'ALL PASSED'); sys.exit(1 if FAILS else 0)
