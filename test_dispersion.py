"""
Dispersion tests for OpenCiceseOil (OpenDrift 1.14.9).

Calls OpenCiceseOil.dispersion_cicese on a stand-in object (no model run), so the
mass bookkeeping of the Roy dispersion sink is checked directly.

Checks:
  1. the sink acts on surface elements (z >= 0) only: submerged elements keep their
     oil mass and their component vector
  2. oil mass is conserved: mass_oil + mass_dispersed unchanged for every element

Run:  python test_dispersion.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy)
"""
import sys, logging, warnings, os, types
from datetime import timedelta
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
if os.environ.get('CICOIL_FILE'):     # test a patched ciceseoil.py without installing it
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)


def fake_model(z, viscosity=1e-5, density=900., mass=100., film=1e-3, energy=200., fbw=0.1, dt_s=900., water_fraction=0.):
    """Stand-in for self in dispersion_cicese: only the attributes the method uses."""
    n = len(z)
    el = types.SimpleNamespace(
        z=np.asarray(z, float), viscosity=np.full(n, viscosity), density=np.full(n, density),
        mass_oil=np.full(n, mass), mass_dispersed=np.zeros(n), ID=np.arange(n),
        oil_film_thickness=np.full(n, film), water_fraction=np.full(n, water_fraction))
    m = types.SimpleNamespace(
        elements=el, time_step=timedelta(seconds=dt_s),
        wave_energy_dissipation=lambda: np.full(n, energy),
        sea_surface_wave_breaking_fraction=lambda: np.full(n, fbw),
        cicese_mass_balance={'mass_components': np.ones((n, 3)) * mass / 3.},
        get_config=lambda key: None)
    return m


# 1 and 2: surface-only sink and mass conservation
m = fake_model(z=[0., -1., 0., -5.])
comp0 = m.cicese_mass_balance['mass_components'].copy()
OpenCiceseOil.dispersion_cicese(m)
el = m.elements
check('surface elements lose oil mass', el.mass_oil[0] < 100. and el.mass_oil[2] < 100.,
      f'{el.mass_oil[0]:.4f}, {el.mass_oil[2]:.4f}')
check('submerged elements keep their oil mass', el.mass_oil[1] == 100. and el.mass_oil[3] == 100.)
check('submerged elements keep their component vector',
      np.array_equal(m.cicese_mass_balance['mass_components'][[1, 3]], comp0[[1, 3]]))
check('surface elements: component vector scaled with the mass loss',
      np.allclose(m.cicese_mass_balance['mass_components'][0].sum(), el.mass_oil[0]))
check('oil mass conserved (mass_oil + mass_dispersed)', np.allclose(el.mass_oil + el.mass_dispersed, 100.))

print('\nFAILED: ' + ', '.join(FAILS) if FAILS else '\nALL PASSED')
sys.exit(1 if FAILS else 0)
