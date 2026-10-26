"""
Dispersion tests for OpenCiceseOil (OpenDrift 1.14.9).

Calls OpenCiceseOil.dispersion_cicese on a stand-in object (no model run), so the
mass bookkeeping of the Roy dispersion sink is checked directly.

Checks:
  1. the sink acts on surface elements (z >= 0) only: submerged elements keep their
     oil mass and their component vector
  2. oil mass is conserved: mass_oil + mass_dispersed unchanged for every element
  3. processes:dispersion_mode (legacy, transport, budget) switches the Roy sink and the Li entrainment
  4. Roy rate (PyGNOME form): kinematic viscosity under the square root (no division by density),
     inverse scaling with emulsion density and slick thickness, cap at 1, and the mass loss in dispersion_cicese

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
from opendrift.models.openoil.ciceseoil import OpenCiceseOil, roy_dispersed_fraction

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

# 3: dispersion_mode (legacy / transport / budget)
o = OpenCiceseOil(weathering_model='noaa', loglevel=50)
check('default dispersion_mode is legacy', o.get_config('processes:dispersion_mode') == 'legacy')
try:
    o.set_config('processes:dispersion_mode', 'nonsense'); bad = False
except Exception:
    bad = True
check('invalid dispersion_mode is rejected', bad)

parent = OpenCiceseOil.__mro__[1]
calls = []
_orig_swm = parent.surface_wave_mixing
parent.surface_wave_mixing = lambda self, dt: calls.append(dt)
try:
    for mode, disp, sink, li in [('legacy', True, True, True), ('legacy', False, False, True),
                                 ('transport', True, False, True), ('transport', False, False, True),
                                 ('budget', True, True, False), ('budget', False, False, False)]:
        o.set_config('processes:dispersion_mode', mode); o.set_config('processes:dispersion', disp)
        del calls[:]
        o.surface_wave_mixing(900.)
        check(f'{mode:9s} processes:dispersion={disp!s:5s}: Roy sink {"on " if sink else "off"}, Li entrainment {"on " if li else "off"}',
              o._roy_sink_enabled() == sink and (len(calls) == 1) == li)
finally:
    parent.surface_wave_mixing = _orig_swm

# 4: Roy rate
c, dt = 2.0, 900.
f_ref = roy_dispersed_fraction(c, 9e-3, 1000., 1e-3, dt)      # 9,000 cSt; rho * h = 1 kg/m2
expected = 2400. * np.exp(-73.682 * np.sqrt(9e-3)) * c * 3.9e-8 * dt
check('rate for rho * h = 1 kg/m2 is C_Roy * c_disp * V_entrain * dt with C_Roy from sqrt(nu)',
      np.isclose(f_ref, expected, rtol=1e-12), f'{float(f_ref):.3e} vs {expected:.3e}')
f_fresh = roy_dispersed_fraction(c, 4e-4, 1000., 1e-3, dt)    # about 400 cSt
check('viscous emulsion (9,000 cSt) disperses far less than fresh oil (400 cSt)',
      f_ref / f_fresh < 0.01, f'ratio {float(f_ref / f_fresh):.2e}')
check('rate scales as 1 / thickness', np.isclose(roy_dispersed_fraction(c, 9e-3, 1000., 5e-4, dt), 2 * f_ref))
check('rate scales as 1 / emulsion density', np.isclose(roy_dispersed_fraction(c, 9e-3, 2000., 1e-3, dt), 0.5 * f_ref))
check('fraction is capped at 1', float(roy_dispersed_fraction(1e9, 1e-5, 1000., 1e-6, dt)) == 1.0)
check('undefined thickness (0) falls back to 1 mm', np.isclose(roy_dispersed_fraction(c, 9e-3, 1000., 0., dt), f_ref))

m = fake_model(z=[0., 0.], viscosity=4e-4, density=1000., film=1e-3)
OpenCiceseOil.dispersion_cicese(m)
f_exp = float(roy_dispersed_fraction(200. ** 0.57 * 0.1, 4e-4, 1000., 1e-3, 900.))
check('dispersion_cicese removes fraction * mass_oil', np.allclose(m.elements.mass_dispersed, 100. * f_exp),
      f'{m.elements.mass_dispersed[0]:.4f} vs {100. * f_exp:.4f}')

print('\nFAILED: ' + ', '.join(FAILS) if FAILS else '\nALL PASSED')
sys.exit(1 if FAILS else 0)
