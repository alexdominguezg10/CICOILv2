"""
Per-component evaporated mass of OpenCiceseOil (OpenDrift 1.14.9) and the volatility bound of the Maya pseudo-component table.

Checks (MAYA, surface oil, evaporation only; a 1 cm film and a thin spreading film):
  1. the per-element sum of comp_*_evaporated equals mass_evaporated (the variables used to hold the remaining mass)
  2. comp_* + comp_*_evaporated equals the initial component mass of each element (nothing else removes oil here)
  3. the cumulative evaporated mass of every component never decreases
  4. the non-volatile residue pseudo-component (Tb about 906 K) evaporates less than 0.1 % of its mass
  5. the total evaporated fraction is below the bound 1 - residue mass fraction (all volatile groups gone), and the thin film evaporates more than the 1 cm film

Run:  python test_evaporation_components.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy)
"""
import sys, logging, warnings, os
from datetime import datetime, timedelta
import numpy as np
logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
if os.environ.get('CICOIL_FILE'):     # test a patched ciceseoil.py without installing it
    import importlib.util, opendrift.models.openoil as _pkg
    _spec = importlib.util.spec_from_file_location('opendrift.models.openoil.ciceseoil', os.environ['CICOIL_FILE'])
    _mod = importlib.util.module_from_spec(_spec); sys.modules[_spec.name] = _mod; _spec.loader.exec_module(_mod)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_constant import Reader as ConstReader

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)


def run(spreading, days):
    o = OpenCiceseOil(weathering_model='cicese', loglevel=50)
    o.set_oiltype('MAYA')
    o.add_reader(ConstReader({
        'x_sea_water_velocity': 0.10, 'y_sea_water_velocity': 0.0, 'x_wind': 5.0, 'y_wind': 0.0,
        'sea_water_temperature': 302.0, 'sea_water_salinity': 36.0, 'sea_floor_depth_below_sea_level': 50.0,
        'sea_surface_height': 0.0, 'ocean_vertical_diffusivity': 0.01, 'sea_surface_wave_significant_height': 1.0,
        'sea_surface_wave_stokes_drift_x_velocity': 0.02, 'sea_surface_wave_stokes_drift_y_velocity': 0.0,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 6.0,
        'sea_surface_wave_mean_period_from_variance_spectral_density_second_frequency_moment': 5.0}))
    o.set_config('drift:vertical_mixing', False)
    for k, v in (('evaporation', True), ('spreading', spreading), ('emulsification', False), ('dispersion', False),
                 ('biodegradation', False), ('handle_released_gas', False), ('subsea_dissolution', False)):
        o.set_config('processes:' + k, v)
    o.seed_elements(lon=-92.1, lat=19.7, number=30, time=datetime(2023, 7, 6, 12), m3_per_hour=100., oiltype='MAYA')
    names = list(o.fluid_properties.composition[9:])
    o.run(duration=timedelta(days=days), time_step=timedelta(hours=1), time_step_output=timedelta(hours=6),
          export_variables=['time', 'lon', 'lat', 'z', 'status', 'mass_oil', 'mass_evaporated'] + ['comp_' + n for n in names] + ['comp_%s_evaporated' % n for n in names])
    return o, names


res = {}
for label, spreading, days in (('1 cm film', False, 3), ('spreading (thin film)', True, 10)):
    print(label)
    o, names = run(spreading, days)
    H = o.result                                       # the comp_* variables are exported explicitly in run()
    g = lambda n: H[n].values.astype('f8')
    ev = np.stack([g('comp_%s_evaporated' % n) for n in names], -1)   # (element, time, component)
    rem = np.stack([g('comp_%s' % n) for n in names], -1)
    mev = g('mass_evaporated')
    ok = np.isfinite(mev) & np.isfinite(ev).all(-1)
    check('%s: sum of comp_*_evaporated = mass_evaporated' % label,
          np.allclose(ev.sum(-1)[ok], mev[ok], rtol=2e-3, atol=1e-3), 'max rel diff %.1e' % np.max(np.abs(ev.sum(-1)[ok] - mev[ok]) / np.maximum(mev[ok], 1e-9)))
    tot = (rem + ev).sum(-1)
    tot0 = np.nanmax(np.where(ok, tot, np.nan), axis=1, keepdims=True)
    rel = ok & (tot > 0)                               # released elements only (before release the arrays are 0)
    T0 = np.broadcast_to(tot0, tot.shape)
    check('%s: remaining + evaporated per element constant in time' % label, np.allclose(tot[rel], T0[rel], rtol=2e-3), 'max rel diff %.1e' % np.max(np.abs(tot[rel] - T0[rel]) / T0[rel]))
    d = np.diff(np.where(ok[..., None], ev, np.nan), axis=1)
    check('%s: cumulative evaporated mass of every component never decreases' % label, np.nanmin(d) > -1e-3, 'min step %.2e kg' % np.nanmin(d))
    ires = names.index('residue') if 'residue' in names else len(names) - 1
    frac_res = ev[..., ires][ok].max() / max((rem + ev)[..., ires][ok].max(), 1e-30)
    check('%s: residue evaporates less than 0.1 %% of its mass' % label, frac_res < 1e-3, '%.2e' % frac_res)
    mf = np.asarray(o.fluid_properties.oil_mass_frac); frac = (ev.sum(-1)[:, -1][ok[:, -1]] / tot0[:, 0][ok[:, -1]]).mean()
    res[label] = frac
    check('%s: evaporated fraction %.3f below the bound %.3f (1 - residue fraction)' % (label, frac, 1 - mf[ires]), frac < 1 - mf[ires])
check('thin film evaporates more than the 1 cm film', res['spreading (thin film)'] > res['1 cm film'], '%.3f against %.3f' % (res['spreading (thin film)'], res['1 cm film']))
print('\nFAILED: ' + ', '.join(FAILS) if FAILS else '\nALL PASSED')
sys.exit(1 if FAILS else 0)
