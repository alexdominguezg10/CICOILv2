"""
Thickness read by the evaporation in slick_thickness_model=field (OpenCiceseOil, OpenDrift 1.14.9).

The slick field holds emulsion volume per area. slick_field:evaporation_thickness chooses what the evaporation reads:
  oil      h_cell * V_oil / V_emulsion  (the oil layer; default)       emulsion   h_cell  (behaviour of commit 66bce53)

Set-up: Maya, 24 h, dispersion off, field 1000 m, constant 5 m/s wind, continuous release. The emulsification process is off and the water fraction of
every element is set to 0.6 before each field step (an emulsion does not form in 24 h under this constant forcing), so that the test acts on the conversion only.
Checks:
  1. default is oil; an invalid value is rejected
  2. emulsion: spillet_thickness equals the cell thickness (sampled from the field) for every surface element
  3. oil: spillet_thickness equals the cell thickness times V_oil / V_emulsion, and spillet_area = V_emulsion / cell thickness
  4. the water fraction reaches > 0.5 (so the test is not trivial) and oil thickness / emulsion thickness is about 1 - w there
  5. an instantaneous release (no scheduled elements) runs with the field and gives a finite thickness (it used to fail in _slick_field_init)
  6. same release, same field: the oil option evaporates more than the emulsion option, and both stay below the Maya bound 1 - residue
Run:  python test_slick_lens_thickness.py   (CICOIL_FILE=/path/ciceseoil.py to test a patched copy; slick_field.py next to it)
"""
import sys, logging, warnings, os
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

REC = []
_orig = OpenCiceseOil._slick_field_step
W_FORCED = 0.6
def _rec(self):
    self.elements.water_fraction[:] = W_FORCED
    _orig(self)
    st = self._slick; el = self.elements; sf = st['sf']; g = st['grid']
    surf = np.where(el.z == 0)[0]
    if len(surf) == 0: return
    vol = sf.emulsion_volume(el.mass_oil[surf], el.water_fraction[surf], el.oil_density[surf], np.asarray(self.local_sea_water_density())[surf])
    i = np.floor((el.lon[surf] - g.lon0) * g.kx / g.dx).astype(int); j = np.floor((el.lat[surf] - g.lat0) * g.ky / g.dx).astype(int)
    ins = (i >= 0) & (i < g.nx) & (j >= 0) & (j < g.ny)
    REC.append(dict(h=sf.sample(g, st['h'], el.lon[surf], el.lat[surf], floor=1e-9)[ins], thick=el.spillet_thickness[surf][ins].astype('f8'),
                    area=el.spillet_area[surf][ins].astype('f8'), voil=(el.mass_oil[surf] / el.oil_density[surf])[ins].astype('f8'),
                    vem=vol[ins].astype('f8'), w=el.water_fraction[surf][ins].astype('f8')))
OpenCiceseOil._slick_field_step = _rec


def run(mode, instant=False):
    REC.clear()
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
    for k, v in (('evaporation', True), ('spreading', False), ('emulsification', False), ('dispersion', False), ('biodegradation', False),
                 ('handle_released_gas', False), ('subsea_dissolution', False), ('slick_thickness_model', 'field')):
        o.set_config('processes:' + k, v)
    o.set_config('slick_field:dx', 1000.)
    if mode is not None: o.set_config('slick_field:evaporation_thickness', mode)
    o.seed_elements(lon=-92.1, lat=19.7, number=30, time=(datetime(2023, 7, 6, 12) if instant else [datetime(2023, 7, 6, 12), datetime(2023, 7, 6, 18)]), m3_per_hour=100., oiltype='MAYA')
    o.run(duration=timedelta(hours=24), time_step=timedelta(hours=1), time_step_output=timedelta(hours=6), export_variables=['time', 'lon', 'lat', 'z', 'status', 'mass_oil', 'mass_evaporated'])
    H = o.result
    mo = H['mass_oil'].values.astype('f8'); me = H['mass_evaporated'].values.astype('f8')
    ok = np.isfinite(me) & np.isfinite(mo)
    frac = (me[:, -1][ok[:, -1]] / np.maximum((mo + me)[:, -1][ok[:, -1]], 1e-30)).mean()
    mf = np.asarray(o.fluid_properties.oil_mass_frac)
    return o, [dict(r) for r in REC], frac, 1 - mf[-1]


o = OpenCiceseOil(weathering_model='cicese', loglevel=50)
check('evaporation_thickness defaults to oil', o.get_config('slick_field:evaporation_thickness') == 'oil')
try: o.set_config('slick_field:evaporation_thickness', 'nonsense'); bad = False
except Exception: bad = True
check('invalid value rejected', bad)

_, rec_e, f_e, bound = run('emulsion')
d = np.concatenate([r['thick'] / r['h'] for r in rec_e])
check('emulsion: spillet_thickness = cell thickness', np.allclose(d, 1.0, rtol=1e-5), 'max |ratio-1| %.1e' % np.max(np.abs(d - 1)))

_, rec_o, f_o, _ = run(None)                                  # default = oil
cat = lambda k: np.concatenate([r[k] for r in rec_o])
exp = cat('h') * cat('voil') / cat('vem')
check('oil: spillet_thickness = cell thickness * V_oil / V_emulsion', np.allclose(cat('thick'), exp, rtol=1e-4), 'max rel diff %.1e' % np.max(np.abs(cat('thick') / exp - 1)))
check('oil: spillet_area = V_emulsion / cell thickness', np.allclose(cat('area'), cat('vem') / cat('h'), rtol=1e-4), 'max rel diff %.1e' % np.max(np.abs(cat('area') / (cat('vem') / cat('h')) - 1)))
w = cat('w'); big = w > 0.5
check('water fraction exceeds 0.5 (test is not trivial)', big.any() and np.allclose(w, W_FORCED), 'max %.2f' % w.max())
if big.any():
    rat = cat('thick')[big] / cat('h')[big]
    check('oil / emulsion thickness is about 1 - w where w > 0.5', np.allclose(rat, 1 - w[big], atol=0.05), 'ratio %.2f..%.2f, 1-w %.2f..%.2f' % (rat.min(), rat.max(), (1 - w[big]).max(), (1 - w[big]).min()))
_, rec_i, f_i, _ = run(None, instant=True)
th = np.concatenate([r['thick'] for r in rec_i])
check('instantaneous release runs with the field and gives a finite, positive thickness', len(th) > 0 and np.isfinite(th).all() and (th > 0).all(), '%d samples, evaporated %.4f' % (len(th), f_i))
check('oil evaporates more than emulsion after 24 h', f_o > f_e, '%.4f against %.4f' % (f_o, f_e))
check('both below the Maya bound %.3f' % bound, max(f_o, f_e) < bound)
print('\nFAILED: ' + ', '.join(FAILS) if FAILS else '\nALL PASSED')
sys.exit(1 if FAILS else 0)
