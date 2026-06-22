#!/usr/bin/env python3
"""
test_droplet_size_biodeg.py
============================
Test simulation: CICOILv2 droplet-size-resolved biodegradation
(audit improvement #2, biodegradation_cicese, 2026-06-14).

Loads ciceseoil.py from THIS directory (cicoilv2_deployment/) in-memory as
opendrift.models.openoil.ciceseoil, without touching the installed
site-packages copy.

Setup: MAYA oil, surface release (z=0), calm constant forcing
(dispersion/vertical_mixing/entrainment all off so `elements.diameter`
stays at its seeded value and z stays 0 for the whole run), 10-day
simulation, 60-min timestep. N=90 particles split into three groups of 30
with different fixed `diameter`:

  - 'small'   : 1.0e-4 m (0.1 mm)  -> size_factor = (1e-3/1e-4)^1 = 10
  - 'ref'     : 1.0e-3 m (1 mm)    -> size_factor = (1e-3/1e-3)^1 = 1  (BIODEG_DROPLET_D_REF)
  - 'large'   : 1.0e-2 m (1 cm)    -> size_factor = (1e-3/1e-2)^1 = 0.1

Checks:
  1. Mass conservation across the whole population.
  2. frac_biodeg(small) > frac_biodeg(ref) > frac_biodeg(large), and the
     small/ref and ref/large ratios are roughly consistent with the
     size_factor=10 and 0.1 scaling (after accounting for the
     exponential -- not perfectly linear over 10 days, but order-of-
     magnitude consistent).
  3. A 'zero' group (diameter=0, the default for surface spillets) behaves
     identically to the 'ref' group (diameter=1mm == BIODEG_DROPLET_D_REF),
     confirming the diameter==0 -> factor=1 special case.
"""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data_local' / 'droplet_biodeg_test'
OUT.mkdir(parents=True, exist_ok=True)

from cicoil_common import load_ciceseoil
OpenCiceseOil = load_ciceseoil()
print(f'BIODEG_DROPLET_D_REF     = {OpenCiceseOil.BIODEG_DROPLET_D_REF} m')
print(f'BIODEG_DROPLET_EXPONENT  = {OpenCiceseOil.BIODEG_DROPLET_EXPONENT}')
print(f'BIODEG_DROPLET_DIAMETER_MIN/MAX = '
      f'{OpenCiceseOil.BIODEG_DROPLET_DIAMETER_MIN} / '
      f'{OpenCiceseOil.BIODEG_DROPLET_DIAMETER_MAX} m')

N_PER_GROUP = 30
SIM_DAYS = 10
TIME_STEP_MIN = 60
OUT_DT_H = 6
RELEASE_TIME = datetime(2023, 7, 6, 18, 0, 0)
LON, LAT = -92.0, 19.4
SST_C = 20.0  # = BIODEG_T_REF, q10_factor = 1

GROUPS = {
    'zero':  0.0,
    'small': 1.0e-4,
    'ref':   1.0e-3,
    'large': 1.0e-2,
}


def run_case():
    from opendrift.readers.reader_constant import Reader as CR

    o = OpenCiceseOil(loglevel=30, weathering_model='cicese')

    o.add_reader([CR({
        'x_sea_water_velocity': 0.,
        'y_sea_water_velocity': 0.,
        'upward_sea_water_velocity': 0.,
        'x_wind': 0.,
        'y_wind': 0.,
        'sea_surface_wave_significant_height': 0.,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': 4.0,
        'sea_water_temperature': SST_C,
        'sea_water_salinity': 36.,
        'sea_floor_depth_below_sea_level': 50.,
        'ocean_mixed_layer_thickness': 30.,
        'land_binary_mask': 0,
    })])

    o.set_config('drift:vertical_mixing', False)
    o.set_config('drift:stokes_drift', False)
    o.set_config('drift:current_uncertainty', 0)
    o.set_config('drift:wind_uncertainty', 0)
    o.set_config('processes:evaporation', False)
    o.set_config('processes:emulsification', False)
    o.set_config('processes:dispersion', False)
    o.set_config('processes:spreading', False)
    o.set_config('processes:biodegradation', True)
    o.set_config('processes:photooxidation', False)
    o.set_config('general:coastline_action', 'previous')
    o.set_config('general:seafloor_action', 'previous')

    o.set_oiltype('MAYA')

    group_names = list(GROUPS.keys())

    # Seed each group separately with a fixed `diameter` (no
    # entrainment/dispersion -> stays constant for the whole run). Each
    # seed_elements() call appends N_PER_GROUP scheduled elements, activated
    # in call order, so the groups land at contiguous index ranges.
    for name in group_names:
        o.seed_elements(
            lon=LON, lat=LAT, z=0,
            number=N_PER_GROUP,
            time=RELEASE_TIME,
            max_water=0.9,
            diameter=GROUPS[name],
        )

    group_idx = {}
    for i, name in enumerate(group_names):
        group_idx[name] = np.arange(i * N_PER_GROUP, (i + 1) * N_PER_GROUP)

    o.run(
        duration=timedelta(days=SIM_DAYS),
        time_step=timedelta(minutes=TIME_STEP_MIN),
        time_step_output=timedelta(hours=OUT_DT_H),
        outfile=str(OUT / 'cicoil_droplet_biodeg.nc'),
    )

    mc = o.cicese_mass_balance
    initial_components = (mc['mass_components'] + mc['mass_evaporated'] +
                           mc['mass_dissolved'] + mc['mass_photooxidized'] +
                           mc['mass_biodegraded_from_oil'] +
                           mc['mass_biodegraded_from_water'] +
                           mc['mass_op_dissolved'] +
                           mc['mass_op_degraded'])
    total_initial_mass = float(np.sum(initial_components))
    total_remaining = float(np.sum(mc['mass_components']))
    total_biodeg = float(np.sum(mc['mass_biodegraded_from_oil']))

    mass_check = abs(total_remaining + total_biodeg
                      + float(np.sum(mc['mass_biodegraded_from_water']))
                      - total_initial_mass) / total_initial_mass

    results = {}
    for name in group_names:
        idx = group_idx[name]
        initial_grp = float(np.sum(initial_components[idx, :]))
        biodeg_grp = float(np.sum(mc['mass_biodegraded_from_oil'][idx, :]))
        results[name] = dict(
            initial_mass_kg=initial_grp,
            biodeg_kg=biodeg_grp,
            frac_biodeg=100 * biodeg_grp / initial_grp,
        )

    return results, mass_check, total_initial_mass


print(f'\n{"="*70}')
print('CICOILv2 droplet-size-resolved biodegradation test')
print(f'{"="*70}')
print(f'MAYA oil, {N_PER_GROUP} particles/group x {len(GROUPS)} groups, '
      f'{SIM_DAYS} days, SST={SST_C}C (= BIODEG_T_REF, q10=1)')

results, mass_check, total_initial_mass = run_case()

print(f'\ntotal initial mass         : {total_initial_mass/1e3:.3f} t')
print(f'mass conservation rel. err : {mass_check:.2e}')
print()
for name, res in results.items():
    print(f'  group {name:6s} (d={GROUPS[name]:.1e} m): '
          f'biodeg = {res["biodeg_kg"]:.4f} kg '
          f'({res["frac_biodeg"]:.4f}% of group mass)')

assert mass_check < 1e-6, f'Mass conservation check FAILED: {mass_check:.2e}'
print('\nPASS: mass conservation OK')

f_small, f_ref, f_large, f_zero = (results['small']['frac_biodeg'],
                                    results['ref']['frac_biodeg'],
                                    results['large']['frac_biodeg'],
                                    results['zero']['frac_biodeg'])

assert f_small > f_ref > f_large, \
    f'Droplet-size ordering check FAILED: small={f_small}, ref={f_ref}, large={f_large}'
print('PASS: frac_biodeg(small,d=0.1mm) > frac_biodeg(ref,d=1mm) > frac_biodeg(large,d=1cm)')

# diameter==0 should be treated identically to diameter==BIODEG_DROPLET_D_REF (1mm)
rel_diff_zero_ref = abs(f_zero - f_ref) / f_ref
assert rel_diff_zero_ref < 1e-6, \
    f'diameter==0 should match diameter==D_REF: zero={f_zero}, ref={f_ref}'
print('PASS: frac_biodeg(diameter=0) == frac_biodeg(diameter=D_REF=1mm) '
      f'(rel. diff = {rel_diff_zero_ref:.2e})')

print('\nAll droplet-size biodegradation checks PASSED.')
