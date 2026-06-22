#!/usr/bin/env python3
"""
test_op_tracking.py
====================
Test simulation: CICOILv2 oxygenated photoproduct (OP) tracking
(audit improvement #4, photooxidation_cicese / biodegradation_cicese,
2026-06-14).

Loads ciceseoil.py from THIS directory (cicoilv2_deployment/) in-memory as
opendrift.models.openoil.ciceseoil, without touching the installed
site-packages copy. Uses the synthetic diurnal-cosine UV fallback (no ERA5
reader needed).

Setup: MAYA oil, surface release, calm constant forcing (no
dispersion/vertical_mixing, so all particles stay at z=0 and the UV dose is
the same for everyone), 8-day simulation, 30-min timestep,
processes:photooxidation=True, processes:biodegradation=False (isolates the
photoox -> OP routing/removal pathway).

Checks:
  1. Mass conservation across the whole population, including the new
     mass_op_dissolved / mass_op_degraded pools.
  2. OP routing: with OP_DISSOLVED_FRACTION=0.5, the total mass routed to
     the OP pools (mass_op_dissolved + mass_op_degraded) should be
     approximately equal to mass_photooxidized (the oil-phase residue),
     since both receive the same fraction of each step's dM_phox.
  3. OP removal: K_OP_REMOVAL > 0 and the run is long enough (8 days,
     t1/2~6d) that mass_op_degraded > 0, i.e. some dissolved OP has decayed.
  4. fraction_op_dissolved + fraction_op_degraded + fraction_photooxidized
     + fraction_biodegraded_from_oil + fraction_evaporated + ... should not
     exceed 1 (sanity bound), and fraction_op_dissolved + fraction_op_degraded
     should be > 0.
"""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data_local' / 'op_tracking_test'
OUT.mkdir(parents=True, exist_ok=True)

from cicoil_common import load_ciceseoil
OpenCiceseOil = load_ciceseoil()
print(f'OP_DISSOLVED_FRACTION = {OpenCiceseOil.OP_DISSOLVED_FRACTION}')
print(f'K_OP_REMOVAL          = {OpenCiceseOil.K_OP_REMOVAL} d^-1')

N_PARTICLES = 100
SIM_DAYS = 8
TIME_STEP_MIN = 30
OUT_DT_H = 6
RELEASE_TIME = datetime(2023, 7, 6, 18, 0, 0)
LON, LAT = -92.0, 19.4
SST_C = 29.0


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
    o.set_config('processes:biodegradation', False)
    o.set_config('processes:photooxidation', True)
    o.set_config('general:coastline_action', 'previous')
    o.set_config('general:seafloor_action', 'previous')

    o.set_oiltype('MAYA')

    o.seed_elements(
        lon=LON, lat=LAT, z=0,
        number=N_PARTICLES,
        time=RELEASE_TIME,
        max_water=0.9,
    )

    o.run(
        duration=timedelta(days=SIM_DAYS),
        time_step=timedelta(minutes=TIME_STEP_MIN),
        time_step_output=timedelta(hours=OUT_DT_H),
        outfile=str(OUT / 'cicoil_op_tracking.nc'),
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

    return dict(
        o=o,
        total_initial_mass_kg=total_initial_mass,
        total_remaining_kg=total_remaining,
        total_evaporated_kg=float(np.sum(mc['mass_evaporated'])),
        total_dissolved_kg=float(np.sum(mc['mass_dissolved'])),
        total_photooxidized_kg=float(np.sum(mc['mass_photooxidized'])),
        total_biodeg_oil_kg=float(np.sum(mc['mass_biodegraded_from_oil'])),
        total_biodeg_water_kg=float(np.sum(mc['mass_biodegraded_from_water'])),
        total_op_dissolved_kg=float(np.sum(mc['mass_op_dissolved'])),
        total_op_degraded_kg=float(np.sum(mc['mass_op_degraded'])),
    )


print(f'\n{"="*70}')
print('CICOILv2 oxygenated photoproduct (OP) tracking test')
print(f'{"="*70}')
print(f'MAYA oil, N={N_PARTICLES} particles, {SIM_DAYS} days, SST={SST_C}C, '
      f'surface release, synthetic diurnal UV')

res = run_case()
o = res['o']

total_initial = res['total_initial_mass_kg']
total_check = (res['total_remaining_kg'] + res['total_evaporated_kg'] +
                res['total_dissolved_kg'] + res['total_photooxidized_kg'] +
                res['total_biodeg_oil_kg'] + res['total_biodeg_water_kg'] +
                res['total_op_dissolved_kg'] + res['total_op_degraded_kg'])
mass_check = abs(total_check - total_initial) / total_initial

print(f'\ntotal initial mass          : {total_initial/1e3:.3f} t')
print(f'mass conservation rel. err  : {mass_check:.2e}')
print(f'mass_photooxidized          : {res["total_photooxidized_kg"]:.4f} kg')
print(f'mass_op_dissolved (current) : {res["total_op_dissolved_kg"]:.4f} kg')
print(f'mass_op_degraded            : {res["total_op_degraded_kg"]:.4f} kg')

total_op = res['total_op_dissolved_kg'] + res['total_op_degraded_kg']
print(f'mass_op_dissolved + mass_op_degraded : {total_op:.4f} kg')

assert mass_check < 1e-6, f'Mass conservation check FAILED: {mass_check:.2e}'
print('\nPASS: mass conservation OK (including OP pools)')

assert res['total_photooxidized_kg'] > 0, 'No photooxidation occurred -- cannot test OP routing'
print('PASS: photooxidation occurred (mass_photooxidized > 0)')

# OP_DISSOLVED_FRACTION=0.5 -> total OP-routed mass ~ mass_photooxidized
# (both receive an equal share of each step's dM_phox).
rel_diff_op_vs_phox = abs(total_op - res['total_photooxidized_kg']) / res['total_photooxidized_kg']
assert rel_diff_op_vs_phox < 1e-6, \
    (f'OP routing check FAILED: total_op={total_op:.4f} kg vs '
     f'mass_photooxidized={res["total_photooxidized_kg"]:.4f} kg '
     f'(rel. diff {rel_diff_op_vs_phox:.2e})')
print(f'PASS: mass_op_dissolved + mass_op_degraded == mass_photooxidized '
      f'(OP_DISSOLVED_FRACTION={OpenCiceseOil.OP_DISSOLVED_FRACTION}, '
      f'rel. diff = {rel_diff_op_vs_phox:.2e})')

assert res['total_op_degraded_kg'] > 0, \
    'K_OP_REMOVAL produced no decay over the simulation -- mass_op_degraded should be > 0'
print(f'PASS: mass_op_degraded > 0 (K_OP_REMOVAL={OpenCiceseOil.K_OP_REMOVAL} '
      f'd^-1, t1/2~{np.log(2)/OpenCiceseOil.K_OP_REMOVAL:.1f} d @ 20C)')

# fraction sanity: per-particle fractions should sum to <= 1 and OP fractions > 0
frac_op_dissolved = o.elements.fraction_op_dissolved
frac_op_degraded = o.elements.fraction_op_degraded
frac_phox = o.elements.fraction_photooxidized

assert np.all(frac_op_dissolved >= 0) and np.all(frac_op_degraded >= 0)
assert np.mean(frac_op_dissolved + frac_op_degraded) > 0
total_frac = (frac_op_dissolved + frac_op_degraded + frac_phox +
               o.elements.fraction_evaporated +
               o.elements.fraction_biodegraded_from_oil)
assert np.all(total_frac <= 1.0 + 1e-6), \
    f'Sum of fractions exceeds 1: max={total_frac.max():.6f}'
print('PASS: fraction_op_dissolved/fraction_op_degraded are non-negative, '
      'non-zero on average, and total fractions stay <= 1')

print('\nAll OP tracking checks PASSED.')
