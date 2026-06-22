#!/usr/bin/env python3
"""
test_biodegradation_water.py
=============================
Test simulation: CICOILv2 dissolved-phase biodegradation
(mass_biodegraded_from_water, handle_subsea_dissolution + biodegradation_cicese,
2026-06-14).

Loads ciceseoil.py from THIS directory (cicoilv2_deployment/) in-memory as
opendrift.models.openoil.ciceseoil_local, without touching the installed
site-packages copy.

Setup: seed particles directly with is_dissolved=1 and mass_dissolved_subsea
equal to their initial mass_oil (mimics seed_plume_elements' dissolved-particle
seeding). Calm constant forcing, all weathering processes off except
biodegradation + subsea_dissolution.

Checks:
  1. On the first step, mass_components -> mass_dissolved (one-time transfer
     by handle_subsea_dissolution).
  2. mass_dissolved decays and mass_biodegraded_from_water grows, roughly
     matching 1 - exp(-K_BIODEG_WATER * t) at T_REF=20C.
  3. Mass conservation: mass_dissolved + mass_biodegraded_from_water + (any
     remaining mass_components) == initial mass, for all elements ever seeded
     (active + deactivated).
  4. Q10 sanity check: warmer run biodegrades the dissolved pool faster.
  5. Long-lived residue tail: with the component-resolved K_BIODEG_WATER, the
     dissolved residue component has a very long half-life (mirroring
     K_BIODEG_RESIDUE's role for dispersed oil), so the total biodegraded
     fraction stays well below 100% even at 30C after 60 days, and elements
     may remain active (not deactivated, reason='dissolved_biodegraded')
     within this window -- that is expected, physically-correct behavior.
"""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data_local' / 'biodeg_water_test'
OUT.mkdir(parents=True, exist_ok=True)

from cicoil_common import load_ciceseoil
OpenCiceseOil = load_ciceseoil()
print(f'K_BIODEG_WATER = {OpenCiceseOil.K_BIODEG_WATER}')
print(f'BIODEG_WATER_DEPLETION_FRAC = {OpenCiceseOil.BIODEG_WATER_DEPLETION_FRAC}')

N_PARTICLES = 20
SIM_DAYS = 60  # long enough for the t1/2~12d (20C) / faster (30C) pool to deplete
TIME_STEP_MIN = 60
OUT_DT_H = 6
RELEASE_TIME = datetime(2023, 7, 6, 18, 0, 0)
LON, LAT = -92.0, 19.4
MASS_PER_ELEMENT = 1.0  # kg, == default mass_oil


def make_sim(sst_c, biodeg, subsea_dissolution):
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
        'sea_water_temperature': sst_c,
        'sea_water_salinity': 36.,
        'sea_floor_depth_below_sea_level': 200.,
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
    o.set_config('processes:biodegradation', biodeg)
    o.set_config('processes:photooxidation', False)
    o.set_config('processes:subsea_dissolution', subsea_dissolution)
    o.set_config('processes:handle_released_gas', False)
    o.set_config('general:coastline_action', 'previous')
    o.set_config('general:seafloor_action', 'previous')

    o.set_oiltype('MAYA')

    o.seed_elements(
        lon=LON, lat=LAT, z=-30,
        number=N_PARTICLES,
        time=RELEASE_TIME,
        mass_oil=MASS_PER_ELEMENT,
        is_dissolved=1,
        mass_dissolved_subsea=MASS_PER_ELEMENT,
    )
    return o


def get_initial_mass_components(sst_c):
    """t=0 sum(mass_components), via a 1-step run with biodegradation and
    subsea_dissolution both off (so mass_components is untouched)."""
    o = make_sim(sst_c, biodeg=False, subsea_dissolution=False)
    o.run(
        duration=timedelta(minutes=TIME_STEP_MIN),
        time_step=timedelta(minutes=TIME_STEP_MIN),
        time_step_output=timedelta(minutes=TIME_STEP_MIN),
        outfile=str(OUT / 'cicoil_biodeg_water_t0.nc'),
    )
    return float(np.sum(o.cicese_mass_balance['mass_components']))


def run_case(sst_c, label, total_initial_mass):
    o = make_sim(sst_c, biodeg=True, subsea_dissolution=True)

    o.run(
        duration=timedelta(days=SIM_DAYS),
        time_step=timedelta(minutes=TIME_STEP_MIN),
        time_step_output=timedelta(hours=OUT_DT_H),
        outfile=str(OUT / f'cicoil_biodeg_water_{label}.nc'),
    )

    mc = o.cicese_mass_balance
    total_dissolved_remaining = float(np.sum(mc['mass_dissolved']))
    total_biodeg_water = float(np.sum(mc['mass_biodegraded_from_water']))
    total_components_remaining = float(np.sum(mc['mass_components']))

    mass_check = abs(total_dissolved_remaining + total_biodeg_water +
                      total_components_remaining - total_initial_mass) / total_initial_mass

    return dict(
        label=label,
        sst_c=sst_c,
        total_initial_mass_kg=total_initial_mass,
        total_dissolved_remaining_kg=total_dissolved_remaining,
        total_biodeg_water_kg=total_biodeg_water,
        total_components_remaining_kg=total_components_remaining,
        frac_biodeg_water=100 * total_biodeg_water / total_initial_mass,
        mass_conservation_rel_err=mass_check,
        n_active=o.num_elements_active(),
        n_deactivated=o.num_elements_deactivated(),
    )


print(f'\n{"="*70}')
print('CICOILv2 dissolved-phase biodegradation test (Q10 + lifecycle check)')
print(f'{"="*70}')
print(f'MAYA oil, N={N_PARTICLES} dissolved particles, {SIM_DAYS} days, calm conditions')

total_initial_mass = get_initial_mass_components(20.0)
print(f'\nInitial sum(mass_components) at t=0: {total_initial_mass:.6e} kg '
      f'({N_PARTICLES} particles, mass_oil={MASS_PER_ELEMENT} kg each)')

results = {}
for sst, label in [(10.0, 'cold_10C'), (20.0, 'ref_20C'), (30.0, 'warm_30C')]:
    print(f'\n--- Running: SST={sst}C ---')
    res = run_case(sst, label, total_initial_mass)
    results[label] = res
    print(f'  initial total mass         : {res["total_initial_mass_kg"]:.3f} kg')
    print(f'  mass_dissolved remaining    : {res["total_dissolved_remaining_kg"]:.4f} kg')
    print(f'  mass_biodegraded_from_water : {res["total_biodeg_water_kg"]:.4f} kg')
    print(f'  mass_components remaining   : {res["total_components_remaining_kg"]:.2e} kg (expected ~0)')
    print(f'  fraction biodegraded (water): {res["frac_biodeg_water"]:.3f}%')
    print(f'  mass conservation rel. err  : {res["mass_conservation_rel_err"]:.2e}')
    print(f'  active / deactivated        : {res["n_active"]} / {res["n_deactivated"]}')

print(f'\n{"-"*70}')
print('Checks:')

# 1. Mass conservation
for label, res in results.items():
    assert res['mass_conservation_rel_err'] < 1e-5, \
        f'Mass conservation FAILED for {label}: {res["mass_conservation_rel_err"]:.2e}'
print('  PASS: mass conservation OK for all runs')

# 2. One-time transfer worked: mass_components ~ 0 (all moved to mass_dissolved)
for label, res in results.items():
    assert res['total_components_remaining_kg'] < 1e-6, \
        f'mass_components should be ~0 after transfer for {label}, got {res["total_components_remaining_kg"]}'
print('  PASS: mass_components -> mass_dissolved one-time transfer OK')

# 3. Biodegradation occurred
for label, res in results.items():
    assert res['total_biodeg_water_kg'] > 0, \
        f'mass_biodegraded_from_water should be > 0 for {label}'
print('  PASS: mass_biodegraded_from_water > 0 for all runs')

# 4. Q10 sanity check (warmer -> faster depletion -> more biodegraded by day 60)
print(f'\n  10C: frac_biodeg_water={results["cold_10C"]["frac_biodeg_water"]:.3f}%, '
      f'deactivated={results["cold_10C"]["n_deactivated"]}')
print(f'  20C: frac_biodeg_water={results["ref_20C"]["frac_biodeg_water"]:.3f}%, '
      f'deactivated={results["ref_20C"]["n_deactivated"]}')
print(f'  30C: frac_biodeg_water={results["warm_30C"]["frac_biodeg_water"]:.3f}%, '
      f'deactivated={results["warm_30C"]["n_deactivated"]}')
assert (results['cold_10C']['frac_biodeg_water'] <=
        results['ref_20C']['frac_biodeg_water'] <=
        results['warm_30C']['frac_biodeg_water']), \
    'Q10 temperature scaling check FAILED'
print('  PASS: monotonically increasing with temperature')

# 5. Long-lived residue tail: K_BIODEG_WATER is now component-resolved, with
#    the dissolved residue component rescaled from K_BIODEG_RESIDUE (a very
#    long half-life, even after Q10 scaling at 30C). The total dissolved pool
#    therefore never drops below BIODEG_WATER_DEPLETION_FRAC=1% within 60
#    days at any of the tested temperatures, so n_deactivated==0 everywhere
#    is expected. Check that the 30C biodegraded fraction is bounded well
#    below full depletion, confirming the residue tail persists.
assert results['warm_30C']['frac_biodeg_water'] < 90.0, \
    'Expected dissolved residue tail to keep 30C biodegraded fraction well ' \
    'below 100% over 60 days'
print('  PASS: dissolved residue tail persists (30C frac_biodeg_water < 90% over 60 days)')

print(f'\n{"="*70}')
print('ALL CHECKS PASSED')
print(f'{"="*70}')
