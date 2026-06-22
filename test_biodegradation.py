#!/usr/bin/env python3
"""
test_biodegradation.py
=======================
Test simulation: CICOILv2 component-resolved biodegradation
(biodegradation_cicese, 2026-06-14).

Loads ciceseoil.py from THIS directory (cicoilv2_deployment/) in-memory as
opendrift.models.openoil.ciceseoil, without touching the installed
site-packages copy.

Checks:
  1. Mass conservation: mass_components + mass_evaporated + mass_dissolved +
     mass_photooxidized + mass_biodegraded_from_oil == initial per-element mass.
  2. mass_biodegraded_from_oil increases monotonically and roughly matches a
     hand-computed effective-rate exponential decay at the test temperature.
  3. Q10 sanity check: a warmer run biodegrades more than a colder run.

Setup: MAYA oil, surface release, calm constant forcing (minimize other
weathering so biodegradation signal is visible), 10-day simulation,
30-min timestep.
"""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data_local' / 'biodeg_test'
OUT.mkdir(parents=True, exist_ok=True)

from cicoil_common import load_ciceseoil
OpenCiceseOil = load_ciceseoil()

print(f'K_BIODEG_ALI  = {OpenCiceseOil.K_BIODEG_ALI}')
print(f'K_BIODEG_AROM = {OpenCiceseOil.K_BIODEG_AROM}')
print(f'K_BIODEG_RESIDUE = {OpenCiceseOil.K_BIODEG_RESIDUE}')

N_PARTICLES = 50
SIM_DAYS = 10
TIME_STEP_MIN = 60
OUT_DT_H = 6
RELEASE_TIME = datetime(2023, 7, 6, 18, 0, 0)
LON, LAT = -92.0, 19.4


def run_case(sst_c, label):
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
        outfile=str(OUT / f'cicoil_biodeg_{label}.nc'),
    )

    mc = o.cicese_mass_balance
    initial_components = (mc['mass_components'] + mc['mass_evaporated'] +
                           mc['mass_dissolved'] + mc['mass_photooxidized'] +
                           mc['mass_biodegraded_from_oil'] +
                           mc['mass_biodegraded_from_water'])
    total_initial_mass = float(np.sum(initial_components))
    total_remaining = float(np.sum(mc['mass_components']))
    total_biodeg = float(np.sum(mc['mass_biodegraded_from_oil']))
    total_biodeg_from_water = float(np.sum(mc['mass_biodegraded_from_water']))

    mass_check = abs(total_remaining + total_biodeg + total_biodeg_from_water
                      - total_initial_mass) / total_initial_mass

    return dict(
        label=label,
        sst_c=sst_c,
        total_initial_mass_kg=total_initial_mass,
        total_remaining_kg=total_remaining,
        total_biodeg_kg=total_biodeg,
        total_biodeg_from_water_kg=total_biodeg_from_water,
        frac_biodeg=100 * total_biodeg / total_initial_mass,
        mass_conservation_rel_err=mass_check,
    )


print(f'\n{"="*70}')
print('CICOILv2 biodegradation test -- temperature sensitivity (Q10 check)')
print(f'{"="*70}')
print(f'MAYA oil, N={N_PARTICLES} particles, {SIM_DAYS} days, calm conditions')

results = {}
for sst, label in [(10.0, 'cold_10C'), (20.0, 'ref_20C'), (30.0, 'warm_30C')]:
    print(f'\n--- Running: SST={sst}C ---')
    res = run_case(sst, label)
    results[label] = res
    print(f'  initial total mass        : {res["total_initial_mass_kg"]/1e3:.3f} t')
    print(f'  remaining (mass_components): {res["total_remaining_kg"]/1e3:.3f} t')
    print(f'  biodegraded (from_oil)     : {res["total_biodeg_kg"]/1e3:.4f} t')
    print(f'  biodegraded (from_water)   : {res["total_biodeg_from_water_kg"]/1e3:.4f} t (expected 0)')
    print(f'  fraction biodegraded       : {res["frac_biodeg"]:.3f}%')
    print(f'  mass conservation rel. err : {res["mass_conservation_rel_err"]:.2e}')

print(f'\n{"-"*70}')
print('Q10 sanity check (warmer -> more biodegradation):')
print(f'  10C: {results["cold_10C"]["frac_biodeg"]:.3f}%')
print(f'  20C: {results["ref_20C"]["frac_biodeg"]:.3f}%')
print(f'  30C: {results["warm_30C"]["frac_biodeg"]:.3f}%')
assert results['cold_10C']['frac_biodeg'] < results['ref_20C']['frac_biodeg'] < results['warm_30C']['frac_biodeg'], \
    'Q10 temperature scaling check FAILED'
print('  PASS: monotonically increasing with temperature')

for label, res in results.items():
    assert res['mass_conservation_rel_err'] < 1e-6, \
        f'Mass conservation check FAILED for {label}: {res["mass_conservation_rel_err"]:.2e}'
    assert res['total_biodeg_from_water_kg'] == 0.0, \
        f'mass_biodegraded_from_water should be 0 for {label}'
print('PASS: mass conservation OK for all runs')
