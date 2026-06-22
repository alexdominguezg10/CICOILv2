#!/usr/bin/env python3
"""
test_emulsification_viscosity.py
=================================
Test: viscosity-stability scaling in emulsification_cicese(), ported from
abkatun_2026/oil_weathering.py emulsification_rate() (Fingas & Fieldhouse
2004): emulsification rate is scaled by

    stability = clip(nu_dry/EMUL_STABILITY_NU_REF, EMUL_STABILITY_MIN, EMUL_STABILITY_MAX)

where nu_dry is the dry-oil (no-emulsion), evaporation-corrected kinematic
viscosity [cSt].

Setup: MAYA oil, 10 particles/case, surface release (z=0), evaporation OFF
(isolates the temperature -> viscosity -> stability path; with evaporation
off, fraction_evaporated=0 so nu_dry == kvis_at_temp(T) exactly), constant
10 m/s wind, emulsification ON, 5-day simulation, 60-min timestep.

For MAYA:
  nu_dry(5C)  ~ 709 cSt  -> stability_default = clip(0.709, 0.3, 3.0) = 0.709
  nu_dry(30C) ~  94 cSt  -> stability_default = clip(0.094, 0.3, 3.0) = 0.3

Checks:
  1. With default EMUL_STABILITY_MIN/MAX (0.3/3.0): cold (higher nu_dry,
     stability=0.709) emulsifies FASTER than warm (stability=0.3) ->
     water_fraction(5C) > water_fraction(30C) at end of run.
  2. With EMUL_STABILITY_MIN = EMUL_STABILITY_MAX = 1.0 (stability == 1
     always, viscosity-independent): water_fraction(5C) == water_fraction(30C)
     to high precision -> recovers the original (pre-port) NOAA Bullwinkle
     behaviour, confirming backward compatibility.
  3. For both temperatures, stability_default < 1.0 here, so
     water_fraction(default) < water_fraction(stability=1) -- the viscosity
     scaling suppresses emulsification relative to the unscaled baseline for
     this (fresh, evaporation-off) MAYA oil at both temperatures.
  4. EMUL_RATE_SCALE: with EMUL_STABILITY_MIN = EMUL_STABILITY_MAX =
     EMUL_RATE_SCALE = 1.0 (full backward compat to the stock NOAA Bullwinkle
     rate), water_fraction(cold) is much LARGER than with the default
     EMUL_RATE_SCALE (0.0142, stability=1.0) -- confirming EMUL_RATE_SCALE
     suppresses the water-uptake rate as documented.
"""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'data_local' / 'emulsification_viscosity_test'
OUT.mkdir(parents=True, exist_ok=True)

from cicoil_common import load_ciceseoil
OpenCiceseOil = load_ciceseoil()
print(f'EMUL_STABILITY_NU_REF/MIN/MAX (default) = '
      f'{OpenCiceseOil.EMUL_STABILITY_NU_REF} / '
      f'{OpenCiceseOil.EMUL_STABILITY_MIN} / '
      f'{OpenCiceseOil.EMUL_STABILITY_MAX}')
print(f'EMUL_RATE_SCALE (default) = {OpenCiceseOil.EMUL_RATE_SCALE}')

N_PARTICLES = 10
SIM_HOURS = 4  # > emul_time=7200s=2h, but well short of water_fraction saturation at max_water
TIME_STEP_MIN = 5
OUT_DT_H = 1
RELEASE_TIME = datetime(2023, 7, 6, 18, 0, 0)
LON, LAT = -92.0, 19.4


def run_case(label, sst_c, stability_min, stability_max, rate_scale=None):
    from opendrift.readers.reader_constant import Reader as CR

    o = OpenCiceseOil(loglevel=30, weathering_model='cicese')
    o.EMUL_STABILITY_MIN = stability_min
    o.EMUL_STABILITY_MAX = stability_max
    if rate_scale is not None:
        o.EMUL_RATE_SCALE = rate_scale

    o.add_reader([CR({
        'x_sea_water_velocity': 0.,
        'y_sea_water_velocity': 0.,
        'upward_sea_water_velocity': 0.,
        'x_wind': 10.,
        'y_wind': 0.,
        'sea_surface_wave_significant_height': 1.0,
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
    o.set_config('processes:emulsification', True)
    o.set_config('processes:dispersion', False)
    o.set_config('processes:spreading', False)
    o.set_config('processes:biodegradation', False)
    o.set_config('processes:photooxidation', False)
    o.set_config('general:coastline_action', 'previous')
    o.set_config('general:seafloor_action', 'previous')

    o.set_oiltype('MAYA')
    o.seed_elements(lon=LON, lat=LAT, z=0, number=N_PARTICLES,
                     time=RELEASE_TIME, max_water=0.9)

    o.run(duration=timedelta(hours=SIM_HOURS),
          time_step=timedelta(minutes=TIME_STEP_MIN),
          time_step_output=timedelta(hours=OUT_DT_H),
          outfile=str(OUT / f'cicoil_{label}.nc'))

    return float(np.mean(o.elements.water_fraction)), \
        float(np.mean(o._oil_viscosity_dry) * 1e6)  # cSt


print(f'\n{"="*70}')
print('CICOILv2 emulsification viscosity-stability scaling test')
print(f'{"="*70}')

wf_cold_default, nu_cold = run_case('cold_default', 5.0, 0.3, 3.0)
wf_warm_default, nu_warm = run_case('warm_default', 30.0, 0.3, 3.0)
wf_cold_unit, _ = run_case('cold_unit', 5.0, 1.0, 1.0)
wf_warm_unit, _ = run_case('warm_unit', 30.0, 1.0, 1.0)
wf_cold_stockrate, _ = run_case('cold_stockrate', 5.0, 1.0, 1.0, rate_scale=1.0)

print(f'\nnu_dry(5C)  = {nu_cold:.2f} cSt  -> stability_default = '
      f'{np.clip(nu_cold/1000., 0.3, 3.0):.4f}')
print(f'nu_dry(30C) = {nu_warm:.2f} cSt  -> stability_default = '
      f'{np.clip(nu_warm/1000., 0.3, 3.0):.4f}')

print(f'\nwater_fraction after {SIM_HOURS} hours:')
print(f'  cold (5C),  default stability               : {wf_cold_default:.6f}')
print(f'  warm (30C), default stability               : {wf_warm_default:.6f}')
print(f'  cold (5C),  stability=1.0, default rate scale: {wf_cold_unit:.6f}')
print(f'  warm (30C), stability=1.0, default rate scale: {wf_warm_unit:.6f}')
print(f'  cold (5C),  stability=1.0, rate_scale=1.0    : {wf_cold_stockrate:.6f}')

# 1. Default stability: cold (higher nu, stability=0.709) emulsifies faster
#    than warm (lower nu, stability=0.3 clipped).
assert wf_cold_default > wf_warm_default, \
    f'Viscosity-stability ordering FAILED: cold={wf_cold_default}, warm={wf_warm_default}'
print('\nPASS: water_fraction(cold, default stability) > '
      'water_fraction(warm, default stability)')

# 2. stability=1.0 (MIN==MAX==1.0) -> viscosity-independent -> cold == warm
rel_diff_unit = abs(wf_cold_unit - wf_warm_unit) / max(wf_cold_unit, 1e-12)
assert rel_diff_unit < 1e-6, \
    f'Backward-compat (stability=1) check FAILED: cold={wf_cold_unit}, warm={wf_warm_unit}'
print('PASS: with EMUL_STABILITY_MIN=MAX=1.0, water_fraction(cold) == '
      f'water_fraction(warm) (rel. diff = {rel_diff_unit:.2e}) '
      '-> recovers viscosity-independent NOAA Bullwinkle behaviour')

# 3. Both stability_default values (<1 here) suppress emulsification
#    relative to the stability=1 baseline.
assert wf_cold_default < wf_cold_unit, \
    f'Cold suppression check FAILED: default={wf_cold_default}, unit={wf_cold_unit}'
assert wf_warm_default < wf_warm_unit, \
    f'Warm suppression check FAILED: default={wf_warm_default}, unit={wf_warm_unit}'
print('PASS: default stability (<1 for both temperatures here) suppresses '
      'emulsification relative to the stability=1 baseline')

# 4. EMUL_RATE_SCALE suppresses the water-uptake rate: with rate_scale=1.0
#    (stock NOAA Bullwinkle rate), water_fraction is much larger than with
#    the default EMUL_RATE_SCALE (0.0142), at fixed stability=1.0.
assert wf_cold_stockrate > wf_cold_unit, \
    (f'EMUL_RATE_SCALE suppression check FAILED: '
     f'rate_scale=1.0 -> {wf_cold_stockrate}, '
     f'default rate_scale ({OpenCiceseOil.EMUL_RATE_SCALE}) -> {wf_cold_unit}')
print(f'PASS: EMUL_RATE_SCALE=1.0 (stock NOAA Bullwinkle rate) gives '
      f'water_fraction={wf_cold_stockrate:.6f} >> default '
      f'EMUL_RATE_SCALE={OpenCiceseOil.EMUL_RATE_SCALE} '
      f'water_fraction={wf_cold_unit:.6f}')

print('\nAll emulsification viscosity-stability checks PASSED.')
