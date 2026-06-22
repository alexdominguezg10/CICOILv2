#!/usr/bin/env python3
"""
test_photoox_era5uv.py
=======================
Test simulation: CIC-OILv2 photooxidation driven by REAL ERA5 UV
(`surface_downward_uv_radiation`, from era5_uv_campeche_*.nc / `uvb`,
preprocessed by preprocess_era5_uv.py), run WITH drift:vertical_mixing,
to contrast against:
  - this session's synthetic-UV result (3.4989% of aromatic mass, day 5)
  - surface_fate_coupled.py's real-UV 0-D result (3.47% of total oil mass,
    day 5)

Loads ciceseoil.py from THIS directory (cicoilv2_deployment/, with the
2026-06-13 K_PHOTO_AROM rescale + f_surface removal + ERA5 UV support) in
memory as opendrift.models.openoil.ciceseoil, without touching the
installed site-packages copy.

Setup: MAYA oil, surface release, constant wind/wave forcing matching the
Jul6-Aug5 2023 Campeche mean conditions (U10=5.76 m/s, Hs=0.538 m,
Tp=4.36 s), 5-day simulation, 30-min timestep, real ERA5 UV reader.
"""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = Path('/Users/alexdominguez/ADominguez/CLAUDE_SPACE/OILSPILL/nohoch_2023/'
           'data_local/cicoil_photoox_test')
OUT.mkdir(parents=True, exist_ok=True)

from preprocess_era5_uv import preprocess_era5_uv  # noqa: E402
from cicoil_common import load_ciceseoil
OpenCiceseOil = load_ciceseoil()
print(f'K_PHOTO_AROM = {mod.OpenCiceseOil.K_PHOTO_AROM}')
print('required_variables has surface_downward_uv_radiation:',
      'surface_downward_uv_radiation' in OpenCiceseOil.required_variables)

# ── Real ERA5 UV reader (preprocessed cf_era5_uv_campeche_2023.nc) ──────────
uv_files = preprocess_era5_uv(raw_files=[
    '/Volumes/DATA_SSD/ADominguez/CLAUDE_SPACE/OILSPILL/photooxidation/'
    'data/era5_uv_campeche/era5_uv_campeche_2023.nc'])
print(f'ERA5 UV CF file: {uv_files[0]}')

# ── Forcing (Jul6-Aug5 2023 Campeche means) ─────────────────────────────────
U10 = 5.76      # m/s
HS  = 0.538     # m
TP  = 4.36      # s
SST_C = 29.0    # deg C

N_PARTICLES = 300
SIM_DAYS = 5
TIME_STEP_MIN = 30
OUT_DT_H = 6
RELEASE_TIME = datetime(2023, 7, 6, 18, 0, 0)
LON, LAT = -92.0, 19.4


def run_case(vertical_mixing, use_era5_uv, label):
    from opendrift.readers.reader_constant import Reader as CR
    from opendrift.readers.reader_netCDF_CF_generic import Reader as CFReader

    o = OpenCiceseOil(loglevel=30, weathering_model='cicese')

    # 'sea_surface_height' is accessed by OceanDrift.vertical_mixing() but
    # missing from OpenCiceseOil.required_variables (which fully overrides
    # the parent dict) — add it (with config fallback) so vertical_mixing works.
    from opendrift.config import CONFIG_LEVEL_BASIC
    o.required_variables = dict(OpenCiceseOil.required_variables)
    o.required_variables['sea_surface_height'] = {'fallback': 0}
    o._add_config({
        'environment:fallback:sea_surface_height': {
            'type': 'float', 'default': 0., 'min': -10., 'max': 10.,
            'units': 'm', 'description': 'fallback sea surface height',
            'level': CONFIG_LEVEL_BASIC,
        }
    })

    if use_era5_uv:
        uv_reader = CFReader(uv_files[0])
        o.add_reader([uv_reader])

    o.add_reader([CR({
        'x_sea_water_velocity': 0.,
        'y_sea_water_velocity': 0.,
        'upward_sea_water_velocity': 0.,
        'x_wind': U10,
        'y_wind': 0.,
        'sea_surface_wave_significant_height': HS,
        'sea_surface_wave_period_at_variance_spectral_density_maximum': TP,
        'sea_water_temperature': SST_C,
        'sea_water_salinity': 36.,
        'sea_floor_depth_below_sea_level': 50.,
        'ocean_mixed_layer_thickness': 30.,
        'land_binary_mask': 0,
    })])

    o.set_config('drift:vertical_mixing', vertical_mixing)
    o.set_config('drift:stokes_drift', True)
    o.set_config('processes:evaporation', True)
    o.set_config('processes:emulsification', True)
    o.set_config('processes:dispersion', True)
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
        outfile=str(OUT / f'cicoil_photoox_{label}.nc'),
    )

    # Reconstruct initial per-component mass from conservation: initial =
    # remaining + evaporated + dissolved + photooxidized + biodegraded.
    mc = o.cicese_mass_balance
    initial_components = (mc['mass_components'] + mc['mass_evaporated'] +
                           mc['mass_dissolved'] + mc['mass_photooxidized'] +
                           mc['mass_biodegraded_from_oil'])
    total_initial_mass = float(np.sum(initial_components))
    initial_aromatic_mass = float(np.sum(initial_components[:, 8:16]))

    photoox_components = mc['mass_photooxidized'][:, 8:16]
    total_photoox = float(np.sum(photoox_components))

    # surface-time fraction from the output trajectory file (z==0)
    try:
        import xarray as xr
        with xr.open_dataset(OUT / f'cicoil_photoox_{label}.nc') as ds:
            f_surface_mean = float(np.nanmean(ds['z'].values == 0))
    except Exception:
        f_surface_mean = np.nan

    return dict(
        label=label,
        vertical_mixing=vertical_mixing,
        use_era5_uv=use_era5_uv,
        total_initial_mass_kg=total_initial_mass,
        initial_aromatic_mass_kg=initial_aromatic_mass,
        total_photoox_kg=total_photoox,
        frac_of_aromatic=100 * total_photoox / initial_aromatic_mass,
        frac_of_total=100 * total_photoox / total_initial_mass,
        f_surface_mean=f_surface_mean,
    )


print(f'\n{"="*70}')
print('CIC-OILv2 photooxidation test — real ERA5 UV vs synthetic UV')
print(f'{"="*70}')
print(f'MAYA oil, N={N_PARTICLES} particles, {SIM_DAYS} days, '
      f'U10={U10} m/s, Hs={HS} m, Tp={TP} s, SST={SST_C}C')
print(f'Release: {RELEASE_TIME} at ({LON}, {LAT})')

results = {}
for vm, use_uv, label in [(True, True, 'era5uv_vertmix_on'),
                            (True, False, 'synthuv_vertmix_on')]:
    print(f'\n--- Running: vertical_mixing={vm}, era5_uv={use_uv} ---')
    res = run_case(vm, use_uv, label)
    results[label] = res
    print(f'  initial total mass     : {res["total_initial_mass_kg"]/1e3:.2f} t')
    print(f'  initial aromatic mass  : {res["initial_aromatic_mass_kg"]/1e3:.3f} t')
    print(f'  photooxidized (G1-G8)  : {res["total_photoox_kg"]/1e3:.4f} t')
    print(f'  fraction of aromatic   : {res["frac_of_aromatic"]:.3f}%')
    print(f'  fraction of total oil  : {res["frac_of_total"]:.4f}%')
    print(f'  mean f_surface (z==0)  : {res["f_surface_mean"]:.3f}')

print(f'\n{"-"*70}')
era5, synth = results['era5uv_vertmix_on'], results['synthuv_vertmix_on']
print(f'Photooxidized (% of aromatic mass), day {SIM_DAYS}, vertical_mixing=True:')
print(f'  Real ERA5 UV       : {era5["frac_of_aromatic"]:.3f}%  '
      f'(f_surface={era5["f_surface_mean"]:.3f})')
print(f'  Synthetic UV       : {synth["frac_of_aromatic"]:.3f}%  '
      f'(f_surface={synth["f_surface_mean"]:.3f})')

print(f'\nFor reference, surface_fate_coupled.py (0-D, real Jul-Aug 2023 ERA5 UV,')
print(f'naph+phen+pyr only) at day {SIM_DAYS}: 3.47% of total oil mass photooxidized')
print(f'(cumulative, see surface_fate_base.csv).')
print(f'Previous session result (synthetic UV, vertical_mixing=True): 3.4989%')
