#!/usr/bin/env python3
"""
run_abkatun_julia.py — Abkatun 2026 CICOILv2 simulation with Julia backend.

Usage:
  python run_abkatun_julia.py                     # 3000 particles, 22 days
  python run_abkatun_julia.py -n 30000            # 30k particles
  python run_abkatun_julia.py -n 100000 -d 30     # 100k particles, 30 days
  python run_abkatun_julia.py --no-julia           # Python-only (no Julia)
  python run_abkatun_julia.py -n 10000 --dt 15     # 15-min timestep
"""

import argparse
import os
import time
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr

# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='Abkatun 2026 CICOILv2 + Julia/CUDA simulation')
    p.add_argument('-n', '--particles', type=int, default=3000,
                   help='Number of particles (default: 3000)')
    p.add_argument('-d', '--days', type=int, default=22,
                   help='Simulation duration in days (default: 22)')
    p.add_argument('--dt', type=int, default=30,
                   help='Timestep in minutes (default: 30)')
    p.add_argument('--out-dt', type=int, default=6,
                   help='Output interval in hours (default: 6)')
    p.add_argument('--lon', type=float, default=-92.236,
                   help='Release longitude (default: -92.236)')
    p.add_argument('--lat', type=float, default=19.274,
                   help='Release latitude (default: 19.274)')
    p.add_argument('--depth', type=float, default=0.0,
                   help='Release depth in m, 0=surface (default: 0)')
    p.add_argument('--oil', type=str, default='MAYA',
                   help='Oil type in ADIOS database (default: MAYA)')
    p.add_argument('--release-date', type=str, default='2026-02-06',
                   help='Release date YYYY-MM-DD (default: 2026-02-06)')
    p.add_argument('--no-julia', action='store_true',
                   help='Disable Julia backend (pure Python weathering)')
    p.add_argument('--no-biodeg', action='store_true',
                   help='Disable biodegradation')
    p.add_argument('--no-photoox', action='store_true',
                   help='Disable photooxidation')
    p.add_argument('--julia-project', type=str,
                   default='/LUSTRE/adomingu/CICOILv2/CICOILPhysics.jl',
                   help='Path to CICOILPhysics.jl')
    p.add_argument('--m3-per-hour', type=float, default=1.0,
                   help='Oil release rate in m3/hour (default: 1.0)')
    return p.parse_args()

# ── Paths ────────────────────────────────────────────────────────────────────

HERE = Path(__file__).resolve().parent
FORCING = HERE / 'forcing'
OUTPUT  = HERE / 'output'
OUTPUT.mkdir(exist_ok=True)

F_OCEAN  = FORCING / 'cmems_cur_glorys_abkatun_2026_jan20_mar25.nc'
F_THETAO = FORCING / 'cmems_thetao_anfc_abkatun_2026_jan25_apr30.nc'
F_SO     = FORCING / 'cmems_so_anfc_abkatun_2026_jan25_apr30.nc'
F_WIND   = FORCING / 'abkatun_2026_wind_feb06_mar25_cf.nc'
F_WAVES  = FORCING / 'abkatun_2026_waves_cf.nc'
F_UV     = FORCING / 'abkatun_2026_era5_uv.nc'


def make_readers():
    from opendrift.readers import reader_netCDF_CF_generic as rg
    readers = []

    if F_OCEAN.exists():
        readers.append(rg.Reader(str(F_OCEAN), name='glorys12_currents'))
        print(f'  ✓ Currents: {F_OCEAN.name}')

    if F_THETAO.exists():
        readers.append(rg.Reader(str(F_THETAO), name='anfc_thetao'))
        print(f'  ✓ Temperature: {F_THETAO.name}')

    if F_SO.exists():
        readers.append(rg.Reader(str(F_SO), name='anfc_so'))
        print(f'  ✓ Salinity: {F_SO.name}')

    if F_WIND.exists():
        readers.append(rg.Reader(str(F_WIND), name='era5_wind'))
        print(f'  ✓ Wind: {F_WIND.name}')

    if F_WAVES.exists():
        ds = xr.open_dataset(F_WAVES)
        drop = [v for v in ('number', 'expver') if v in ds.coords or v in ds.variables]
        if drop:
            ds = ds.drop_vars(drop)
        if 'valid_time' in ds.dims and 'time' not in ds.dims:
            ds = ds.rename({'valid_time': 'time'})
        readers.append(rg.Reader(ds, name='era5_waves'))
        print(f'  ✓ Waves: {F_WAVES.name}')

    if F_UV.exists():
        ds = xr.open_dataset(F_UV)
        drop = [v for v in ('number', 'expver') if v in ds.coords or v in ds.variables]
        if drop:
            ds = ds.drop_vars(drop)
        for v in ds.data_vars:
            if 'uv' in v.lower() or 'uv' in ds[v].attrs.get('long_name', '').lower():
                ds[v].attrs['standard_name'] = 'surface_downward_uv_radiation'
        readers.append(rg.Reader(ds, name='era5_uv'))
        print(f'  ✓ UV: {F_UV.name}')

    missing = [f for f in [F_OCEAN, F_THETAO, F_SO, F_WIND, F_WAVES, F_UV]
               if not f.exists()]
    if missing:
        print(f'\n  ⚠ Missing:')
        for f in missing:
            print(f'    {f}')

    return readers


def main():
    args = parse_args()

    release_time = datetime.strptime(args.release_date, '%Y-%m-%d')
    use_julia = not args.no_julia
    label = f'abkatun_{args.days}d_{args.particles}p'
    if use_julia:
        label += '_julia'
    else:
        label += '_python'

    print('='*60)
    print('Abkatun 2026 — CICOILv2 Simulation')
    print('='*60)
    print(f'  Particles:  {args.particles:,}')
    print(f'  Duration:   {args.days} days')
    print(f'  Timestep:   {args.dt} min')
    print(f'  Output:     every {args.out_dt} h')
    print(f'  Release:    {args.lon}°E, {args.lat}°N, z={args.depth}m')
    print(f'  Oil:        {args.oil}, {args.m3_per_hour} m³/h')
    print(f'  Date:       {release_time}')
    print(f'  Backend:    {"Julia" if use_julia else "Python"}')
    if use_julia:
        print(f'  Julia proj: {args.julia_project}')
    print()

    from opendrift.models.openoil.ciceseoil import OpenCiceseOil
    from opendrift.config import CONFIG_LEVEL_BASIC

    o = OpenCiceseOil(loglevel=20, weathering_model='cicese')

    o.required_variables = dict(OpenCiceseOil.required_variables)
    o.required_variables['sea_surface_height'] = {'fallback': 0}
    o._add_config({
        'environment:fallback:sea_surface_height': {
            'type': 'float', 'default': 0., 'min': -10., 'max': 10.,
            'units': 'm', 'description': 'fallback sea surface height',
            'level': CONFIG_LEVEL_BASIC,
        }
    })

    print('Loading forcing data...')
    readers = make_readers()
    if readers:
        o.add_reader(readers)
    print()

    o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 30.0)
    o.set_config('environment:fallback:sea_water_temperature', 25.0)
    o.set_config('environment:fallback:sea_water_salinity', 36.0)
    o.set_config('environment:fallback:ocean_mixed_layer_thickness', 20.0)

    o.set_config('drift:vertical_mixing', True)
    o.set_config('drift:stokes_drift', True)
    o.set_config('general:coastline_action', 'previous')
    o.set_config('general:seafloor_action', 'previous')

    o.set_config('processes:evaporation', True)
    o.set_config('processes:emulsification', True)
    o.set_config('processes:dispersion', True)
    o.set_config('processes:biodegradation', not args.no_biodeg)
    o.set_config('processes:photooxidation', not args.no_photoox)

    if use_julia:
        o.set_config('processes:julia_weathering', True)
        o.set_config('julia:project_path', args.julia_project)
        print(f'Julia weathering: ENABLED')
    else:
        print(f'Julia weathering: DISABLED (Python only)')
    print()

    o.set_config('seed:m3_per_hour', args.m3_per_hour)
    o.set_oiltype(args.oil)
    o.seed_elements(
        lon=args.lon, lat=args.lat, z=args.depth,
        number=args.particles,
        time=release_time,
        max_water=0.9,
    )

    outfile = OUTPUT / f'{label}.nc'
    print(f'Output: {outfile}')
    print(f'Running {args.days}-day simulation with {args.particles:,} particles...')
    print()

    t0 = time.time()

    o.run(
        duration=timedelta(days=args.days),
        time_step=timedelta(minutes=args.dt),
        time_step_output=timedelta(hours=args.out_dt),
        outfile=str(outfile),
    )

    wall_time = time.time() - t0

    # ── Mass balance report ──────────────────────────────────────────────────
    mc = o.cicese_mass_balance
    total_init = float(np.sum(
        mc['mass_components'] + mc['mass_evaporated'] +
        mc['mass_dissolved'] + mc['mass_photooxidized'] +
        mc['mass_biodegraded_from_oil'] + mc['mass_biodegraded_from_water'] +
        mc['mass_op_dissolved'] + mc['mass_op_degraded']
    ))

    print()
    print('='*60)
    print(f'Mass Balance Report  ({label})')
    print('='*60)
    for name, key in [
        ('Remaining (oil)',    'mass_components'),
        ('Evaporated',         'mass_evaporated'),
        ('Dissolved (subsea)', 'mass_dissolved'),
        ('Photooxidized',      'mass_photooxidized'),
        ('Biodegraded (oil)',  'mass_biodegraded_from_oil'),
        ('Biodegraded (water)','mass_biodegraded_from_water'),
        ('OP dissolved',       'mass_op_dissolved'),
        ('OP degraded',        'mass_op_degraded'),
    ]:
        val = float(np.sum(mc[key]))
        pct = 100 * val / total_init if total_init > 0 else 0
        print(f'  {name:25s}: {pct:8.4f}%')

    total_accounted = sum(float(np.sum(mc[k])) for k in [
        'mass_components', 'mass_evaporated', 'mass_dissolved',
        'mass_photooxidized', 'mass_biodegraded_from_oil',
        'mass_biodegraded_from_water', 'mass_op_dissolved', 'mass_op_degraded',
    ])
    rel_err = abs(total_accounted - total_init) / total_init if total_init > 0 else 0
    print(f'  {"Mass conservation err":25s}: {rel_err:.2e}')
    print(f'  {"Wall time":25s}: {wall_time:.1f} s ({wall_time/60:.1f} min)')

    np.savez(
        OUTPUT / f'{label}_summary.npz',
        total_init=total_init,
        mass_remaining=float(np.sum(mc['mass_components'])),
        mass_evaporated=float(np.sum(mc['mass_evaporated'])),
        mass_biodegraded=float(np.sum(mc['mass_biodegraded_from_oil'])),
        mass_photooxidized=float(np.sum(mc['mass_photooxidized'])),
        mass_check=rel_err,
        wall_time=wall_time,
        n_particles=args.particles,
        sim_days=args.days,
        backend='julia' if use_julia else 'python',
    )

    print(f'\nDone. Output: {outfile}')


if __name__ == '__main__':
    main()
