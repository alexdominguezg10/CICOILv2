#!/usr/bin/env python3
"""
preprocess_era5_uv.py
======================
CF-preprocess ERA5 Bay-of-Campeche UV files (era5_uv_campeche_*.nc, variable
`uvb`, hourly-accumulated J m**-2) so `reader_netCDF_CF_generic` can supply
them to CIC-OILv2 as `surface_downward_uv_radiation` (W m**-2).

- valid_time -> time
- uvb (J m**-2, hourly accum) -> surface_downward_uv_radiation (W m**-2, /3600)
- standard_name = 'surface_downward_uv_radiation'
- drops the scalar 'number' coord and 'expver' aux coord (not needed by
  reader_netCDF_CF_generic and can confuse variable detection)

latitude/longitude are already CF-named in the source files, no rename needed.
"""

from pathlib import Path

import xarray as xr

RAW_DIR = Path('/Volumes/DATA_SSD/ADominguez/CLAUDE_SPACE/OILSPILL/photooxidation/'
               'data/era5_uv_campeche')
OUT_DIR = Path(__file__).resolve().parent / 'data_local' / 'era5_uv_cf'


def preprocess_era5_uv(raw_files=None, out_dir=OUT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)

    if raw_files is None:
        raw_files = sorted(RAW_DIR.glob('era5_uv_campeche_*.nc'))

    out_files = []
    for raw in raw_files:
        raw = Path(raw)
        out_f = out_dir / f'cf_{raw.name}'
        if out_f.exists():
            out_files.append(str(out_f))
            continue

        print(f'  ERA5 UV preprocessing: {raw.name} -> {out_f.name}')
        ds = xr.open_dataset(raw)

        if 'valid_time' in ds.dims and 'time' not in ds.dims:
            ds = ds.rename({'valid_time': 'time'})

        drop_vars = [v for v in ('number', 'expver') if v in ds.coords or v in ds.variables]
        if drop_vars:
            ds = ds.drop_vars(drop_vars)

        # J m**-2 (hourly accumulated) -> W m**-2
        uv = ds['uvb'] / 3600.0
        uv.attrs = {
            'units': 'W m-2',
            'long_name': 'Surface downward UV radiation',
            'standard_name': 'surface_downward_uv_radiation',
        }
        ds = ds.drop_vars('uvb')
        ds['surface_downward_uv_radiation'] = uv

        ds.to_netcdf(out_f)
        ds.close()
        out_files.append(str(out_f))

    print(f'  ERA5 UV CF files ready: {len(out_files)}')
    return out_files


if __name__ == '__main__':
    files = preprocess_era5_uv()
    for f in files:
        print(f)
