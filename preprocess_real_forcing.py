#!/usr/bin/env python3
"""
preprocess_real_forcing.py
===========================
CF-ready readers for the NohochAlfa 2023 real-forcing CIC-OILv2 run:
  - Mercator GLORYS12 wide (currents uo/vo + SST thetao + salinity so +
    sea surface height zos + mixed layer depth mlotst), Jun1-Aug5 2023
  - ERA5 wide wind (CF-preprocessed, x_wind/y_wind), Jun-Aug 2023
  - ERA5 waves (swh, pp1d -> sea_surface_wave_significant_height /
    sea_surface_wave_period_at_variance_spectral_density_maximum), Jul-Aug 2023

All datasets are opened with xarray, given CF standard_name attributes where
missing, and passed directly to reader_netCDF_CF_generic.Reader(ds) -- no
copies of the (large) source files are written to disk.
"""

from pathlib import Path

import xarray as xr

DATA = Path('/Volumes/DATA_SSD/ADominguez/CLAUDE_SPACE/OILSPILL/DATA/NohochAlfa_2023')
WAVES_DIR = Path('/Users/alexdominguez/ADominguez/CLAUDE_SPACE/OILSPILL/photooxidation/data')


def mercator_reader():
    from opendrift.readers import reader_netCDF_CF_generic as rg
    ds = xr.open_dataset(DATA / 'mercator' / 'mercator_glorys12_wide_2023-06_2023-08.nc')
    ds['thetao'].attrs['standard_name'] = 'sea_water_temperature'
    ds['mlotst'].attrs['standard_name'] = 'ocean_mixed_layer_thickness'
    return rg.Reader(ds, name='mercator_glorys12_wide')


def wind_reader():
    from opendrift.readers import reader_netCDF_CF_generic as rg
    files = sorted((DATA / 'era5_cf').glob('cf_era5_wide_2023-0*.nc'))
    return rg.Reader([str(f) for f in files], name='era5_wide_wind')


def wave_reader():
    from opendrift.readers import reader_netCDF_CF_generic as rg
    ds7 = xr.open_dataset(WAVES_DIR / 'era5_waves_2023_07.nc')
    ds8 = xr.open_dataset(WAVES_DIR / 'era5_waves_2023_08.nc')
    ds = xr.concat([ds7, ds8], dim='valid_time')
    ds = ds.rename({'valid_time': 'time'})
    drop = [v for v in ('number', 'expver', 'step', 'surface')
            if v in ds.coords or v in ds.variables]
    ds = ds.drop_vars(drop)
    ds['swh'].attrs['standard_name'] = 'sea_surface_wave_significant_height'
    ds['pp1d'].attrs['standard_name'] = 'sea_surface_wave_period_at_variance_spectral_density_maximum'
    return rg.Reader(ds, name='era5_waves')


def all_readers():
    return [mercator_reader(), wind_reader(), wave_reader()]


if __name__ == '__main__':
    for r in all_readers():
        print(r.name, '->', sorted(r.variables))
        print('   time:', r.start_time, '-', r.end_time)
