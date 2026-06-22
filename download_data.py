"""
CICOILv2.0 — Data Download Utilities
======================================
Functions to download and pre-process all forcing data required for a
CICOILv2.0 simulation:

  - Mercator Ocean GLORYS12 (ocean currents, T, S, SSH, MLD)
  - ERA5 (10-m wind, surface solar radiation, 2-m temperature)
  - ERA5 wave fields (Hs wind-sea, Hs swell, peak period, Stokes drift)
  - CROCO native output (via NEMO-optimized reader, no download needed)

Usage example
-------------
>>> from cicoilv2.download_data import download_glorys12, download_era5_winds, download_era5_waves
>>> download_glorys12(
...     lon_min=-96, lon_max=-86, lat_min=17, lat_max=23,
...     date_start='2023-07-01', date_end='2023-08-05',
...     outfile='mercator_glorys12_nohoch.nc'
... )
>>> download_era5_winds(
...     lon_min=-96, lon_max=-86, lat_min=17, lat_max=23,
...     year=2023, months=[7, 8],
...     outfile='era5_winds_nohoch.nc'
... )
>>> download_era5_waves(
...     lon_min=-96, lon_max=-86, lat_min=17, lat_max=23,
...     year=2023, months=[7, 8],
...     outfile='era5_waves_nohoch.nc'
... )

Notes
-----
Requires:
  - copernicusmarine (pip install copernicusmarine)
  - cdsapi          (pip install cdsapi, with ~/.cdsapirc configured)
  - xarray, netCDF4, numpy
"""

from __future__ import annotations

import os
import warnings
import logging
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Tuple, Union, Optional

import numpy as np

logger = logging.getLogger(__name__)

# ── Mercator GLORYS12 dataset IDs ─────────────────────────────────────────────
# Variables are split across multiple ANFC datasets; GLORYS12 reanalysis is one file.
_GLORYS12_ID = 'cmems_mod_glo_phy_my_0.083deg_P1D-m'

# ANFC near-real-time (split by variable group)
_ANFC_IDS = {
    'currents':    'cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m',
    'temperature': 'cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m',
    'salinity':    'cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m',
    'ssh':         'cmems_mod_glo_phy-ssh_anfc_0.083deg_P1D-m',
    'mld':         'cmems_mod_glo_phy-mld_anfc_0.083deg_P1D-m',
}

# CF standard names for OpenDrift / reader_netCDF_CF_generic
_GLORYS12_VARS = {
    'uo': 'x_sea_water_velocity',
    'vo': 'y_sea_water_velocity',
    'thetao': 'sea_water_temperature',
    'so': 'sea_water_salinity',
    'zos': 'sea_surface_height',
    'mlotst': 'ocean_mixed_layer_thickness',
}


def download_glorys12(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    date_start: str,
    date_end: str,
    outfile: Union[str, Path] = 'glorys12.nc',
    depth_min: float = 0.5,
    depth_max: float = 200.0,
) -> Path:
    """Download GLORYS12 reanalysis ocean forcing for CICOILv2.0.

    Uses the *single* GLORYS12 dataset ID (all variables in one file).
    Suitable for hindcast simulations (up to ~2021 depending on availability).
    For near-real-time simulations use :func:`download_mercator_anfc`.

    Parameters
    ----------
    lon_min, lon_max : float
        Longitude bounds (degrees East).
    lat_min, lat_max : float
        Latitude bounds (degrees North).
    date_start, date_end : str
        ISO-8601 date strings, e.g. ``'2023-07-01'``.
    outfile : str or Path
        Output NetCDF file path.
    depth_min, depth_max : float
        Depth range in metres (positive down). Default 0.5–200 m.

    Returns
    -------
    Path
        Path to the downloaded NetCDF file.

    Notes
    -----
    Requires ``copernicusmarine`` package and valid CMEMS credentials.
    Set credentials via ``copernicusmarine login`` or environment variables
    ``COPERNICUSMARINE_SERVICE_USERNAME`` / ``COPERNICUSMARINE_SERVICE_PASSWORD``.
    """
    try:
        import copernicusmarine as cm
    except ImportError:
        raise ImportError(
            "copernicusmarine is required. Install with: pip install copernicusmarine"
        )

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    variables = list(_GLORYS12_VARS.keys())

    logger.info(
        f"Downloading GLORYS12 {date_start}–{date_end} "
        f"lon=[{lon_min},{lon_max}] lat=[{lat_min},{lat_max}]"
    )

    cm.subset(
        dataset_id=_GLORYS12_ID,
        variables=variables,
        minimum_longitude=lon_min,
        maximum_longitude=lon_max,
        minimum_latitude=lat_min,
        maximum_latitude=lat_max,
        start_datetime=f"{date_start}T00:00:00",
        end_datetime=f"{date_end}T23:59:59",
        minimum_depth=depth_min,
        maximum_depth=depth_max,
        output_filename=str(outfile),
    )

    logger.info(f"GLORYS12 saved to {outfile} ({outfile.stat().st_size/1e6:.1f} MB)")
    return outfile


def download_mercator_anfc(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    date_start: str,
    date_end: str,
    outfile: Union[str, Path] = 'mercator_anfc.nc',
    depth_min: float = 0.5,
    depth_max: float = 200.0,
    tmpdir: Optional[Union[str, Path]] = None,
) -> Path:
    """Download Mercator ANFC near-real-time ocean forcing for CICOILv2.0.

    ANFC variables are split across five separate dataset IDs (currents,
    temperature, salinity, SSH, MLD). This function downloads each to a
    temporary file, merges them with xarray, and saves one combined NetCDF.

    Parameters
    ----------
    lon_min, lon_max, lat_min, lat_max : float
        Spatial bounds.
    date_start, date_end : str
        ISO-8601 date strings.
    outfile : str or Path
        Merged output NetCDF path.
    depth_min, depth_max : float
        Depth range in metres.
    tmpdir : str or Path, optional
        Directory for temporary per-variable files. Defaults to same dir as outfile.

    Returns
    -------
    Path
        Path to merged NetCDF file.
    """
    try:
        import copernicusmarine as cm
        import xarray as xr
    except ImportError:
        raise ImportError(
            "copernicusmarine and xarray are required. "
            "Install with: pip install copernicusmarine xarray"
        )

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    tmpdir = Path(tmpdir) if tmpdir else outfile.parent

    var_groups = {
        'currents':    (['uo', 'vo'],  _ANFC_IDS['currents']),
        'temperature': (['thetao'],    _ANFC_IDS['temperature']),
        'salinity':    (['so'],        _ANFC_IDS['salinity']),
        'ssh':         (['zos'],       _ANFC_IDS['ssh']),
        'mld':         (['mlotst'],    _ANFC_IDS['mld']),
    }

    tmp_files = []
    datasets = []
    try:
        for grp, (variables, ds_id) in var_groups.items():
            tmp = tmpdir / f'_anfc_{grp}.nc'
            logger.info(f"  Downloading ANFC {grp} ({ds_id}) …")
            cm.subset(
                dataset_id=ds_id,
                variables=variables,
                minimum_longitude=lon_min, maximum_longitude=lon_max,
                minimum_latitude=lat_min,  maximum_latitude=lat_max,
                start_datetime=f"{date_start}T00:00:00",
                end_datetime=f"{date_end}T23:59:59",
                minimum_depth=depth_min,   maximum_depth=depth_max,
                output_filename=str(tmp),
            )
            tmp_files.append(tmp)
            datasets.append(xr.open_dataset(tmp))

        # Merge and save
        logger.info("  Merging ANFC variable groups …")
        merged = xr.merge(datasets)
        merged.to_netcdf(outfile)
        logger.info(f"ANFC merged → {outfile} ({outfile.stat().st_size/1e6:.1f} MB)")
    finally:
        for ds in datasets:
            ds.close()
        for tmp in tmp_files:
            if tmp.exists():
                tmp.unlink()

    return outfile


def download_era5_winds(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    year: int,
    months: List[int],
    outfile: Union[str, Path] = 'era5_winds.nc',
    time_step: str = '1h',
) -> Path:
    """Download ERA5 wind and atmospheric surface fields.

    Downloads:
    - ``u10`` / ``v10`` — 10-m U/V wind components (m s⁻¹)
    - ``msl`` — mean sea level pressure (Pa)
    - ``t2m`` — 2-m temperature (K)
    - ``ssr`` — surface solar radiation downwards (J m⁻²)

    The output file is CF-preprocessed: ``valid_time`` dimension renamed to
    ``time``, and CF ``standard_name`` attributes added for OpenDrift compatibility.

    Parameters
    ----------
    lon_min, lon_max, lat_min, lat_max : float
        Spatial extent (degrees).
    year : int
        Year to download.
    months : list of int
        Months to download (1–12).
    outfile : str or Path
        Output path.
    time_step : str
        Temporal resolution: ``'1h'`` (hourly) or ``'6h'`` (6-hourly).

    Returns
    -------
    Path
        Path to the CF-preprocessed NetCDF file.
    """
    try:
        import cdsapi
        import xarray as xr
    except ImportError:
        raise ImportError(
            "cdsapi and xarray required. Install: pip install cdsapi xarray\n"
            "Then configure ~/.cdsapirc — see https://cds.climate.copernicus.eu/api-how-to"
        )

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    raw = outfile.with_suffix('.raw.nc')

    c = cdsapi.Client()
    logger.info(f"ERA5 winds {year} months={months}, lon=[{lon_min},{lon_max}]")

    c.retrieve('reanalysis-era5-single-levels', {
        'product_type':  'reanalysis',
        'variable': [
            '10m_u_component_of_wind',
            '10m_v_component_of_wind',
            'mean_sea_level_pressure',
            '2m_temperature',
            'surface_solar_radiation_downwards',
        ],
        'year':  str(year),
        'month': [f'{m:02d}' for m in months],
        'day':   [f'{d:02d}' for d in range(1, 32)],
        'time':  [f'{h:02d}:00' for h in range(0, 24,
                   1 if time_step == '1h' else 6)],
        'area':  [lat_max, lon_min, lat_min, lon_max],  # N W S E
        'format': 'netcdf',
        'download_format': 'unarchived',
    }, str(raw))

    # CF preprocessing
    ds = xr.open_dataset(raw)
    ds = _preprocess_era5_cf(ds)
    ds.to_netcdf(outfile)
    ds.close()
    raw.unlink()

    logger.info(f"ERA5 winds → {outfile} ({outfile.stat().st_size/1e6:.1f} MB)")
    return outfile


def download_era5_waves(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    year: int,
    months: List[int],
    outfile: Union[str, Path] = 'era5_waves.nc',
    time_step: str = '1h',
) -> Path:
    """Download ERA5 wave fields for CICOIL wave-submergence and Stokes drift.

    Downloads:
    - ``shww`` — significant wave height of wind waves (m)  ← used for f_surface
    - ``swh``  — significant wave height of combined wind+swell (m)
    - ``pp1d`` — peak wave period (s)
    - ``p1ww`` — mean period of wind waves (s)
    - ``ust``  — eastward surface Stokes drift (m s⁻¹)
    - ``vst``  — northward surface Stokes drift (m s⁻¹)

    Parameters
    ----------
    lon_min, lon_max, lat_min, lat_max : float
        Spatial extent.
    year : int
        Year.
    months : list of int
        Months.
    outfile : str or Path
        Output path.
    time_step : str
        ``'1h'`` or ``'3h'``.

    Returns
    -------
    Path
        CF-preprocessed NetCDF.
    """
    try:
        import cdsapi
        import xarray as xr
    except ImportError:
        raise ImportError("cdsapi and xarray required.")

    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    raw = outfile.with_suffix('.raw.nc')

    c = cdsapi.Client()
    logger.info(f"ERA5 waves {year} months={months}")

    c.retrieve('reanalysis-era5-single-levels', {
        'product_type':  'reanalysis',
        'variable': [
            'significant_height_of_wind_waves',
            'significant_height_of_combined_wind_waves_and_swell',
            'peak_wave_period',
            'mean_wave_period_of_wind_waves',
            'eastward_surface_stokes_drift',
            'northward_surface_stokes_drift',
        ],
        'year':  str(year),
        'month': [f'{m:02d}' for m in months],
        'day':   [f'{d:02d}' for d in range(1, 32)],
        'time':  [f'{h:02d}:00' for h in range(0, 24,
                   1 if time_step == '1h' else 3)],
        'area':  [lat_max, lon_min, lat_min, lon_max],
        'format': 'netcdf',
        'download_format': 'unarchived',
    }, str(raw))

    ds = xr.open_dataset(raw)
    ds = _preprocess_era5_cf(ds)
    ds.to_netcdf(outfile)
    ds.close()
    raw.unlink()

    logger.info(f"ERA5 waves → {outfile} ({outfile.stat().st_size/1e6:.1f} MB)")
    return outfile


def _preprocess_era5_cf(ds):
    """Add CF standard_name attributes and normalise time dimension.

    OpenDrift's reader_netCDF_CF_generic maps variables by standard_name.
    CDS API v2 uses ``valid_time`` instead of ``time``; this function renames it.

    Parameters
    ----------
    ds : xarray.Dataset
        Raw ERA5 dataset from CDS.

    Returns
    -------
    xarray.Dataset
        CF-compliant dataset ready for OpenDrift.
    """
    import xarray as xr

    # Rename valid_time → time (CDS API v2 change)
    if 'valid_time' in ds.dims:
        ds = ds.rename({'valid_time': 'time'})
    if 'valid_time' in ds.coords and 'time' not in ds.coords:
        ds = ds.rename({'valid_time': 'time'})

    # CF standard_name mapping
    cf_names = {
        'u10':  'x_wind',
        'v10':  'y_wind',
        'msl':  'air_pressure_at_sea_level',
        't2m':  'air_temperature',
        'ssr':  'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time',
        'shww': 'sea_surface_wind_wave_significant_height',
        'swh':  'sea_surface_wave_significant_height',
        'pp1d': 'sea_surface_wave_period_at_variance_spectral_density_maximum',
        'p1ww': 'sea_surface_wind_wave_mean_period',
        'ust':  'sea_surface_wave_stokes_drift_x_velocity',
        'vst':  'sea_surface_wave_stokes_drift_y_velocity',
    }
    for var, std_name in cf_names.items():
        if var in ds:
            ds[var].attrs['standard_name'] = std_name

    # Rename lat/lon if needed
    rename_map = {}
    if 'latitude' in ds.dims and 'lat' not in ds.dims:
        rename_map['latitude'] = 'lat'
    if 'longitude' in ds.dims and 'lon' not in ds.dims:
        rename_map['longitude'] = 'lon'
    if rename_map:
        ds = ds.rename(rename_map)

    return ds


def preprocess_era5_zip(zip_nc_path: Union[str, Path],
                         outfile: Optional[Union[str, Path]] = None) -> Path:
    """Unpack a ZIP-disguised ERA5 file (old CDS API format).

    Older CDS downloads save ZIP archives with ``.nc`` extension. Each ZIP
    contains two inner NC files (instant + accumulated variables). This function
    extracts and merges them.

    Parameters
    ----------
    zip_nc_path : str or Path
        The ``.nc`` file that is actually a ZIP archive.
    outfile : str or Path, optional
        Output merged NetCDF. Defaults to ``zip_nc_path`` with ``_cf`` suffix.

    Returns
    -------
    Path
        Merged CF NetCDF file.
    """
    import zipfile
    import tempfile
    import xarray as xr

    zip_nc_path = Path(zip_nc_path)
    if outfile is None:
        outfile = zip_nc_path.with_stem(zip_nc_path.stem + '_cf')

    outfile = Path(outfile)

    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(zip_nc_path) as zf:
            zf.extractall(tmpdir)
        inner = sorted(Path(tmpdir).glob('*.nc'))
        datasets = [xr.open_dataset(p) for p in inner]
        merged = xr.merge(datasets)
        merged = _preprocess_era5_cf(merged)
        merged.to_netcdf(outfile)
        for ds in datasets:
            ds.close()

    logger.info(f"ERA5 ZIP unpacked and merged → {outfile}")
    return outfile


def check_domain_coverage(
    forcing_files: dict,
    spill_lon: float,
    spill_lat: float,
    sim_days: int,
    buffer_deg: float = 3.0,
) -> dict:
    """Verify that all forcing files cover the required spatial and temporal domain.

    Parameters
    ----------
    forcing_files : dict
        Mapping of label → file path, e.g.
        ``{'ocean': 'glorys12.nc', 'wind': 'era5_winds.nc'}``.
    spill_lon, spill_lat : float
        Spill location (degrees).
    sim_days : int
        Simulation duration in days.
    buffer_deg : float
        Minimum domain margin around spill location (degrees). Default 3.

    Returns
    -------
    dict
        For each file: ``{'ok': bool, 'lon_ok': bool, 'lat_ok': bool,
        'time_ok': bool, 'n_times': int}``.
    """
    try:
        import xarray as xr
    except ImportError:
        raise ImportError("xarray required.")

    results = {}
    for label, fpath in forcing_files.items():
        fpath = Path(fpath)
        if not fpath.exists():
            results[label] = {'ok': False, 'error': 'File not found'}
            continue
        try:
            ds = xr.open_dataset(fpath)
            lon = ds.get('lon', ds.get('longitude', None))
            lat = ds.get('lat', ds.get('latitude', None))
            t   = ds.get('time', ds.get('valid_time', None))

            lon_ok = lat_ok = t_ok = False
            n_times = 0
            if lon is not None:
                lon_ok = (float(lon.min()) <= spill_lon - buffer_deg and
                          float(lon.max()) >= spill_lon + buffer_deg)
            if lat is not None:
                lat_ok = (float(lat.min()) <= spill_lat - buffer_deg and
                          float(lat.max()) >= spill_lat + buffer_deg)
            if t is not None:
                n_times = len(t)
                t_ok = n_times >= sim_days * 1  # at least 1 step/day

            results[label] = {
                'ok':      lon_ok and lat_ok and t_ok,
                'lon_ok':  lon_ok,
                'lat_ok':  lat_ok,
                'time_ok': t_ok,
                'n_times': n_times,
                'lon_range': (float(lon.min()), float(lon.max())) if lon is not None else None,
                'lat_range': (float(lat.min()), float(lat.max())) if lat is not None else None,
            }
            ds.close()
        except Exception as e:
            results[label] = {'ok': False, 'error': str(e)}

    # Print summary
    print("\nForcing domain check:")
    print(f"  Spill: ({spill_lon:.2f}°E, {spill_lat:.2f}°N), "
          f"{sim_days}d sim, ±{buffer_deg}° buffer required")
    for label, r in results.items():
        status = '✓' if r.get('ok') else '✗'
        print(f"  {status} {label}: {r}")
    return results


# ── CROCO helper ───────────────────────────────────────────────────────────────

def list_croco_files(croco_dir: Union[str, Path],
                     pattern: str = '*.nc') -> list:
    """List CROCO output NetCDF files matching a pattern.

    Parameters
    ----------
    croco_dir : str or Path
        Directory containing CROCO output.
    pattern : str
        Glob pattern. Default ``'*.nc'``.

    Returns
    -------
    list of Path
        Sorted list of matching files.
    """
    croco_dir = Path(croco_dir)
    files = sorted(croco_dir.glob(pattern))
    logger.info(f"Found {len(files)} CROCO files in {croco_dir}")
    return files
