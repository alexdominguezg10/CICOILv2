"""Ambient profile for the TAMOC near-field model, built from OpenDrift readers.

Port of CICOIL v1 ``vertical_profile.create_profile`` (written for the OpenDrift 1.9 reader API, ``reader.get_variables``)
to OpenDrift 1.14.9 (``get_variables_interpolated``).  The physics of the step is unchanged: temperature, salinity and the two
horizontal velocity components are read at one position and time at every depth level of the reader, down to the seafloor
(first level with a missing value), and written as a TAMOC ambient-profile netCDF file that ``run_tamoc.run_tamoc`` reads.

Install: copy next to ``ciceseoil.py`` (``opendrift/models/openoil/``).  The first part (``sample_profile``) needs only OpenDrift
and NumPy; ``create_profile`` also needs TAMOC (imported inside the function, so importing this module never requires it).

Differences from the v1 function (all deliberate):
  * several readers are accepted (T and S often come from a different file than the currents); each variable is read from the first
    reader that offers it, the depth levels are the union of the valid levels, and a variable is interpolated in depth between its own;
  * the raw block of ``reader.get_variables`` is read (not ``get_variables_interpolated``, which fills missing values from neighbouring
    cells and levels in OpenDrift 1.14.9 and would hide the seafloor); the cell nearest to the position and the time nearest to
    ``time`` are used, as in v1;
  * the profile stops at the first level where ANY of the four variables is missing (v1 tested only salinity);
  * fewer than two valid levels raises ValueError (v1 returned an empty profile);
  * a level at z = 0 is added (copy of the shallowest value) when the shallowest reader level is below the surface (``add_surface``);
  * ``create_profile`` closes the netCDF file and returns its path (v1 returned the open dataset).
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)

VARIABLES = ('sea_water_temperature', 'sea_water_salinity', 'x_sea_water_velocity', 'y_sea_water_velocity')


def _as_list(readers):
    return list(readers) if isinstance(readers, (list, tuple)) else [readers]


def _reader_for(readers, variable):
    for r in readers:
        if variable in r.variables:
            return r
    raise ValueError("no reader provides '%s' (readers: %s)" % (variable, [getattr(r, 'name', repr(r)) for r in readers]))


def _raw_column(reader, variable, time, lon, lat):
    """Raw (not interpolated) column of one variable at the grid cell nearest to (lon, lat) and the reader time nearest to ``time``.

    Uses ``reader.get_variables`` (the raw block, NaN kept where the model has no data), as v1 did.  OpenDrift's
    ``get_variables_interpolated`` must NOT be used here: it fills missing values from neighbouring cells and levels, so the seafloor
    and land can no longer be recognised.  Returns (depth m positive down, shallowest first; values), both restricted to the valid
    levels above the first missing one."""
    z = np.asarray(getattr(reader, 'z', []), dtype=float).ravel()
    if z.size == 0:
        raise ValueError("reader '%s' has no depth levels (reader.z)" % getattr(reader, 'name', reader))
    z = -np.sort(np.abs(z))                                   # OpenDrift z is negative below the surface; shallowest first
    x, y = reader.lonlat2xy(np.array([float(lon)]), np.array([float(lat)]))
    x0, y0 = float(np.squeeze(x)), float(np.squeeze(y))
    block = reader.get_variables([variable], time, x=np.array([x0]), y=np.array([y0]), z=z)
    data = np.asarray(np.ma.filled(block[variable], np.nan), dtype=float)
    if data.ndim != 3:
        raise ValueError("variable '%s' of reader '%s' is not 3-D (z, y, x): shape %s" % (variable, getattr(reader, 'name', reader), data.shape))
    ix = int(np.argmin(np.abs(np.asarray(block['x'], dtype=float) - x0)))
    iy = int(np.argmin(np.abs(np.asarray(block['y'], dtype=float) - y0)))
    col = data[:, iy, ix]
    zb = -np.abs(np.asarray(block['z'], dtype=float).ravel())
    order = np.argsort(-zb)                                   # shallowest first
    zb, col = zb[order], col[order]
    n_valid = int(np.argmax(~np.isfinite(col))) if not np.all(np.isfinite(col)) else len(col)
    return -zb[:n_valid], col[:n_valid]


def sample_profile(readers, lon, lat, time, z_max=None, add_surface=True):
    """Read T, S, u, v at (lon, lat, time) at every depth level of the readers, from the surface down to the seafloor.

    readers     one OpenDrift reader or a list of readers (3-D z-level readers, e.g. GLORYS12, HYCOM; sigma-coordinate readers such as
                reader_croco_native are not supported here)
    z_max       optional, deepest level to keep (m, positive)
    add_surface add a level at z = 0 when the shallowest level is below the surface
    Each variable is read at the grid cell nearest to (lon, lat) and the reader time nearest to ``time`` (no horizontal or time
    interpolation, as in v1).  The column of a variable ends at its first missing value (seafloor); the profile ends at the shallowest
    of these ends.  Where readers have different levels the output uses the union of the valid levels and each variable is
    interpolated linearly in depth between its own levels.
    Returns a dict of 1-D arrays: ``depth`` (m, positive down, increasing), ``temperature`` (K), ``salinity`` (psu),
    ``ua`` and ``va`` (m s-1).  Temperature is taken as degC (OpenDrift convention) and converted to K; a mean above 150 is
    assumed to be in K already and is left unchanged (with a warning).
    """
    readers = _as_list(readers)
    cols = {}
    for v in VARIABLES:
        d, c = _raw_column(_reader_for(readers, v), v, time, lon, lat)
        if len(d) < 1:
            raise ValueError('no valid %s at lon=%s lat=%s: the position is on land or outside the data; use the nearest ocean point' % (v, lon, lat))
        cols[v] = (d, c)
    deepest = min(d[-1] for d, _ in cols.values())            # shallowest end of the four columns
    depth = np.unique(np.round(np.concatenate([d[d <= deepest + 1e-9] for d, _ in cols.values()]), 6))
    if z_max is not None:
        depth = depth[depth <= abs(z_max) + 1e-9]
    if len(depth) < 2:
        raise ValueError('fewer than 2 valid depth levels at lon=%s lat=%s (%d found): the position is on land or in water too '
                         'shallow for this grid; use the nearest ocean point' % (lon, lat, len(depth)))
    out = {v: np.interp(depth, cols[v][0], cols[v][1]) for v in VARIABLES}
    if add_surface and depth[0] > 0.0:
        depth = np.concatenate([[0.0], depth])
        out = {v: np.concatenate([[a[0]], a]) for v, a in out.items()}
    temperature = out['sea_water_temperature'].copy()
    if np.nanmean(temperature) > 150.0:
        logger.warning('sea_water_temperature looks like Kelvin already (mean %.1f); not adding 273.15', np.nanmean(temperature))
    else:
        temperature = temperature + 273.15
    return dict(depth=depth, temperature=temperature, salinity=out['sea_water_salinity'],
                ua=out['x_sea_water_velocity'], va=out['y_sea_water_velocity'])


def create_profile(nc_file, readers, lon, lat, time, source='OpenDrift readers', summary='Ambient profile from OpenDrift readers',
                   sea_name='Gulf of Mexico', z_max=None, add_surface=True):
    """Write a TAMOC ambient-profile netCDF file (T, S, pressure, ua, va against depth) and return its path.

    Same file layout as the v1 function: variables z, temperature, salinity, pressure, ua, va (m, K, psu, Pa, m/s, m/s), pressure from
    ``tamoc.ambient.compute_pressure``.  ``time`` is a datetime (the profile's reference time)."""
    from netCDF4 import date2num
    from tamoc import ambient

    p = sample_profile(readers, lon, lat, time, z_max=z_max, add_surface=add_surface)
    pressure = ambient.compute_pressure(p['depth'], p['temperature'], p['salinity'], 0)
    p_time = date2num(time, units='seconds since 1970-01-01 00:00:00 0:00', calendar='julian')
    nc = ambient.create_nc_db(nc_file, summary, source, sea_name, lat, lon, p_time)
    data = np.vstack((p['depth'], p['temperature'], p['salinity'], pressure, p['ua'], p['va'])).transpose()
    symbols = ['z', 'temperature', 'salinity', 'pressure', 'ua', 'va']
    units = ['m', 'K', 'psu', 'Pa', 'm/s', 'm/s']
    comments = ['OpenDrift reader data'] * len(symbols)
    nc = ambient.fill_nc_db(nc, data, symbols, units, comments, z_col=0)
    nc.close()
    logger.info('TAMOC ambient profile written: %s (%d levels, 0..%.1f m)', nc_file, len(p['depth']), p['depth'][-1])
    return nc_file
