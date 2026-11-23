#!/usr/bin/env python3
"""Tests for vertical_profile.py (OpenDrift 1.14.9) on synthetic readers; no real model data.

Synthetic ocean: currents in one file, temperature/salinity in another (different depth levels), linear fields, NaN at the seafloor
below 30 m, and one all-NaN ("land") column.  Checks: values per level, union of levels across readers, surface level, K conversion,
stop at the seafloor, z_max, errors (land point, missing variable), and (only if TAMOC is importable) the TAMOC profile file.
Run:  python test_vertical_profile.py   (CICOIL_FILE not needed; vertical_profile.py is imported from this folder)
"""
import os, sys, tempfile, logging, warnings
from datetime import datetime
import numpy as np, xarray as xr

logging.disable(logging.WARNING); warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vertical_profile as vp
from opendrift.readers import reader_netCDF_CF_generic as rg

FAILS = []
def check(name, cond, detail=''):
    print(('  OK    ' if cond else '  FAIL  ') + name + (f'  [{detail}]' if detail else ''))
    if not cond: FAILS.append(name)

tmp = tempfile.mkdtemp()
t = np.array(['2026-02-06T00:00', '2026-02-06T12:00'], dtype='datetime64[ns]')
lon = np.linspace(-93.0, -91.5, 7); lat = np.linspace(18.5, 20.0, 7)
dep_uv = np.array([0.5, 5., 10., 20., 30.])           # currents: 5 levels
dep_ts = np.array([0.5, 5., 10., 20., 30., 50.])      # T, S: one more, deeper level
LAND = 5                                              # lon index of the all-NaN column
def field(dep, fn):
    a = np.broadcast_to(fn(dep)[None, :, None, None], (2, len(dep), len(lat), len(lon))).astype('f4').copy()
    a[:, dep >= 30.0, :, :] = np.nan                  # seafloor at ~25 m: levels 30 m and deeper are missing
    a[:, :, :, LAND] = np.nan                         # land column
    return a
def write(name, dep, variables):
    ds = xr.Dataset({k: (('time', 'depth', 'latitude', 'longitude'), v, {'standard_name': s}) for k, (v, s) in variables.items()},
                    coords={'time': t, 'depth': ('depth', dep, {'standard_name': 'depth', 'positive': 'down', 'units': 'm'}),
                            'latitude': ('latitude', lat, {'standard_name': 'latitude', 'units': 'degrees_north'}),
                            'longitude': ('longitude', lon, {'standard_name': 'longitude', 'units': 'degrees_east'})})
    p = os.path.join(tmp, name); ds.to_netcdf(p); return p
p_uv = write('uv.nc', dep_uv, {'uo': (field(dep_uv, lambda z: 0.10 + 0.01 * z), 'eastward_sea_water_velocity'),
                               'vo': (field(dep_uv, lambda z: -0.05 + 0 * z), 'northward_sea_water_velocity')})
p_ts = write('ts.nc', dep_ts, {'thetao': (field(dep_ts, lambda z: 28.0 - 0.2 * z), 'sea_water_potential_temperature'),
                               'so': (field(dep_ts, lambda z: 36.0 + 0.01 * z), 'sea_water_salinity')})
r_uv = rg.Reader(p_uv, name='uv'); r_ts = rg.Reader(p_ts, name='ts')
time = datetime(2026, 2, 6, 6); LON, LAT = float(lon[2]), float(lat[3])

print('1. T/S and currents from two different readers')
P = vp.sample_profile([r_uv, r_ts], LON, LAT, time)
d = P['depth']
check('depth starts at 0 (surface level added) and increases', d[0] == 0.0 and np.all(np.diff(d) > 0), str(np.round(d, 2)))
check('stops above the seafloor (all depths < 30 m)', d[-1] < 30.0, f'deepest {d[-1]} m')
check('at least levels 0, 0.5, 5, 10 present', len(d) >= 4, f'{len(d)} levels')
i5 = int(np.argmin(np.abs(d - 5.0))); check('level 5 m is exactly 5 m', abs(d[i5] - 5.0) < 1e-6)
check('temperature at 5 m = 28 - 0.2*5 degC, in K', abs(P['temperature'][i5] - (28.0 - 1.0 + 273.15)) < 1e-3, f"{P['temperature'][i5]:.3f}")
check('salinity at 5 m = 36.05', abs(P['salinity'][i5] - 36.05) < 1e-4, f"{P['salinity'][i5]:.4f}")
check('ua at 5 m = 0.15 and va = -0.05', abs(P['ua'][i5] - 0.15) < 1e-4 and abs(P['va'][i5] + 0.05) < 1e-4, f"ua {P['ua'][i5]:.4f} va {P['va'][i5]:.4f}")
check('surface level copies the shallowest value', abs(P['temperature'][0] - P['temperature'][1]) < 1e-9 and d[1] > 0)

print('1b. currents on different depth levels than T/S (linear fields: interpolation must be exact)')
dep_uv2 = np.array([1., 8., 15., 25., 30.])
p_uv2 = write('uv2.nc', dep_uv2, {'uo': (field(dep_uv2, lambda z: 0.10 + 0.01 * z), 'eastward_sea_water_velocity'),
                                  'vo': (field(dep_uv2, lambda z: -0.05 + 0 * z), 'northward_sea_water_velocity')})
Q = vp.sample_profile([rg.Reader(p_uv2, name='uv2'), r_ts], LON, LAT, time)
dq = Q['depth']
check('union of levels from both readers, up to the shallowest seafloor (20 m)', abs(dq[-1] - 20.0) < 1e-6 and 8.0 in np.round(dq, 6) and 5.0 in np.round(dq, 6),
      str(np.round(dq, 2)))
check('ua interpolated at the T/S-only level 5 m = 0.10 + 0.01*5', abs(Q['ua'][int(np.argmin(np.abs(dq - 5)))] - 0.15) < 1e-4)
check('temperature at the current-only level 8 m = 28 - 0.2*8 degC', abs(Q['temperature'][int(np.argmin(np.abs(dq - 8)))] - (28 - 1.6 + 273.15)) < 1e-3)

print('2. options and errors')
P2 = vp.sample_profile([r_uv, r_ts], LON, LAT, time, add_surface=False)
check('add_surface=False: first level is the shallowest reader level', abs(P2['depth'][0] - 0.5) < 1e-6, str(P2['depth'][:2]))
P3 = vp.sample_profile([r_uv, r_ts], LON, LAT, time, z_max=10.0)
check('z_max=10 limits the profile', P3['depth'][-1] <= 10.0 + 1e-9, f"deepest {P3['depth'][-1]}")
try:
    vp.sample_profile([r_uv, r_ts], float(lon[LAND]), LAT, time); check('land column raises ValueError', False)
except ValueError as e:
    check('land column raises ValueError (position on land)', 'on land' in str(e), str(e)[:70])
try:
    vp.sample_profile([r_uv], LON, LAT, time); check('missing temperature raises ValueError', False)
except ValueError as e:
    check('missing temperature raises ValueError', 'sea_water_temperature' in str(e), str(e)[:60])
p_all = write('all.nc', dep_ts, {'uo': (field(dep_ts, lambda z: 0.10 + 0.01 * z), 'eastward_sea_water_velocity'),
                                 'vo': (field(dep_ts, lambda z: -0.05 + 0 * z), 'northward_sea_water_velocity'),
                                 'thetao': (field(dep_ts, lambda z: 28.0 - 0.2 * z), 'sea_water_potential_temperature'),
                                 'so': (field(dep_ts, lambda z: 36.0 + 0.01 * z), 'sea_water_salinity')})
R = vp.sample_profile(rg.Reader(p_all, name='all'), LON, LAT, time)      # a single reader, not a list
check('a single reader (not a list) with all four variables works and matches', len(R['depth']) == len(P['depth']) and np.allclose(R['temperature'], P['temperature']),
      f"{len(R['depth'])} levels")

print('3. TAMOC profile file (only if TAMOC is importable)')
try:
    import tamoc  # noqa: F401
    out = os.path.join(tmp, 'profile.nc')
    path = vp.create_profile(out, [r_uv, r_ts], LON, LAT, time)
    from tamoc import ambient
    prof = ambient.Profile(path, chem_names='all')
    T, S, Pr = prof.get_values(10.0, ['temperature', 'salinity', 'pressure'])
    check('create_profile file loads in tamoc.ambient.Profile', os.path.exists(path))
    check('profile at 10 m: T 299.15 K, S 36.10, p ~ 1 atm + 10 m of water', abs(T - 299.15) < 0.05 and abs(S - 36.10) < 0.01 and 190000 < Pr < 210000,
          f'T {T:.2f} S {S:.2f} p {Pr:.0f}')
except ImportError:
    print('  SKIP  TAMOC not importable here: the TAMOC writer is tested on Ixachi (sbatch), not locally')

print('\n' + ('ALL PASSED' if not FAILS else 'FAILED: ' + ', '.join(FAILS)))
sys.exit(1 if FAILS else 0)
