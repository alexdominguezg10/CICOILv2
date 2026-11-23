#!/usr/bin/env python3
"""Diagnosis of reader_croco_native on a REAL croco_avg.nc (read-only).  Found by check_readers.py: KeyError 'x_grid_index'.
Prints the reader's domain, the dataset dimensions/variables, then tries every variable one by one at one wet point and at a 5x5 grid,
with FULL tracebacks, so the failing variable / code line is visible.   CROCO_NATIVE=/path/croco_avg.nc python diag_croco_native.py"""
import os, sys, traceback, warnings, logging
import numpy as np
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
p = os.environ['CROCO_NATIVE']
from opendrift.readers.reader_croco_native import Reader
import opendrift
print('opendrift', opendrift.__version__, '| file', p, flush=True)
r = Reader(p, name='croco_native')
print('variables:', list(r.variables))
print('time:', r.start_time, '..', r.end_time, '| nsteps', getattr(r, 'times', None) is not None and len(r.times))
print('xmin/xmax/ymin/ymax:', r.xmin, r.xmax, r.ymin, r.ymax, '| proj4:', getattr(r, 'proj4', None))
ds = getattr(r, 'Dataset', None)
if ds is not None:
    print('Dataset dims:', dict(ds.sizes))
    for v in ('u', 'v', 'temp', 'salt', 'zeta', 'mask_rho', 'h', 'lon_rho', 'lat_rho', 's_rho'):
        if v in ds:
            print(f'  {v:9s} dims={ds[v].dims} shape={ds[v].shape} chunks={getattr(ds[v].data, "chunks", None) is not None}')
    print('coords:', list(ds.coords))
t = r.start_time + (r.end_time - r.start_time) / 2
xs = np.linspace(r.xmin, r.xmax, 7)[1:-1]; ys = np.linspace(r.ymin, r.ymax, 7)[1:-1]
X, Y = np.meshgrid(xs, ys); lon, lat = r.xy2lonlat(X.ravel(), Y.ravel())
lon = np.asarray(lon, float); lat = np.asarray(lat, float)
print(f'probe lon {lon.min():.2f}..{lon.max():.2f} lat {lat.min():.2f}..{lat.max():.2f}', flush=True)
for npts, label in ((1, 'one point'), (len(lon), '5x5 points')):
    for v in list(r.variables):
        try:
            res = r.get_variables_interpolated([v], time=t, lon=lon[:npts], lat=lat[:npts], z=np.zeros(npts))
            env = res[0] if isinstance(res, tuple) else res
            a = np.asarray(env[v], float)
            print(f'OK   {label:11s} {v:32s} finite {np.isfinite(a).mean():.2f}')
        except Exception:
            print(f'FAIL {label:11s} {v}'); traceback.print_exc()
