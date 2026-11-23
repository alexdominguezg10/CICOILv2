#!/usr/bin/env python3
"""Read-only: print dimensions, variables and depth levels of the forcing files (xarray only, no OpenDrift).
usage: python inspect_forcing.py FILE [FILE ...]   or  python inspect_forcing.py DIR  (all *.nc in DIR, no recursion)"""
import sys, os, glob
import numpy as np, xarray as xr
paths = []
for a in sys.argv[1:]:
    paths += sorted(glob.glob(os.path.join(a, '*.nc'))) if os.path.isdir(a) else [a]
for p in paths:
    print('=' * 100); print(p, f'({os.path.getsize(p) / 1e6:.1f} MB)')
    ds = xr.open_dataset(p)
    print('dims:', dict(ds.sizes))
    for v in ds.data_vars:
        print(f'  {v:28s} dims={ds[v].dims} standard_name={ds[v].attrs.get("standard_name")} units={ds[v].attrs.get("units")}')
    for c in ('depth', 'z', 'lev', 'level', 'deptht', 'depthu'):
        if c in ds.coords or c in ds.variables:
            vals = np.asarray(ds[c].values).ravel()
            print(f'  depth coord {c}: n={vals.size}, values {np.round(vals[:12], 2)}{" ..." if vals.size > 12 else ""} attrs={dict(ds[c].attrs)}')
    for c in ('time', 'valid_time'):
        if c in ds.coords:
            t = ds[c].values; print(f'  {c}: {t.size} steps {str(t[0])[:16]} .. {str(t[-1])[:16]}')
    for c in ('longitude', 'lon', 'latitude', 'lat'):
        if c in ds.coords:
            v = ds[c].values; print(f'  {c}: {v.min():.3f}..{v.max():.3f} (n={v.size})')
    ds.close()
