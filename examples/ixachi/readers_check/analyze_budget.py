#!/usr/bin/env python3
"""Where does the oil mass go?  Read-only analysis of one OpenDrift/CICOIL output file (xarray only).
usage: python analyze_budget.py FILE.nc
Prints per output step: elements by status, sum of mass_oil of the active elements, and what happened to the ones that left."""
import sys
import numpy as np, xarray as xr
ds = xr.open_dataset(sys.argv[1])
print('dims', dict(ds.sizes)); print('variables', [v for v in ds.data_vars][:40])
st = ds['status']
print('status attrs:', {k: v for k, v in st.attrs.items() if k in ('flag_values', 'flag_meanings')})
t = ds['time'].values
mass = ds['mass_oil'].values                    # (trajectory, time)
stv = st.values
n = ds.sizes['trajectory']
print(f'trajectories: {n}')
first_valid = np.array([np.argmax(np.isfinite(mass[i])) if np.isfinite(mass[i]).any() else -1 for i in range(n)])
last_valid = np.array([len(t) - 1 - np.argmax(np.isfinite(mass[i][::-1])) if np.isfinite(mass[i]).any() else -1 for i in range(n)])
print('elements never appearing (no finite mass at any output step):', int((first_valid < 0).sum()))
print('first output step at which elements appear: ', np.unique(first_valid, return_counts=True))
print('\ntime step | n finite mass | n status 0 | other status codes | sum mass_oil (kg) | sum mass_evaporated | sum mass_dispersed')
me = ds['mass_evaporated'].values if 'mass_evaporated' in ds else None
md = ds['mass_dispersed'].values if 'mass_dispersed' in ds else None
for k in range(len(t)):
    f = np.isfinite(mass[:, k]); s = stv[:, k]
    other = {int(v): int(c) for v, c in zip(*np.unique(s[np.isfinite(s) & (s != 0)], return_counts=True))}
    print(f'{str(t[k])[:16]} | {int(f.sum()):5d} | {int(np.nansum(s == 0)):5d} | {other} | {np.nansum(mass[:, k]):12.0f} | '
          f'{(np.nansum(me[:, k]) if me is not None else float("nan")):12.0f} | {(np.nansum(md[:, k]) if md is not None else float("nan")):10.0f}')
# elements that disappear before the last step: when, with what status and how much mass
gone = np.where((last_valid >= 0) & (last_valid < len(t) - 1))[0]
print(f'\nelements whose data end before the last output step: {len(gone)}')
if len(gone):
    lv = last_valid[gone]; m_last = np.array([mass[i, last_valid[i]] for i in gone]); s_last = np.array([stv[i, last_valid[i]] for i in gone])
    print('  by last valid step:', dict(zip(*[x.tolist() for x in np.unique(lv, return_counts=True)])))
    print('  last status code:', dict(zip(*[x.tolist() for x in np.unique(s_last[np.isfinite(s_last)], return_counts=True)])))
    print(f'  mass_oil of those elements at their last valid step: {np.nansum(m_last):.0f} kg')
for v in ('z', 'lon', 'lat'):
    if v in ds:
        a = ds[v].values[:, -1]; print(f'{v} at the last step: min {np.nanmin(a):.3f} max {np.nanmax(a):.3f}')
