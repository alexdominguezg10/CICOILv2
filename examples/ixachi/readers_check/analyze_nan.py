#!/usr/bin/env python3
"""Which elements keep status 0 but lose their mass (NaN)?  Read-only.  usage: python analyze_nan.py FILE.nc
Groups elements by origin_marker (one per seed_elements call = one TAMOC particle class, or the dissolved batch) and by is_bubble / is_dissolved."""
import sys
import numpy as np, xarray as xr
ds = xr.open_dataset(sys.argv[1])
t = ds['time'].values; n = ds.sizes['trajectory']
mass = ds['mass_oil'].values; st = ds['status'].values
om = ds['origin_marker'].values[:, 0] if ds['origin_marker'].ndim == 2 else ds['origin_marker'].values
def first(v): return v[np.isfinite(v).argmax()] if np.isfinite(v).any() else np.nan
ib = np.array([first(ds['is_bubble'].values[i]) if 'is_bubble' in ds else np.nan for i in range(n)])
idis = np.array([first(ds['is_dissolved'].values[i]) if 'is_dissolved' in ds else np.nan for i in range(n)])
diam = np.array([first(ds['diameter'].values[i]) if 'diameter' in ds else np.nan for i in range(n)])
m0 = np.array([first(mass[i]) for i in range(n)])
last = len(t) - 1
print(f'{n} trajectories, {len(t)} output steps; last step {str(t[-1])[:16]}')
print('\norigin_marker | n | bubble | dissolved | diameter(mm) | mass_oil at first appearance (kg) | finite mass at end | status0 at end | status0 & NaN mass at end | never appear')
for m in np.unique(om[np.isfinite(om)]):
    sel = om == m
    fin_end = np.isfinite(mass[sel, last]); s0 = (st[sel, last] == 0); never = ~np.isfinite(mass[sel]).any(axis=1)
    print(f'{int(m):6d} | {int(sel.sum()):3d} | {np.nanmax(ib[sel]):.0f} | {np.nanmax(idis[sel]):.0f} | {np.nanmedian(diam[sel]) * 1e3:7.2f} | {np.nansum(m0[sel]):12.0f} | '
          f'{int(fin_end.sum()):3d} | {int(s0.sum()):3d} | {int((s0 & ~fin_end).sum()):3d} | {int(never.sum())}')
# the odd ones: status 0 but NaN mass at the end - what else is NaN for them?
odd = np.where((st[:, last] == 0) & ~np.isfinite(mass[:, last]))[0]
print(f'\nelements with status 0 and NaN mass_oil at the last step: {len(odd)}')
if len(odd):
    for v in ('lon', 'lat', 'z', 'density', 'viscosity', 'water_fraction', 'mass_evaporated', 'mass_dispersed', 'moving'):
        if v in ds:
            a = ds[v].values[odd, last]; print(f'  {v:16s} NaN for {int((~np.isfinite(a)).sum())}/{len(odd)}')
    i = odd[0]
    print(f'  example trajectory {i}: origin_marker {om[i]}, bubble {ib[i]}, dissolved {idis[i]}, diameter {diam[i] * 1e3:.2f} mm')
    print('   mass_oil over time :', np.round(mass[i], 1).tolist())
    print('   status over time   :', st[i].tolist())
    for v in ('z', 'mass_evaporated', 'mass_dispersed'):
        if v in ds: print(f'   {v:16s}:', np.round(ds[v].values[i], 2).tolist())
