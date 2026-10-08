#!/usr/bin/env python3
"""
Unit test for readers/reader_croco_native.py on a tiny synthetic CROCO-like file (no real model data needed).

u (on xi_u points) increases linearly with the xi index and v (on eta_v points) with the eta index, so the
exact de-staggered value at a rho node is known: u(node i) = i - 0.5, v(node j) = j - 0.5 (interior nodes).
  python test_croco_native_reader.py
"""
import os, sys, tempfile, logging
from datetime import datetime
import numpy as np, xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'readers'))
from reader_croco_native import Reader

logging.disable(logging.WARNING)
NT, NS, NE, NX = 2, 4, 6, 8
lon1 = -93.0 + 0.1 * np.arange(NX); lat1 = 18.0 + 0.1 * np.arange(NE)
lon, lat = np.meshgrid(lon1, lat1)
s = (np.arange(NS) - NS + 0.5) / NS; Cs = s.copy()
i_u = np.arange(NX - 1, dtype='f4'); j_v = np.arange(NE - 1, dtype='f4')
u = np.broadcast_to(i_u, (NT, NS, NE, NX - 1)).copy()
v = np.broadcast_to(j_v[:, None], (NT, NS, NE - 1, NX)).copy()
mask = np.ones((NE, NX), 'f4'); mask[0, :] = 0           # southern row is land
ds = xr.Dataset({
    'u': (('time', 's_rho', 'eta_rho', 'xi_u'), u, {'units': 'm s-1'}),
    'v': (('time', 's_rho', 'eta_v', 'xi_rho'), v, {'units': 'm s-1'}),
    'zeta': (('time', 'eta_rho', 'xi_rho'), np.zeros((NT, NE, NX), 'f4')),
    'temp': (('time', 's_rho', 'eta_rho', 'xi_rho'), np.full((NT, NS, NE, NX), 20.0, 'f4')),
    'salt': (('time', 's_rho', 'eta_rho', 'xi_rho'), np.full((NT, NS, NE, NX), 36.0, 'f4')),
    'h': (('eta_rho', 'xi_rho'), np.full((NE, NX), 100.0, 'f4')),
    'mask_rho': (('eta_rho', 'xi_rho'), mask),
    'angle': (('eta_rho', 'xi_rho'), np.zeros((NE, NX), 'f4')),
    'lon_rho': (('eta_rho', 'xi_rho'), lon), 'lat_rho': (('eta_rho', 'xi_rho'), lat),
    's_rho': ('s_rho', s), 'Cs_rho': ('s_rho', Cs), 'hc': ((), np.float32(10.0)), 'Vtransform': ((), np.float32(2.0)),
    'time': ('time', np.array([0.0, 3600.0]), {'units': 'seconds since 2026-10-07 00:00:00'}),
})
path = os.path.join(tempfile.mkdtemp(), 'croco_synth.nc'); ds.to_netcdf(path)

t = datetime(2026, 10, 7, 0, 30)
i0, j0 = 4, 3                                             # interior rho node (wet)
lo, la = np.array([lon1[i0]]), np.array([lat1[j0]])
req = ['x_sea_water_velocity', 'y_sea_water_velocity']

def sample(reader):
    e, _ = reader.get_variables_interpolated(req, time=t, lon=lo, lat=la, z=np.zeros(1))
    return float(e[req[0]][0]), float(e[req[1]][0])

ok = True
def check(name, cond, msg=''):
    global ok; ok &= bool(cond); print(('PASS ' if cond else 'FAIL ') + name + (' ' + msg if msg else ''))

r = Reader(path, name='synthetic')
check('reader opens, destaggered flag set', r.destaggered)
check('mask_u derived from mask_rho, shape (NE, NX-1)', tuple(r.mask_u.shape) == (NE, NX - 1))
check('mask_v derived from mask_rho, shape (NE-1, NX)', tuple(r.mask_v.shape) == (NE - 1, NX))
check('mask_u is 0 next to land row, 1 inside', r.mask_u.values[0].sum() == 0 and r.mask_u.values[1].sum() == NX - 1)
uu, vv = sample(r)
check('de-staggered u at rho node = i-0.5', abs(uu - (i0 - 0.5)) < 1e-4, f'got {uu:.4f}, expected {i0 - 0.5}')
check('de-staggered v at rho node = j-0.5', abs(vv - (j0 - 0.5)) < 1e-4, f'got {vv:.4f}, expected {j0 - 0.5}')
us, vs = sample(Reader(path, name='stock', destagger=False))
check('destagger=False reproduces the stock half-cell offset', abs(us - uu) > 0.1 or abs(vs - vv) > 0.1, f'stock u,v = {us:.3f},{vs:.3f}')
print('ALL PASSED' if ok else 'FAILURES'); sys.exit(0 if ok else 1)
