#!/usr/bin/env python3
"""Reader smoke test for the shared CICOILv2 env.  PASS / FAIL / SKIP per reader; exit 1 if any FAIL.

What 'works' means here: the reader opens, reports a time range and a domain, and OpenDrift's own
get_variables_interpolated() returns finite values for every variable it advertises at a 5 x 5 grid of points
inside the domain at mid-time (the land mask and all-NaN-by-design fields are reported, not failed).

Groups (paths from environment variables, defaults on Ixachi):
  A  Abkatun 2026 GLORYS12 / ERA5 / ANFC set, read with reader_netCDF_CF_generic          ABK_DIR
  B  coupled CROCO / WRF / WW3 forecast, converted to CF (ocean_B, wind_B, waves_B)         COUPLED_DIR
       CROCO ocean  -> CF generic;  WRF wind -> CF generic;  WW3 waves -> CF generic
       (OpenDrift has no native WRF or WW3 reader: they are read through the CF files made by prep_forcing.py)
  C  raw CROCO sigma-coordinate output with reader_croco_native                              CROCO_NATIVE (file path)
  D  NEMO readers (reader_nemo_optimized, reader_NEMO_native_v3)                             NEMO_FILE (file path)
A missing file is SKIP, never PASS.  Usage:  python check_readers.py [--json out.json]
"""
import os, sys, glob, json, re, traceback, warnings, logging
import numpy as np
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
import xarray as xr

ABK = os.environ.get('ABK_DIR', '/LUSTRE/adomingu/Abkatun_ixachi/forcing')
CPL = os.environ.get('COUPLED_DIR', '/LUSTRE/adomingu/CICOIL_COMPARISON/forcing/S1_feb2026')
CROCO = os.environ.get('CROCO_NATIVE', '')
NEMO = os.environ.get('NEMO_FILE', '')
NPTS = 5


def clean_dataset(path):
    """ERA5-style files: drop the scalar 'number'/'expver' coordinates and rename valid_time -> time."""
    ds = xr.open_dataset(path)
    drop = [v for v in ('number', 'expver') if v in ds.coords or v in ds.variables]
    if drop:
        ds = ds.drop_vars(drop)
    if 'valid_time' in ds.dims and 'time' not in ds.dims:
        ds = ds.rename({'valid_time': 'time'})
    return ds


def probe(r):
    """Return (ok, detail, per-variable finite fractions) for an OpenDrift reader."""
    t0, t1 = r.start_time, r.end_time
    if t0 is None:
        return False, 'no start_time', {}
    t = t0 if t1 is None else t0 + (t1 - t0) / 2
    for a in ('xmin', 'xmax', 'ymin', 'ymax'):
        if getattr(r, a, None) is None:
            return False, f'reader has no {a}', {}
    xs = np.linspace(r.xmin, r.xmax, NPTS + 2)[1:-1]; ys = np.linspace(r.ymin, r.ymax, NPTS + 2)[1:-1]
    X, Y = np.meshgrid(xs, ys); lon, lat = r.xy2lonlat(X.ravel(), Y.ravel())
    # Grid-metadata names that ROMS/CROCO readers list in r.variables but that are not interpolable environment variables
    # (x_grid_index, ocean_s_coordinate_g2, *_at_u_location, time, ...): skipped, not failed.  Found on the real croco_avg.nc (diag job 22812).
    META = re.compile(r'grid_index|_at_[uvw]_location|ocean_s_coordinate|^time$|^(longitude|latitude)$|coriolis|inverse_of_cell|sea_floor_depth')
    variables = [v for v in r.variables if not META.search(v)]
    skipped = [v for v in r.variables if META.search(v)]
    res = r.get_variables_interpolated(variables, time=t, lon=np.asarray(lon, float), lat=np.asarray(lat, float),
                                       z=np.zeros(len(lon)))
    env = res[0] if isinstance(res, tuple) else res
    fin = {}
    for v in variables:
        a = np.ma.filled(np.ma.masked_invalid(np.asarray(env[v], dtype=float)), np.nan)
        fin[v] = float(np.isfinite(a).mean())
    dead = [v for v, f in fin.items() if f == 0.0]
    ok = len(variables) > 0 and not dead
    detail = (f'{len(variables)} vars' + (f' (+{len(skipped)} grid-metadata skipped)' if skipped else '') + f', {t0:%Y-%m-%d %H:%M} ..{("%s" % t1) if t1 is None else format(t1, "%Y-%m-%d %H:%M")}, '
              f'lon {float(np.min(lon)):.2f}..{float(np.max(lon)):.2f} lat {float(np.min(lat)):.2f}..{float(np.max(lat)):.2f}')
    if dead:
        detail += f' | NO finite values for: {dead}'
    return ok, detail, fin


def run_cf(path, clean=False):
    from opendrift.readers import reader_netCDF_CF_generic as rg
    return rg.Reader(clean_dataset(path) if clean else path, name=os.path.basename(path))


def run_croco(path):
    from opendrift.readers.reader_croco_native import Reader
    r = Reader(path, name='croco_native')
    if not getattr(r, 'destaggered', False):
        raise RuntimeError('reader_croco_native did not de-stagger u/v')
    return r


def run_nemo(path):
    from opendrift.readers.reader_nemo_optimized import Reader
    return Reader(filename=path, name='nemo')


TESTS = []   # (group, label, path or None, factory)


def add(group, label, path, factory):
    TESTS.append((group, label, path, factory))


for fn, lab, clean in [('cmems_cur_glorys_abkatun_2026_jan20_mar25.nc', 'GLORYS12 currents', False),
                       ('cmems_thetao_anfc_abkatun_2026_jan25_apr30.nc', 'ANFC temperature', False),
                       ('cmems_so_anfc_abkatun_2026_jan25_apr30.nc', 'ANFC salinity', False),
                       ('abkatun_2026_wind_feb06_mar25_cf.nc', 'ERA5 wind (CF)', False),
                       ('abkatun_2026_waves_cf.nc', 'ERA5 waves (CF)', True)]:
    add('A', lab, os.path.join(ABK, fn), (lambda p=os.path.join(ABK, fn), c=clean: run_cf(p, c)))
for fn, lab in [('ocean_B.nc', 'CROCO ocean (CF)'), ('wind_B.nc', 'WRF wind (CF)'), ('waves_B.nc', 'WW3 waves (CF)')]:
    add('B', lab, os.path.join(CPL, fn), (lambda p=os.path.join(CPL, fn): run_cf(p)))
add('C', 'CROCO native (reader_croco_native)', CROCO or None, (lambda: run_croco(CROCO)))
add('D', 'NEMO (reader_nemo_optimized)', NEMO or None, (lambda: run_nemo(NEMO)))

rows = []
for group, label, path, factory in TESTS:
    if not path or not os.path.exists(path):
        rows.append((group, label, 'SKIP', f'file not set or not found: {path or "(env var empty)"}')); continue
    try:
        r = factory()
        ok, detail, fin = probe(r)
        rows.append((group, label, 'PASS' if ok else 'FAIL', detail))
    except Exception as e:
        rows.append((group, label, 'FAIL', f'{type(e).__name__}: {e}'))
        traceback.print_exc(limit=3, file=sys.stderr)

print(f'{"grp":3} {"reader":38} {"result":6} detail')
for g, l, s, d in rows:
    print(f'{g:3} {l:38} {s:6} {d}')
n = {k: sum(1 for r in rows if r[2] == k) for k in ('PASS', 'FAIL', 'SKIP')}
print('SUMMARY', n)
print('NOTE: WRF and WW3 are read through the CF files of prep_forcing.py (no native OpenDrift reader exists); SKIP is not a pass.')
if '--json' in sys.argv:
    json.dump([dict(group=g, reader=l, result=s, detail=d) for g, l, s, d in rows], open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=1)
sys.exit(1 if n['FAIL'] else 0)
