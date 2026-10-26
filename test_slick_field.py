"""Checks for slick_field.py (ported from PARCELS_GPU test_slick.py, NumPy only):   python test_slick_field.py"""
import sys, numpy as np
import slick_field as S

be = "numpy"
xp = np
ok_all = True
host = np.asarray


def check(name, ok, info=""):
    global ok_all
    ok_all &= bool(ok); print(("PASS " if ok else "FAIL ") + name, info)


# 1. an isolated blob spreads as Fay's regime II (Multiphase_Plume verify_single_blob case): front radius against R^4 = 32 V D0 t / pi
V, delta, nu_w = 0.1, 0.241, 9.093e-7
D0 = S.fay_D0(V, delta, nu_w)
dx, L = 0.2, 60.0
n = int(L / dx)
grid = S.Grid(0.0, 0.0, dx, n, n)
X, Y, lon, lat = grid.centres()
lon0 = grid.lon0 + 0.5 * n * dx / grid.kx; lat0 = grid.lat0 + 0.5 * n * dx / grid.ky
t0 = 16.6
m = S.lens_stencil(V, D0, t0, dx)
h, out = S.deposit_lens(xp, grid, xp.asarray([lon0]), xp.asarray([lat0]), xp.asarray([V]), D0, t0, m)
V_in = float(host(xp.sum(h))) * dx * dx
check("1a lens deposit conserves the volume (%.6f of %.6f m3)" % (V_in, V), abs(V_in - V) / V < 1e-5)
t = t0; res = []; targets = [100.0, 200.0, 400.0, 600.0]
hmax_b = V / (dx * dx)
while t < 600.0:
    dt = min(60.0, 600.0 - t)
    hm = float(host(xp.max(h)))
    nsub = S.nsub_needed(dt, dx, 0.0, D0, hm)
    h, _ = S.step_field(xp, h, xp.zeros_like(h), xp.zeros_like(h), D0, dx, dt, nsub)
    t += dt
    if t >= targets[0] - 1e-9:
        targets.pop(0) if len(targets) > 1 else targets.__setitem__(0, 1e99)
        hh = host(h); r = np.hypot(X - 0.5 * n * dx, Y - 0.5 * n * dx)
        Rm = float(r[hh > 1e-3 * hh.max()].max()); Ra = S.barenblatt_R(V, D0, t)
        res.append((t, Rm, Ra))
errs = [abs(Rm - Ra) / Ra for t_, Rm, Ra in res]
check("1b front radius against Fay at t = %s s: relative errors %s" % (", ".join("%.0f" % r[0] for r in res), ", ".join("%.1f%%" % (100 * e) for e in errs)), max(errs) < 0.05)
V_end = float(host(xp.sum(h))) * dx * dx
check("1c volume after the march (relative change %.1e)" % ((V_end - V) / V), abs(V_end - V) / V < (1e-5 if xp is np else 1e-4))

# 2. translation in a uniform current: mass, centroid, positivity, numerical diffusion of the MUSCL/Heun scheme
dx = 10.0; n = 160
grid2 = S.Grid(0.0, 0.0, dx, n, n); X2, Y2, _, _ = grid2.centres()
sig = 40.0; x0, y0 = 400.0, 500.0
h0 = np.exp(-((X2 - x0) ** 2 + (Y2 - y0) ** 2) / (2 * sig ** 2)) * 1e-3
u0, v0, T = 0.5, 0.2, 1500.0
nsub = S.nsub_needed(T, dx, 0.6, 0.0, 0.0)
h2, _ = S.step_field(xp, xp.asarray(h0), xp.full(h0.shape, u0, dtype=xp.float32), xp.full(h0.shape, v0, dtype=xp.float32), 0.0, dx, T, nsub)
hh = host(h2).astype(float)
M0, M1 = h0.sum() * dx * dx, hh.sum() * dx * dx
cx, cy = (hh * X2).sum() / hh.sum(), (hh * Y2).sum() / hh.sum()
var = (hh * ((X2 - cx) ** 2)).sum() / hh.sum()
Kn = (var - sig ** 2) / (2 * T)
check("2a mass conserved in translation (relative change %.1e)" % ((M1 - M0) / M0), abs(M1 - M0) / M0 < 1e-5)
check("2b centroid moved by (%.1f, %.1f) m, expected (%.1f, %.1f) m" % (cx - x0, cy - y0, u0 * T, v0 * T), abs(cx - x0 - u0 * T) < 0.1 * dx and abs(cy - y0 - v0 * T) < 0.1 * dx)
check("2c no negative thickness (min %.1e)" % hh.min(), hh.min() >= 0.0)
check("2d numerical diffusivity of the advection %.3f m2/s (first-order upwind would add u dx / 2 = %.2f m2/s)" % (Kn, 0.5 * u0 * dx), Kn < 0.25 * 0.5 * u0 * dx)

# 3. the lens deposit makes the peak thickness converge with the grid (Multiphase_Plume: points gave 0.48, 0.83, 1.88 mm; lenses 0.48, 0.53, 0.52 mm)
V1, D01, tinj = 0.1, S.fay_D0(0.1, 0.241, 9.093e-7), 16.6
peak_lens, peak_pt = [], []
for dxx in (1.0, 0.5, 0.25):
    nn = int(40 / dxx); g = S.Grid(0.0, 0.0, dxx, nn, nn)
    lo = g.lon0 + 0.5 * nn * dxx / g.kx + 0.37 * dxx / g.kx; la = g.lat0 + 0.5 * nn * dxx / g.ky + 0.21 * dxx / g.ky
    hl, _ = S.deposit_lens(xp, g, xp.asarray([lo]), xp.asarray([la]), xp.asarray([V1]), D01, tinj, S.lens_stencil(V1, D01, tinj, dxx))
    hp, _ = S.deposit(xp, g, xp.asarray([lo]), xp.asarray([la]), xp.asarray([V1]))
    peak_lens.append(float(host(xp.max(hl))) * 1e3); peak_pt.append(float(host(xp.max(hp))) * 1e3)
check("3 peak thickness of a deposited blob at dx = 1, 0.5, 0.25 m: lens %s mm (converges), point %s mm" % (", ".join("%.2f" % p for p in peak_lens), ", ".join("%.2f" % p for p in peak_pt)),
      (max(peak_lens) - min(peak_lens)) / max(peak_lens) < 0.15 and peak_pt[2] > 2.5 * peak_pt[0])

# 4. diagnostics on an anisotropic Gaussian slick: length and breadth are 4 sigma of the major and minor axes, class areas are step functions of h
dx = 10.0; n = 200
g4 = S.Grid(0.0, 0.0, dx, n, n); X4, Y4, _, _ = g4.centres()
sx, sy = 400.0, 100.0
hg = 1.0e-3 * np.exp(-((X4 - 1000.0) ** 2 / (2 * sx ** 2) + (Y4 - 1000.0) ** 2 / (2 * sy ** 2)))
d = host(S.diagnostics(xp, xp.asarray(hg, dtype=xp.float32), xp.asarray(X4, dtype=xp.float32), xp.asarray(Y4, dtype=xp.float32), dx))
check("4a volume of the diagnostics = sum of h dx^2 (%.3f vs %.3f m3)" % (d[0], hg.sum() * dx * dx), abs(d[0] - hg.sum() * dx * dx) / d[0] < 1e-4)
vis = (hg >= S.BONN[0])
Lth = 4 * sx * 0.97; Bth = 4 * sy * 0.97
check("4b length %.0f m and breadth %.0f m of an ellipse with sigma 400 and 100 m (4 sigma = 1600 and 400, truncated by the visibility threshold)" % (d[7], d[8]), abs(d[7] / 1600 - 1) < 0.1 and abs(d[8] / 400 - 1) < 0.1)
a_true = (hg >= 200e-6).sum() * dx * dx
check("4c area of the thickest class: %.0f m2 (direct count %.0f m2)" % (d[6], a_true), abs(d[6] - a_true) < 1.0)

# 5. NumPy and JAX give the same field (only when run with jax: compare with a NumPy reference computed here)
if be == "jax":
    ns = S.nsub_needed(600.0, dx, 0.6, 1e-3, 1e-3)
    href, _ = S.step_field(np, h0.astype("f4"), np.full(h0.shape, u0, dtype="f4"), np.full(h0.shape, v0, dtype="f4"), 1e-3, dx, 600.0, ns)
    hjax, _ = S.step_field(xp, xp.asarray(h0, dtype=xp.float32), xp.full(h0.shape, u0, dtype=xp.float32), xp.full(h0.shape, v0, dtype=xp.float32), 1e-3, dx, 600.0, ns)
    err = float(np.max(np.abs(host(hjax) - href)) / np.max(href))
    check("5 JAX equals NumPy in %d stable sub-steps (max relative difference %.1e)" % (ns, err), err < 1e-4)

print("ALL PASS" if ok_all else "SOME FAILED"); sys.exit(0 if ok_all else 1)

