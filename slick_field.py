# Ported from PARCELS_GPU/parcels_gpu/slick.py (same author, 2026-10-16; after Multiphase_Plume/eulerian_slick.py, 2026-08-21) for CICOILv2, 2026-10-26.
# Unchanged except this header and the added helper functions at the end (emulsion_volume, sample); the NumPy path is the one used by CICOILv2 (`xp` is numpy).
"""Gridded surface-slick thickness on top of the Lagrangian particles (NumPy / CuPy / JAX): the "true surface expression" of the slick.

Model (after Multiphase_Plume/eulerian_slick.py, 2026-08-21):   dh/dt + div(h u) = D0 laplacian(h^2),   h(x, y, t) in m, one value per cell.
  * u = surface drift (current + windage + Stokes drift), sampled on the grid with drift3d.velocity3 (the same velocity the particles feel).
  * Spreading is the porous-medium equation with m = 2: D = 2 D0 h, so a thick patch spreads fast and a thin film slowly, and neighbouring patches coalesce
    without adding their thicknesses.  D0 is calibrated against Fay's gravity-viscous regime II (R^4 = 32 V D0 t / pi):
    D0 = pi k2^4 (delta g)^(2/3) V^(1/3) nu^(-1/3) / 32, k2 = 1.45, V = a characteristic volume (total surface volume at the start).  D0 keeps a residual V^(1/3)
    dependence (Multiphase_Plume documents it: a factor 10 in V moves the radius by about 20%).
  * Advection: flux form, second-order MUSCL with the minmod limiter and Heun time stepping (the first-order upwind scheme would add a numerical
    diffusivity of u dx / 2, about the size of the physical K_h); sub-stepped for CFL <= 0.4 (and the explicit stability limit of the spreading).
  * Coupling with the particles (one-way: the particles are not changed): the slick is made of the oil volume of the particles that are on the surface.
    Each step h is multiplied by the ratio  B / A  (A: emulsion volume of the particles on the surface at the start of the step; B: the same particles after weathering,
    evaporation and entrainment, i.e. a spatially uniform loss), moved and spread, and the volume of the particles that reach the surface in the step (release,
    resurfacing droplets) is added at their position with cloud-in-cell weights.  Volume that reaches land cells is removed (stranded) and counted.
Diagnostics: area by Bonn Agreement thickness class, maximum thickness, extent and width.
"""
import numpy as np

G = 9.81
K2 = 1.45
KX_M, KY_M = 111320.0, 110574.0            # metres per degree (longitude scaled by cos(lat0))
# Bonn Agreement appearance classes, lower thickness bounds in m: sheen, rainbow, metallic, discontinuous true colour, true colour
BONN = (0.04e-6, 0.3e-6, 5.0e-6, 50.0e-6, 200.0e-6)
BONN_NAMES = ("sheen", "rainbow", "metallic", "discontinuous", "true colour")


def fay_D0(V, delta, nu_w, k2=K2):
    """Porous-medium coefficient D0 that makes the front of an isolated blob of volume V follow Fay's regime II."""
    return np.pi * k2 ** 4 * (delta * G) ** (2.0 / 3.0) * V ** (1.0 / 3.0) * nu_w ** (-1.0 / 3.0) / 32.0


class Grid:
    """Regular grid in local metres: x = (lon - lon0) kx, y = (lat - lat0) ky, cell (j, i) centred at ((i + 1/2) dx, (j + 1/2) dx)."""
    def __init__(self, lon0, lat0, dx, nx, ny):
        self.lon0, self.lat0, self.dx, self.nx, self.ny = lon0, lat0, dx, nx, ny
        self.kx = KX_M * np.cos(np.deg2rad(lat0 + 0.5 * ny * dx / KY_M)); self.ky = KY_M

    def centres(self):
        x = (np.arange(self.nx) + 0.5) * self.dx; y = (np.arange(self.ny) + 0.5) * self.dx
        X, Y = np.meshgrid(x, y)
        return X, Y, self.lon0 + X / self.kx, self.lat0 + Y / self.ky


def bin_sum(xp, idx, w, size):
    if xp.__name__.startswith("jax"): return xp.zeros(size, dtype=w.dtype).at[idx].add(w)
    return xp.bincount(idx, weights=w, minlength=size).astype(w.dtype)


def deposit(xp, grid, lon, lat, vol):
    """Cloud-in-cell deposit of volumes (m3) of points at (lon, lat): thickness increment (ny, nx) in m, and the volume that fell outside the grid."""
    fx = (lon - grid.lon0) * grid.kx / grid.dx - 0.5; fy = (lat - grid.lat0) * grid.ky / grid.dx - 0.5
    i0 = xp.floor(fx); j0 = xp.floor(fy); wx = fx - i0; wy = fy - j0
    i0 = i0.astype("int32"); j0 = j0.astype("int32")
    out = xp.zeros(grid.nx * grid.ny, dtype=vol.dtype); inside = vol * 0.0
    for di, dj, w in ((0, 0, (1 - wx) * (1 - wy)), (1, 0, wx * (1 - wy)), (0, 1, (1 - wx) * wy), (1, 1, wx * wy)):
        i, j = i0 + di, j0 + dj
        ok = (i >= 0) & (i < grid.nx) & (j >= 0) & (j < grid.ny)
        out = out + bin_sum(xp, xp.clip(j, 0, grid.ny - 1) * grid.nx + xp.clip(i, 0, grid.nx - 1), xp.where(ok, w * vol, 0.0), grid.nx * grid.ny)
        inside = inside + xp.where(ok, w * vol, 0.0)
    return out.reshape(grid.ny, grid.nx) / (grid.dx * grid.dx), xp.sum(vol) - xp.sum(inside)


def barenblatt_R(V, D0, t):
    """Front radius of an isolated blob of volume V at age t: R^4 = 32 V D0 t / pi (Fay regime II)."""
    return (32.0 * V * D0 * t / np.pi) ** 0.25


def deposit_lens(xp, grid, lon, lat, vol, D0, t_inj, m):
    """Deposit each volume as the Barenblatt lens of age t_inj (the similarity solution the spreading then propagates) instead of a point.  Multiphase_Plume found
    (2026-08-21) that a point deposit makes the peak thickness depend on the grid (0.48, 0.83, 1.88 mm as dx halved), while a lens converges (0.48, 0.53, 0.52).
    The lens is evaluated on a (2m+1)^2 stencil of cells around the position and rescaled to conserve the volume exactly; if its radius is smaller than a cell
    the result is the cloud-in-cell weights (the stencil then holds one to four cells).  Returns (thickness (ny, nx), volume outside the grid)."""
    X0 = (lon - grid.lon0) * grid.kx; Y0 = (lat - grid.lat0) * grid.ky
    ic = xp.floor(X0 / grid.dx).astype("int32"); jc = xp.floor(Y0 / grid.dx).astype("int32")
    R = barenblatt_R(xp.maximum(vol, 1e-30), D0, t_inj)
    s = (D0 * t_inj) ** -0.5; kap = 1.0 / 16.0
    idxs, ws, oks = [], [], []
    for dj in range(-m, m + 1):
        for di in range(-m, m + 1):
            i, j = ic + di, jc + dj
            r2 = ((i + 0.5) * grid.dx - X0) ** 2 + ((j + 0.5) * grid.dx - Y0) ** 2
            w = xp.maximum(s * (kap * R ** 2 * s - kap * r2 * s), 0.0)
            idxs.append(xp.clip(j, 0, grid.ny - 1) * grid.nx + xp.clip(i, 0, grid.nx - 1)); ws.append(w)
            oks.append((i >= 0) & (i < grid.nx) & (j >= 0) & (j < grid.ny))
    wsum = sum(ws)
    # cloud-in-cell fallback where the lens falls between cell centres (wsum == 0): weights of the four cells around the point
    fx = X0 / grid.dx - 0.5; fy = Y0 / grid.dx - 0.5
    i0 = xp.floor(fx); j0 = xp.floor(fy); wx = fx - i0; wy = fy - j0; i0 = i0.astype("int32"); j0 = j0.astype("int32")
    out = xp.zeros(grid.nx * grid.ny, dtype=vol.dtype); inside = vol * 0.0
    use_lens = wsum > 0
    for idx, w, ok in zip(idxs, ws, oks):
        wt = xp.where(use_lens, w / xp.maximum(wsum, 1e-30), 0.0)
        out = out + bin_sum(xp, idx, xp.where(ok, wt * vol, 0.0), grid.nx * grid.ny); inside = inside + xp.where(ok, wt * vol, 0.0)
    for di, dj, w in ((0, 0, (1 - wx) * (1 - wy)), (1, 0, wx * (1 - wy)), (0, 1, (1 - wx) * wy), (1, 1, wx * wy)):
        i, j = i0 + di, j0 + dj
        ok = (i >= 0) & (i < grid.nx) & (j >= 0) & (j < grid.ny)
        wt = xp.where(use_lens, 0.0, w)
        out = out + bin_sum(xp, xp.clip(j, 0, grid.ny - 1) * grid.nx + xp.clip(i, 0, grid.nx - 1), xp.where(ok, wt * vol, 0.0), grid.nx * grid.ny)
        inside = inside + xp.where(ok, wt * vol, 0.0)
    return out.reshape(grid.ny, grid.nx) / (grid.dx * grid.dx), xp.sum(vol) - xp.sum(inside)


def lens_stencil(V_max, D0, t_inj, dx):
    """Half-width m (cells) of the lens stencil for the largest volume."""
    return int(np.ceil(barenblatt_R(V_max, D0, t_inj) / dx)) + 1


def _minmod(xp, a, b):
    return xp.where(a * b > 0, xp.where(xp.abs(a) < xp.abs(b), a, b), 0.0)


def _tendency(xp, h, u, v, dx):
    """-div(h u) by second-order MUSCL / minmod fluxes; no flux through the outer boundary (zero normal velocity at the faces of the domain edge is NOT
    imposed: the flux out of the domain is allowed and the caller counts it as lost)."""
    hp = xp.pad(h, 1, mode="edge")
    sx = _minmod(xp, hp[1:-1, 1:-1] - hp[1:-1, :-2], hp[1:-1, 2:] - hp[1:-1, 1:-1])                  # slope in x of every cell
    sy = _minmod(xp, hp[1:-1, 1:-1] - hp[:-2, 1:-1], hp[2:, 1:-1] - hp[1:-1, 1:-1])
    hL = h + 0.5 * sx; hR = h - 0.5 * sx                                                              # values at the east and west face of each cell
    uf = 0.5 * (u[:, :-1] + u[:, 1:])                                                                 # velocity at the interior faces (ny, nx-1)
    Fx = uf * xp.where(uf > 0, hL[:, :-1], hR[:, 1:])
    Fx = xp.concatenate([xp.maximum(-u[:, :1], 0.0) * 0.0 + u[:, :1] * hR[:, :1] * (u[:, :1] < 0), Fx, u[:, -1:] * hL[:, -1:] * (u[:, -1:] > 0)], axis=1)
    hN = h + 0.5 * sy; hS = h - 0.5 * sy
    vf = 0.5 * (v[:-1, :] + v[1:, :])
    Fy = vf * xp.where(vf > 0, hN[:-1, :], hS[1:, :])
    Fy = xp.concatenate([v[:1, :] * hS[:1, :] * (v[:1, :] < 0), Fy, v[-1:, :] * hN[-1:, :] * (v[-1:, :] > 0)], axis=0)
    return -((Fx[:, 1:] - Fx[:, :-1]) + (Fy[1:, :] - Fy[:-1, :])) / dx


def advect(xp, h, u, v, dx, dt):
    """One Heun step (second-order strong-stability-preserving) of dh/dt = -div(h u)."""
    h1 = xp.maximum(h + dt * _tendency(xp, h, u, v, dx), 0.0)
    return xp.maximum(0.5 * h + 0.5 * (h1 + dt * _tendency(xp, h1, u, v, dx)), 0.0)


def spread(xp, h, D0, dx, dt):
    """One explicit step of dh/dt = D0 laplacian(h^2) (5-point stencil, zero gradient at the boundary)."""
    h2 = xp.pad(h * h, 1, mode="edge")
    lap = (h2[1:-1, :-2] + h2[1:-1, 2:] + h2[:-2, 1:-1] + h2[2:, 1:-1] - 4.0 * h2[1:-1, 1:-1]) / (dx * dx)
    return xp.maximum(h + dt * D0 * lap, 0.0)


def loop(xp, n, body, init):
    if xp.__name__.startswith("jax"):
        import jax
        return jax.lax.fori_loop(0, n, lambda i, c: body(c), init)
    for _ in range(n): init = body(init)
    return init


def nsub_needed(dt, dx, umax, D0, hmax, cfl=0.4, safety=0.2):
    """Number of equal sub-steps of dt that satisfy both the advective CFL and the explicit spreading limit dt <= safety dx^2 / (8 D0 hmax)."""
    n_adv = umax * dt / (cfl * dx); n_dif = dt * 8.0 * D0 * hmax / (safety * dx * dx)
    return int(np.ceil(max(n_adv, n_dif, 1.0)))


def step_field(xp, h, u, v, D0, dx, dt, nsub, land=None):
    """Advect (u, v in m/s on the grid) and spread h for dt seconds in nsub sub-steps.  `land`: boolean grid of land cells (velocity 0, volume removed
    at the end of the step).  Returns (h, stranded volume m3)."""
    if land is not None: u = xp.where(land, 0.0, u); v = xp.where(land, 0.0, v)
    ds = dt / nsub
    h = loop(xp, nsub, lambda c: spread(xp, advect(xp, c, u, v, dx, ds), D0, dx, ds), h)
    if land is None: return h, xp.sum(h) * 0.0
    return xp.where(land, 0.0, h), xp.sum(xp.where(land, h, 0.0)) * dx * dx


def diagnostics(xp, h, X, Y, dx):
    """Per-step slick numbers: volume (m3), max thickness (m), area (m2) of each Bonn class (>= lower bound), and the extent from the second moments of the
    volume of the visible slick (h >= sheen bound): length = 4 sigma_major, breadth = 4 sigma_minor (the full axes of an elliptical uniform patch), mean width =
    area / length.  Returns a vector of 3 + 5 + 4 numbers: [volume, hmax, A_sheen..A_true, length, breadth, mean_width, visible_area]."""
    vis = h >= BONN[0]
    w = xp.where(vis, h, 0.0); W = xp.maximum(xp.sum(w), 1e-30)
    mx, my = xp.sum(w * X) / W, xp.sum(w * Y) / W
    cxx, cyy, cxy = xp.sum(w * (X - mx) ** 2) / W, xp.sum(w * (Y - my) ** 2) / W, xp.sum(w * (X - mx) * (Y - my)) / W
    tr, det = cxx + cyy, cxx * cyy - cxy * cxy
    disc = xp.sqrt(xp.maximum(0.25 * tr * tr - det, 0.0))
    l1, l2 = xp.maximum(0.5 * tr + disc, 0.0), xp.maximum(0.5 * tr - disc, 0.0)
    length, breadth = 4.0 * xp.sqrt(l1), 4.0 * xp.sqrt(l2)
    area = xp.sum(vis) * dx * dx
    cls = [xp.sum(h >= b) * dx * dx for b in BONN]
    return xp.stack([xp.sum(h) * dx * dx, xp.max(h)] + cls + [length, breadth, area / xp.maximum(length, 1.0), area]).astype(xp.float32)


# --- from PARCELS_GPU/parcels_gpu/slick_couple.py (sample; emulsion volume rewritten for CICOILv2 element arrays) ---
def emulsion_volume(mass_oil, water_fraction, oil_density, sea_water_density):
    """m3 of oil + entrained water of every element: mass of oil / (1 - water) / emulsion density."""
    w = np.minimum(water_fraction, 0.95)
    dens = w * sea_water_density + (1.0 - w) * oil_density
    return mass_oil / (1.0 - w) / dens


def sample(grid, h, lon, lat, floor=1.0e-9):
    """Thickness of the grid cell that holds each element (nearest cell); `floor` for the elements outside the grid."""
    i = np.floor((lon - grid.lon0) * grid.kx / grid.dx).astype("int32"); j = np.floor((lat - grid.lat0) * grid.ky / grid.dx).astype("int32")
    ok = (i >= 0) & (i < grid.nx) & (j >= 0) & (j < grid.ny)
    return np.where(ok, np.maximum(h[np.clip(j, 0, grid.ny - 1), np.clip(i, 0, grid.nx - 1)], floor), floor)
