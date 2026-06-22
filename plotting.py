"""
CICOILv2.0 — Plotting Routines
================================
Diagnostic and publication-quality figures for CICOILv2.0 simulations.

Available functions
-------------------
plot_mass_budget(budget, ...)
    Stacked area chart of oil fate fractions over time.

plot_trajectory_map(nc_file, ...)
    Particle trajectory map overlaid on coastline.

plot_photoox_timeseries(budget, ...)
    Photooxidation + evaporation time series, with and without UV.

plot_aromatic_budget(o, ...)
    Per-aromatic-group depletion (evaporation + photooxidation) bar chart.

plot_surface_fraction_map(nc_file, ...)
    Gridded particle density / surface oil mass at a chosen timestep.

plot_f_surface_sensitivity(...)
    f_surface sensitivity to Hs and Tp.

compare_weathering_modes(budgets, labels, ...)
    Side-by-side comparison of NOAA vs CICESE vs CICESE+photoox budgets.

Usage example
-------------
>>> from cicoilv2.plotting import plot_mass_budget, plot_trajectory_map
>>> import matplotlib.pyplot as plt
>>> o.run(...)
>>> budget = o.get_oil_budget()
>>> fig = plot_mass_budget(budget, title='Nohoch Alfa Jul 2023', filename='budget.png')
>>> fig = plot_trajectory_map('nohoch_cicese_24h.nc', spill_lon=-92.1, spill_lat=19.7)
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Optional, Union, List, Dict

import numpy as np

logger = logging.getLogger(__name__)

# ── Style defaults ─────────────────────────────────────────────────────────────
_COLORS = {
    'surface':      '#1565c0',   # blue
    'evaporated':   '#90caf9',   # light blue
    'photooxidized':'#f57f17',   # amber
    'dispersed':    '#37474f',   # dark grey
    'submerged':    '#0d47a1',   # dark blue
    'biodegraded':  '#4a148c',   # purple
    'stranded':     '#212121',   # near-black
    'gas':          '#ff8f00',   # orange
    'dissolved':    '#558b2f',   # green
}

_LABELS = {
    'surface':      'At surface',
    'evaporated':   'Evaporated',
    'photooxidized':'Photooxidized',
    'dispersed':    'Dispersed',
    'submerged':    'Submerged',
    'biodegraded':  'Biodegraded',
    'stranded':     'Stranded',
    'gas':          'Released gas',
    'dissolved':    'Dissolved (subsea)',
}

# ── Mass budget ────────────────────────────────────────────────────────────────

def plot_mass_budget(
    budget: dict,
    title: str = 'CICOILv2.0 — Oil fate budget',
    timestep_hours: float = 2.0,
    unit: str = 'tonnes',
    show_photoox: bool = True,
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (12, 5),
    dpi: int = 150,
):
    """Stacked-area mass budget plot.

    Parameters
    ----------
    budget : dict
        Output of ``OpenCiceseOil.get_oil_budget()``.
    title : str
        Figure title.
    timestep_hours : float
        Hours between output steps (used for x-axis).
    unit : str
        ``'tonnes'`` (default) or ``'kg'``.
    show_photoox : bool
        Include photooxidation wedge if present in budget.
    filename : str or Path, optional
        Save figure to file if provided.
    figsize, dpi : tuple, int
        Figure dimensions and resolution.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    scale = 1e-3 if unit == 'tonnes' else 1.0
    unit_label = 't' if unit == 'tonnes' else 'kg'

    initial = budget['mass_total'][0]
    n_steps = len(budget['mass_total'])
    t = np.arange(n_steps) * timestep_hours

    # Build wedge order (bottom → top)
    keys = ['dispersed', 'submerged', 'dissolved', 'gas',
            'biodegraded', 'surface', 'evaporated']
    if show_photoox and 'mass_photooxidized' in budget:
        keys.append('photooxidized')
    if 'mass_stranded' in budget and np.any(budget['mass_stranded'] > 0):
        keys.append('stranded')

    stacks = []
    colors = []
    labels = []
    for k in keys:
        key_full = f'mass_{k}'
        if key_full in budget and np.any(np.asarray(budget[key_full]) > 0):
            stacks.append(np.asarray(budget[key_full]) * scale)
            colors.append(_COLORS.get(k, '#999'))
            labels.append(_LABELS.get(k, k))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    fig.suptitle(title, fontsize=11, fontweight='bold')

    # Left: stacked area (mass)
    if stacks:
        ax1.stackplot(t, *stacks, colors=colors, labels=labels, alpha=0.88)
    ax1.set_xlabel('Time (h)')
    ax1.set_ylabel(f'Mass ({unit_label})')
    ax1.set_title('Oil fate — mass')
    ax1.set_xlim([0, t[-1]])
    ax1.set_ylim([0, initial * scale * 1.02])
    handles, lbls = ax1.get_legend_handles_labels()
    ax1.legend(handles[::-1], lbls[::-1], fontsize=8, loc='upper left',
               framealpha=0.8)
    ax1.grid(alpha=0.25)

    # Right: percentage area
    ax2.stackplot(
        t,
        *[100 * s / (initial * scale) for s in stacks],
        colors=colors, labels=labels, alpha=0.88,
    )
    ax2.set_xlabel('Time (h)')
    ax2.set_ylabel('Fraction (%)')
    ax2.set_title('Oil fate — fraction')
    ax2.set_xlim([0, t[-1]])
    ax2.set_ylim([0, 102])
    ax2.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=100, decimals=0))
    ax2.grid(alpha=0.25)

    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Budget plot saved → {filename}")
    return fig


# ── Trajectory map ─────────────────────────────────────────────────────────────

def plot_trajectory_map(
    nc_file: Union[str, Path],
    spill_lon: float,
    spill_lat: float,
    title: str = 'CICOILv2.0 — Particle trajectories',
    time_index: int = -1,
    color_by: str = 'status',
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (10, 8),
    dpi: int = 150,
):
    """Plot particle trajectories over a Cartopy map.

    Parameters
    ----------
    nc_file : str or Path
        OpenDrift output NetCDF.
    spill_lon, spill_lat : float
        Spill source location (marked with a star).
    title : str
        Figure title.
    time_index : int
        Timestep to show particle positions (default -1 = final).
    color_by : str
        ``'status'`` (active/stranded/evaporated) or ``'z'`` (depth).
    filename : str or Path, optional
        Save to file.
    figsize, dpi : tuple, int

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import xarray as xr

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        HAS_CARTOPY = True
    except ImportError:
        HAS_CARTOPY = False
        warnings.warn("Cartopy not available — falling back to plain lon/lat plot.")

    ds = xr.open_dataset(nc_file)
    lon = ds['lon'].values      # (trajectory, time)
    lat = ds['lat'].values
    z   = ds['z'].values
    status = ds['status'].values

    # All positions for trail
    lon_flat = lon.ravel()
    lat_flat = lat.ravel()

    # Final positions
    lon_final = lon[:, time_index]
    lat_final = lat[:, time_index]
    st_final  = status[:, time_index]
    z_final   = z[:, time_index]
    ds.close()

    # Domain
    pad = 2.0
    lon_min = np.nanmin(lon_flat) - pad
    lon_max = np.nanmax(lon_flat) + pad
    lat_min = np.nanmin(lat_flat) - pad
    lat_max = np.nanmax(lat_flat) + pad

    if HAS_CARTOPY:
        proj = ccrs.PlateCarree()
        fig = plt.figure(figsize=figsize)
        ax  = fig.add_subplot(111, projection=proj)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=proj)
        ax.add_feature(cfeature.LAND, facecolor='#e0ddd7', zorder=2)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.7, zorder=3)
        ax.add_feature(cfeature.BORDERS, linewidth=0.4, linestyle=':', zorder=3)
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray',
                          alpha=0.5, linestyle='--')
        gl.top_labels = False; gl.right_labels = False
        transform = ccrs.PlateCarree()
    else:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim([lon_min, lon_max])
        ax.set_ylim([lat_min, lat_max])
        ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')
        transform = None

    # Particle trails (thin, alpha)
    plot_kw = dict(transform=transform) if HAS_CARTOPY else {}
    for i in range(0, lon.shape[0], max(1, lon.shape[0] // 500)):
        ax.plot(lon[i, :], lat[i, :], '-', lw=0.3, alpha=0.15, color='#1565c0',
                **plot_kw)

    # Final positions coloured by status
    status_colors = {0: '#1565c0', 1: '#b71c1c', 2: '#f9a825',
                     3: '#37474f', 4: '#4a148c'}
    status_labels = {0: 'Active', 1: 'Stranded', 2: 'Evaporated',
                     3: 'Dispersed', 4: 'Biodegraded'}
    for st, col in status_colors.items():
        mask = st_final == st
        if mask.sum() > 0:
            ax.scatter(lon_final[mask], lat_final[mask],
                       s=6, c=col, alpha=0.7, label=status_labels[st],
                       zorder=4, **plot_kw)

    # Spill source
    ax.scatter([spill_lon], [spill_lat], s=200, c='red', marker='*',
               zorder=10, label='Spill source', **plot_kw)

    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='lower right', framealpha=0.85)

    plt.tight_layout()
    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Trajectory map saved → {filename}")
    return fig


# ── Photooxidation time series ─────────────────────────────────────────────────

def plot_photoox_timeseries(
    budget_with: dict,
    budget_without: Optional[dict] = None,
    timestep_hours: float = 2.0,
    title: str = 'CICOILv2.0 — Photooxidation impact',
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (13, 5),
    dpi: int = 150,
):
    """Photooxidation time series — optionally compared to a run without UV.

    Parameters
    ----------
    budget_with : dict
        Budget from a simulation with ``processes:photooxidation = True``.
    budget_without : dict, optional
        Budget from the same simulation with photooxidation off.
    timestep_hours : float
        Output interval in hours.
    title, filename, figsize, dpi : various
        Standard figure parameters.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    n = len(budget_with['mass_total'])
    t = np.arange(n) * timestep_hours
    initial = budget_with['mass_total'][0]
    scale = 1e-3  # kg → t

    phox = np.asarray(budget_with.get('mass_photooxidized',
                                      np.zeros(n))) * scale
    evap = np.asarray(budget_with['mass_evaporated']) * scale
    surf = np.asarray(budget_with['mass_surface']) * scale

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    fig.suptitle(title, fontsize=11, fontweight='bold')

    # Left: absolute masses
    ax = axes[0]
    ax.plot(t, surf, '-', color=_COLORS['surface'],      lw=2, label='Surface')
    ax.plot(t, evap, '-', color=_COLORS['evaporated'],   lw=2, label='Evaporated')
    ax.plot(t, phox, '-', color=_COLORS['photooxidized'],lw=2, label='Photooxidized')
    if budget_without:
        surf2 = np.asarray(budget_without['mass_surface']) * scale
        evap2 = np.asarray(budget_without['mass_evaporated']) * scale
        ax.plot(t, surf2, '--', color=_COLORS['surface'],    lw=1.5,
                alpha=0.6, label='Surface (no UV)')
        ax.plot(t, evap2, '--', color=_COLORS['evaporated'], lw=1.5,
                alpha=0.6, label='Evaporated (no UV)')
    ax.set_xlabel('Time (h)'); ax.set_ylabel('Mass (t)')
    ax.set_title('Absolute mass')
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # Middle: photooxidized fraction
    ax = axes[1]
    phox_pct = 100 * np.asarray(
        budget_with.get('mass_photooxidized', np.zeros(n))) / (initial * 1e-3)
    ax.fill_between(t, phox_pct, alpha=0.6, color=_COLORS['photooxidized'],
                    label='Photooxidized %')
    ax.set_xlabel('Time (h)'); ax.set_ylabel('Fraction (%)')
    ax.set_title('Photooxidized fraction of initial')
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # Right: f_surface bar (informational — use stored value if available)
    ax = axes[2]
    times_day = t / 24.0
    # Synthetic diurnal UV pattern (for illustration)
    cos_z = np.cos(np.radians((t % 24 - 12) * 15))
    I_UV_synth = np.maximum(cos_z, 0) * 50.0
    ax.fill_between(t, I_UV_synth, alpha=0.55, color='#f9a825', label='UV (synthetic)')
    ax.set_xlabel('Time (h)'); ax.set_ylabel('UV irradiance (W m⁻²)')
    ax.set_title('Synthetic diurnal UV forcing')
    ax.legend(fontsize=8); ax.grid(alpha=0.25)
    # Mark local noon
    for noon_h in np.arange(12, t[-1]+1, 24):
        ax.axvline(noon_h, color='red', lw=0.8, alpha=0.4, linestyle='--')

    plt.tight_layout()
    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Photoox time series → {filename}")
    return fig


# ── Aromatic group budget ──────────────────────────────────────────────────────

def plot_aromatic_budget(
    o,
    title: str = 'CICOILv2.0 — Aromatic group depletion',
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (10, 5),
    dpi: int = 150,
):
    """Stacked bar showing evaporation + photooxidation per aromatic group.

    Requires a completed CICESE-mode simulation with
    ``processes:photooxidation = True``.

    Parameters
    ----------
    o : OpenCiceseOil
        Completed simulation object.
    title, filename, figsize, dpi : various

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    groups = [f'G{i}' for i in range(1, 9)]
    group_labels = [
        'G1\nbenzene', 'G2\ntoluene', 'G3\nxylenes',
        'G4\nnaphthalene', 'G5\nacenaph.',
        'G6\nphenanthrene', 'G7\nchrysene', 'G8\npyrene+',
    ]
    x = np.arange(8)

    evap_arom  = o.cicese_mass_balance['mass_evaporated'][:, 8:16].sum(axis=0) / 1000
    phox_arom  = o.cicese_mass_balance.get('mass_photooxidized',
                  np.zeros_like(evap_arom))
    if phox_arom.ndim == 2:
        phox_arom = phox_arom[:, 8:16].sum(axis=0) / 1000
    remain_arom = o.cicese_mass_balance['mass_components'][:, 8:16].sum(axis=0) / 1000

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    fig.suptitle(title, fontsize=11, fontweight='bold')

    # Left: absolute
    ax1.bar(x, remain_arom, label='Remaining', color='#1565c0', alpha=0.8)
    ax1.bar(x, evap_arom,   bottom=remain_arom,
            label='Evaporated', color=_COLORS['evaporated'], alpha=0.8)
    ax1.bar(x, phox_arom,   bottom=remain_arom + evap_arom,
            label='Photooxidized', color=_COLORS['photooxidized'], alpha=0.8)
    ax1.set_xticks(x); ax1.set_xticklabels(group_labels, fontsize=8)
    ax1.set_ylabel('Mass (t)'); ax1.set_title('Aromatic group fate')
    ax1.legend(fontsize=8); ax1.grid(axis='y', alpha=0.25)

    # Right: fraction of initial
    initial = evap_arom + phox_arom + remain_arom
    initial[initial == 0] = 1  # avoid div-by-zero
    ax2.bar(x, 100 * remain_arom / initial, label='Remaining',
            color='#1565c0', alpha=0.8)
    ax2.bar(x, 100 * evap_arom / initial,   bottom=100 * remain_arom / initial,
            label='Evaporated', color=_COLORS['evaporated'], alpha=0.8)
    ax2.bar(x, 100 * phox_arom / initial,
            bottom=100 * (remain_arom + evap_arom) / initial,
            label='Photooxidized', color=_COLORS['photooxidized'], alpha=0.8)
    ax2.set_xticks(x); ax2.set_xticklabels(group_labels, fontsize=8)
    ax2.set_ylabel('Fraction of initial (%)'); ax2.set_title('Aromatic depletion fractions')
    ax2.legend(fontsize=8); ax2.grid(axis='y', alpha=0.25)

    plt.tight_layout()
    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Aromatic budget → {filename}")
    return fig


# ── f_surface sensitivity ──────────────────────────────────────────────────────

def plot_f_surface_sensitivity(
    Hs_range: tuple = (0.1, 3.0),
    Tp_range: tuple = (2.0, 12.0),
    w_rise: float = 0.27,
    n: int = 50,
    title: str = 'Wave surface-exposure fraction (f_surface)',
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (7, 5),
    dpi: int = 150,
):
    """Heatmap of f_surface as a function of Hs and Tp.

    Parameters
    ----------
    Hs_range, Tp_range : tuple of (min, max)
        Ranges for significant wave height (m) and peak period (s).
    w_rise : float
        Rise velocity of surface slick (m s⁻¹). Default 0.27 m s⁻¹.
    n : int
        Grid resolution.
    title, filename, figsize, dpi : various

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    Hs_arr = np.linspace(*Hs_range, n)
    Tp_arr = np.linspace(*Tp_range, n)
    Hs_grid, Tp_grid = np.meshgrid(Hs_arr, Tp_arr)

    # f_surface formula (Delvigne & Sweeney 1988)
    g = 9.81
    H_rms  = Hs_grid / np.sqrt(2)
    lam    = g * Tp_grid**2 / (2 * np.pi)
    k_wave = 2 * np.pi / lam
    z_mix  = 4 * H_rms**2 * k_wave
    t_rise = z_mix / w_rise
    f_surf = Tp_grid / (Tp_grid + t_rise)
    f_surf[Hs_grid <= 0] = 1.0

    fig, ax = plt.subplots(figsize=figsize)
    cf = ax.contourf(Hs_arr, Tp_arr, f_surf, levels=20, cmap='RdYlGn')
    cs = ax.contour(Hs_arr, Tp_arr, f_surf,
                    levels=[0.5, 0.7, 0.8, 0.9, 0.95],
                    colors='k', linewidths=0.8)
    ax.clabel(cs, fmt='%.2f', fontsize=8)
    plt.colorbar(cf, ax=ax, label='f_surface')

    # Mark typical Cantarell July conditions
    ax.plot(0.54, 4.4, 'w*', ms=14, label='Cantarell Jul (ERA5 mean)', zorder=5)
    ax.legend(fontsize=9, loc='upper right')

    ax.set_xlabel('Significant wave height Hs (m)')
    ax.set_ylabel('Peak period Tp (s)')
    ax.set_title(f'{title}\n(w_rise = {w_rise} m s⁻¹)', fontsize=10)
    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"f_surface sensitivity → {filename}")
    return fig


# ── Gridded density ────────────────────────────────────────────────────────────

def plot_surface_density(
    nc_file: Union[str, Path],
    time_index: int = -1,
    spill_lon: Optional[float] = None,
    spill_lat: Optional[float] = None,
    bin_deg: float = 0.05,
    weight_var: Optional[str] = 'mass_oil',
    title: str = 'CICOILv2.0 — Surface oil density',
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (10, 8),
    dpi: int = 150,
):
    """Gridded surface oil mass density map at a chosen timestep.

    Parameters
    ----------
    nc_file : str or Path
        OpenDrift output NetCDF.
    time_index : int
        Timestep to plot. Default -1 (final).
    spill_lon, spill_lat : float, optional
        Mark spill location with a star.
    bin_deg : float
        Grid cell size in degrees. Default 0.05° (~5 km).
    weight_var : str, optional
        Element variable to use as weights (e.g. ``'mass_oil'``).
        ``None`` for particle counts.
    title, filename, figsize, dpi : various

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import xarray as xr

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        HAS_CARTOPY = True
    except ImportError:
        HAS_CARTOPY = False

    ds  = xr.open_dataset(nc_file)
    lon = ds['lon'].values[:, time_index]
    lat = ds['lat'].values[:, time_index]
    st  = ds['status'].values[:, time_index]
    z   = ds['z'].values[:, time_index]

    weights = None
    if weight_var and weight_var in ds:
        weights = ds[weight_var].values[:, time_index]
    ds.close()

    # Keep surface active particles
    mask = (st == 0) & (z >= -0.1) & np.isfinite(lon) & np.isfinite(lat)
    lon_s = lon[mask]; lat_s = lat[mask]
    w_s   = weights[mask] / 1000 if weights is not None else None  # kg → t

    lon_bins = np.arange(lon_s.min() - bin_deg, lon_s.max() + 2*bin_deg, bin_deg)
    lat_bins = np.arange(lat_s.min() - bin_deg, lat_s.max() + 2*bin_deg, bin_deg)
    H, lon_e, lat_e = np.histogram2d(
        lon_s, lat_s, bins=[lon_bins, lat_bins], weights=w_s)

    cmap = plt.cm.YlOrRd.copy()
    cmap.set_under('none')
    norm = mcolors.LogNorm(vmin=max(H[H > 0].min(), 1e-3), vmax=H.max())

    if HAS_CARTOPY:
        proj = ccrs.PlateCarree()
        fig  = plt.figure(figsize=figsize)
        ax   = fig.add_subplot(111, projection=proj)
        ax.set_extent([lon_e[0], lon_e[-1], lat_e[0], lat_e[-1]], crs=proj)
        ax.add_feature(cfeature.LAND, facecolor='#d6d0c8', zorder=2)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8, zorder=3)
        gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5,
                          linestyle='--')
        gl.top_labels = False; gl.right_labels = False
        pm = ax.pcolormesh(
            lon_e[:-1] + bin_deg/2, lat_e[:-1] + bin_deg/2, H.T,
            cmap=cmap, norm=norm, transform=proj, zorder=4)
    else:
        fig, ax = plt.subplots(figsize=figsize)
        pm = ax.pcolormesh(
            lon_e[:-1] + bin_deg/2, lat_e[:-1] + bin_deg/2, H.T,
            cmap=cmap, norm=norm)
        ax.set_xlabel('Longitude (°E)'); ax.set_ylabel('Latitude (°N)')

    if spill_lon is not None:
        kw = dict(transform=ccrs.PlateCarree()) if HAS_CARTOPY else {}
        ax.scatter([spill_lon], [spill_lat], s=200, c='blue',
                   marker='*', zorder=10, **kw)

    label = f'Surface oil ({weight_var}, t cell⁻¹)' if weight_var else 'Particles cell⁻¹'
    plt.colorbar(pm, ax=ax, label=label, shrink=0.8, pad=0.02)
    ax.set_title(title, fontsize=11, fontweight='bold')
    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Density map → {filename}")
    return fig


# ── Weathering modes comparison ────────────────────────────────────────────────

def compare_weathering_modes(
    budgets: List[dict],
    labels: List[str],
    timestep_hours: float = 2.0,
    title: str = 'CICOILv2.0 — Weathering mode comparison',
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (14, 5),
    dpi: int = 150,
):
    """Side-by-side oil fate comparison for multiple weathering configurations.

    Typical use: NOAA vs CICESE vs CICESE+photoox.

    Parameters
    ----------
    budgets : list of dict
        One budget dict per configuration.
    labels : list of str
        Legend labels (same length as ``budgets``).
    timestep_hours : float
        Output interval.
    title, filename, figsize, dpi : various

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    n_conf = len(budgets)
    fates  = ['surface', 'evaporated', 'dispersed', 'biodegraded', 'photooxidized']
    line_styles = ['-', '--', '-.', ':']

    fig, axes = plt.subplots(1, len(fates), figsize=figsize, sharey=False)
    fig.suptitle(title, fontsize=11, fontweight='bold')

    for i, (bud, lbl) in enumerate(zip(budgets, labels)):
        n   = len(bud['mass_total'])
        t   = np.arange(n) * timestep_hours
        ini = bud['mass_total'][0]
        ls  = line_styles[i % len(line_styles)]

        for j, fate in enumerate(fates):
            ax  = axes[j]
            key = f'mass_{fate}'
            if key in bud:
                pct = 100 * np.asarray(bud[key]) / ini
                ax.plot(t, pct, ls, lw=2, color=_COLORS.get(fate, '#555'),
                        alpha=0.85, label=lbl)
            ax.set_xlabel('Time (h)')
            ax.set_ylabel('Fraction (%)')
            ax.set_title(_LABELS.get(fate, fate).capitalize())
            ax.legend(fontsize=7, loc='best')
            ax.grid(alpha=0.2)

    plt.tight_layout()
    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Mode comparison → {filename}")
    return fig


# ── UV / photoox rate constant reference plot ──────────────────────────────────

def plot_photoox_rate_constants(
    filename: Optional[Union[str, Path]] = None,
    figsize: tuple = (8, 4),
    dpi: int = 150,
):
    """Bar chart of CICOILv2.0 photooxidation rate constants for aromatic groups.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    K_PHOTO = np.array([2.9e-8, 2.0e-8, 5.0e-8, 2.0e-7,
                         3.5e-7, 1.0e-6, 2.5e-6, 5.0e-6])
    t_half  = np.log(2) / (K_PHOTO * 50.0 * 3600 * 8)  # days (8 h UV/day)

    group_labels = [
        'G1\nbenzene', 'G2\ntoluene', 'G3\nxylenes',
        'G4\nnaphthalene', 'G5\nacenaph.',
        'G6\nphenanthrene', 'G7\nchrysene', 'G8\npyrene+',
    ]
    x = np.arange(8)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    ax1.bar(x, K_PHOTO * 1e6, color=[plt.cm.Reds(0.3 + 0.7*i/7) for i in range(8)])
    ax1.set_xticks(x); ax1.set_xticklabels(group_labels, fontsize=8)
    ax1.set_ylabel('k (×10⁻⁶ m² J⁻¹)')
    ax1.set_title('Photoox rate constants')
    ax1.grid(axis='y', alpha=0.3)

    ax2.bar(x, t_half, color=[plt.cm.RdYlGn(0.8 - 0.7*i/7) for i in range(8)])
    ax2.set_xticks(x); ax2.set_xticklabels(group_labels, fontsize=8)
    ax2.set_ylabel('Half-life (days, 8 h UV/day)')
    ax2.set_title('Photoox half-lives')
    ax2.grid(axis='y', alpha=0.3)

    fig.suptitle('CICOILv2.0 — Aromatic photooxidation parameters', fontweight='bold')
    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=dpi, bbox_inches='tight')
        logger.info(f"Rate constants → {filename}")
    return fig
