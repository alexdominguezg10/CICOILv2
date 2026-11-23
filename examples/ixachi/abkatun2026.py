#!/usr/bin/env python3
"""Abkatun 2026 example for CICOILv2 (shared env).  Spill 2026-02-06 at 19.274 N, 92.236 W.

  --mode 2d : surface release (z = 0), drift:vertical_mixing OFF  -> oil stays on the surface (no wave entrainment)
  --mode 3d : subsea release at the source depth (default 19 m, ATL24 bathymetry), drift:vertical_mixing ON
              -> droplets rise, wave entrainment and vertical mixing act.  --droplets selects how the release is built:
                 tamoc (default)  the CICOIL v1 route, see below
                 default          CICOIL's own 0.5-5 mm uniform droplet range (seed_elements, no TAMOC)
                 lognormal:D50_mm:SIGMA   lognormal sizes (seed_elements, no TAMOC)

The CICOIL v1 route for TAMOC data (same steps as v1's examples/example_cicoil-tamoc.py, ported to OpenDrift 1.14.9):
  1. live-oil composition from the ADIOS record:  o.fluid_properties.get_live_composition(GOR)
  2. TAMOC ambient profile (T, S, u, v against depth) read from the OpenDrift readers at the source: vertical_profile.create_profile
     (v1 module vertical_profile.py, rewritten for OpenDrift 1.14.9 and added to CICOILv2 on branch feat/vertical-profile)
  3. bent plume model:  run_tamoc.run_tamoc(...) -> <name>_plume.nc, <name>_particles.nc
  4. Plume(plume_file).get_particle_properties(particles_file, tracked=False)
     (v1's example used tracked=True, which needs per-particle track files <name>_plume000.nc ...; TAMOC 4.1.2 does not write them
      (found on Ixachi, job 22816).  tracked=False is the path Plume supports without them: the dissolved mass is taken from the
      plume only and the far-field dissolution during the rise is not accounted for.)
  5. o.seed_plume_elements(lon, lat, plume, [t_start, t_end], number=N, z_uncertainty=...): one group of elements per TAMOC
     size class, number and mass proportional to the class mass flow, depth drawn around the plume exit depth, horizontal
     spread = plume width, seeded over the release interval.
  The TAMOC release parameters below (flow, GOR, hole diameter, jet temperature) are NOT recorded for Abkatun 2026: the flow, GOR
  and hole come from the Abkatun TAMOC sensitivity script; the jet temperature is ambient at the source + --tj-delta (assumption).

Honest limits (read before presenting the results):
  * At 19-30 m depth the TAMOC plume reaches the surface; the far field then starts from where TAMOC says the droplets leave the plume.
  * When a diameter is given, CICOILv2 keeps it constant for the whole run (seed_elements sets keep_droplet_diameter).
  * --forcing coupled (default): CROCO + WRF + WW3 coupled-forecast CF files (S1_feb2026, 1-24 Feb 2026, read in place).
    --forcing glorys: GLORYS12 + ERA5 + ANFC (inputs folder); these files have ONE depth level, so no vertical profile and no 3d TAMOC.
    The 2d/instantaneous runs use a total volume of 1000 m3 (not the real Abkatun volume; --volume-m3).
  * Mass per element = volume x density / N (the model does NOT use a fixed 185.6 kg); the script prints it.

Usage:  python abkatun2026.py --mode 2d|3d --forcing-dir DIR --out DIR [-n 2000] [--days 10] [--seed 1] [--dispersion-mode legacy|transport|budget]
"""
import argparse, os, sys, tempfile, warnings, logging
from datetime import datetime, timedelta
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--mode', required=True, choices=['2d', '3d'])
ap.add_argument('--forcing-dir', default=os.environ.get('ABK_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'inputs', 'abkatun2026')))
ap.add_argument('--out', required=True, help='output directory')
ap.add_argument('--forcing', choices=['coupled', 'glorys'], default='coupled',
                help='coupled = CROCO ocean + WRF wind + WW3 waves (CF files, 10 depth levels, needed for the TAMOC profile); '
                     'glorys = small GLORYS12/ANFC/ERA5 set (ONE depth level: no vertical profile, 3d tamoc cannot use it)')
ap.add_argument('--coupled-dir', default=os.environ.get('COUPLED_DIR', '/LUSTRE/adomingu/CICOIL_COMPARISON/forcing/S1_feb2026'))
ap.add_argument('-n', '--particles', type=int, default=2000)
ap.add_argument('--days', type=float, default=10.0)
ap.add_argument('--dt', type=int, default=30, help='time step, minutes')
ap.add_argument('--out-dt', type=int, default=6, help='output interval, hours')
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--loglevel', type=int, default=30, help='OpenDrift log level (20 = info: shows missing-data and deactivation messages)')
ap.add_argument('--volume-m3', type=float, default=1000.0, help='total oil volume released instantaneously')
ap.add_argument('--lon', type=float, default=-92.236); ap.add_argument('--lat', type=float, default=19.274)
ap.add_argument('--depth', type=float, default=19.0, help='3d: source depth in m (ATL24 18.5 m ~ 19; GEBCO ~30)')
ap.add_argument('--release-date', default='2026-02-06')
ap.add_argument('--dispersion-mode', default='legacy', choices=['legacy', 'transport', 'budget'])
ap.add_argument('--droplets', default='tamoc', help='3d only: tamoc | default | lognormal:D50_mm:SIGMA')
ap.add_argument('--gor', type=float, default=150.0, help='TAMOC gas-oil ratio, scf/bbl')
ap.add_argument('--q-bbl-day', type=float, default=10000.0, help='TAMOC oil flow rate, bbl/day')
ap.add_argument('--hole', type=float, default=0.10, help='TAMOC hole diameter, m')
ap.add_argument('--bins', type=int, default=12, help='TAMOC number of size bins for droplets and bubbles (v1 example: 12)')
ap.add_argument('--tj-delta', type=float, default=5.0, help='TAMOC jet temperature = ambient at the source + this (K); ASSUMPTION')
ap.add_argument('--z-unc', type=float, default=2.0, help='seed_plume_elements z_uncertainty (m); v1 example used 10 m for a 1500 m release')
ap.add_argument('--release-hours', type=float, default=24.0, help='3d tamoc: duration of the release (mass = flow x duration)')
ap.add_argument('--profile-lon', type=float, default=None, help='longitude where the TAMOC ambient profile is read (default: the source)')
ap.add_argument('--profile-lat', type=float, default=None, help='latitude where the TAMOC ambient profile is read (default: the source)')
ap.add_argument('--bio', action='store_true', help='biodegradation on'); ap.add_argument('--photo', action='store_true', help='photooxidation on')
a = ap.parse_args()
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
os.makedirs(a.out, exist_ok=True)
FD = os.path.abspath(a.forcing_dir)


# ------------------------------------------------------------------ forcing readers (same files as Abkatun_ixachi/run_abkatun_julia.py)
def make_readers_coupled():
    """CROCO ocean (T, S, u, v on 10 depth levels, 0..200 m), WRF wind and WW3 waves/Stokes drift, as CF files made by prep_forcing.py."""
    from opendrift.readers import reader_netCDF_CF_generic as rg
    readers = []
    for fn, nm in (('ocean_B.nc', 'croco'), ('wind_B.nc', 'wrf'), ('waves_B.nc', 'ww3')):
        pth = os.path.join(a.coupled_dir, fn)
        if not os.path.exists(pth):
            sys.exit(f'FAIL: {pth} not found (--coupled-dir)')
        readers.append(rg.Reader(pth, name=nm))
    return readers


def make_readers():
    import xarray as xr
    from opendrift.readers import reader_netCDF_CF_generic as rg
    files = dict(ocean='cmems_cur_glorys_abkatun_2026_jan20_mar25.nc', thetao='cmems_thetao_anfc_abkatun_2026_jan25_apr30.nc',
                 so='cmems_so_anfc_abkatun_2026_jan25_apr30.nc', wind='abkatun_2026_wind_feb06_mar25_cf.nc',
                 waves='abkatun_2026_waves_cf.nc', uv='abkatun_2026_era5_uv.nc')
    readers, missing = [], []
    for key, fn in files.items():
        p = os.path.join(FD, fn)
        if not os.path.exists(p):
            missing.append(fn); continue
        if key in ('waves', 'uv'):       # ERA5 files: scalar number/expver coordinates break the reader
            ds = xr.open_dataset(p)
            drop = [v for v in ('number', 'expver') if v in ds.coords or v in ds.variables]
            if drop: ds = ds.drop_vars(drop)
            if 'valid_time' in ds.dims and 'time' not in ds.dims: ds = ds.rename({'valid_time': 'time'})
            if key == 'uv':
                for v in ds.data_vars:
                    if 'uv' in v.lower() or 'uv' in ds[v].attrs.get('long_name', '').lower():
                        ds[v].attrs['standard_name'] = 'surface_downward_uv_radiation'
            readers.append(rg.Reader(ds, name=key))
        else:
            readers.append(rg.Reader(p, name=key))
    if 'cmems_cur_glorys_abkatun_2026_jan20_mar25.nc' in missing:
        sys.exit(f'FAIL: ocean currents file not found in {FD}')
    if missing:
        print('WARNING: missing forcing files (fallbacks used):', missing)
    return readers


# ------------------------------------------------------------------ TAMOC -> CICOIL, the v1 route
def tamoc_v1_route(o, readers, t0):
    """Steps 1-5 of the v1 route (see the docstring).  Returns the number of elements scheduled."""
    import numpy as np
    from datetime import timedelta
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        from opendrift.models.openoil.vertical_profile import sample_profile, create_profile
    except ImportError:            # branch feat/vertical-profile not installed in the env yet: use the copy next to this script
        sys.path.insert(0, here)
        from vertical_profile import sample_profile, create_profile
    from opendrift.models.openoil.tamoc_plume import Plume
    import opendrift.models.openoil.run_tamoc as tamoc
    import adios_db.computation.physical_properties as props

    # RUN-TIME PATCH (CICOILv2 source untouched).  CICOILv2 b032140 handle_gas_particles() passes integer indices (np.where) to OpenDrift
    # 1.14.9 deactivate_elements(), which wants a boolean mask: with gas bubbles among the elements it stops in the first step with
    # "operands could not be broadcast together with shapes (4,) (19,)" (found with this example, Ixachi job 22817; reproduced locally).
    # The v2 end-to-end test T4 avoids it by setting processes:handle_released_gas False.  This is the boolean-mask version of the same function.
    def _handle_gas_particles_mask(self):
        bubbles = self.elements.is_bubble == 1
        if np.any(bubbles):
            self.deactivate_elements(bubbles, reason='released_gas')
    OpenCiceseOil.handle_gas_particles = _handle_gas_particles_mask
    print('run-time patch applied: OpenCiceseOil.handle_gas_particles uses a boolean mask (see the comment in abkatun2026.py)', flush=True)

    tdir = os.path.join(a.out, 'tamoc'); os.makedirs(tdir, exist_ok=True)
    name = os.path.join(tdir, f'tamoc_abkatun2026_{a.release_date}')
    plon = a.lon if a.profile_lon is None else a.profile_lon
    plat = a.lat if a.profile_lat is None else a.profile_lat

    # 1. live-oil composition (GOR sets how much gas the live oil carries)
    live_comp, composition, chemdata, chemunits = o.fluid_properties.get_live_composition(GOR=a.gor)

    # 2. ambient profile from the readers (T, S, u, v against depth) at the source
    zmax = abs(a.depth) + 20.0      # profile deep enough for TAMOC; the coupled CF files hold filled values below the real seafloor
    prof = sample_profile(readers, plon, plat, t0, z_max=zmax)
    print(f'TAMOC ambient profile at {plon:.3f}, {plat:.3f}: {len(prof["depth"])} levels, 0..{prof["depth"][-1]:.1f} m, '
          f'T {prof["temperature"][0] - 273.15:.1f}..{prof["temperature"][-1] - 273.15:.1f} C', flush=True)
    if prof['depth'][-1] < abs(a.depth):
        sys.exit(f'FAIL: the ambient profile ends at {prof["depth"][-1]:.1f} m (seafloor of the nearest grid cell) but the release is at '
                 f'{abs(a.depth):.1f} m; read the profile at a deeper ocean point with --profile-lon/--profile-lat')
    profile_nc = os.path.join(tdir, 'v_profile.nc')
    create_profile(profile_nc, readers, plon, plat, t0, z_max=zmax)

    # 3. bent plume model.  total_flow = mass flow of the live fluid, kg/s, from the oil rate at the ambient temperature of the source
    T_src = float(np.interp(abs(a.depth), prof['depth'], prof['temperature']))
    rho_oil = float(props.Density(o.fluid_properties.record).at_temp(T_src, 'K'))
    flow_kg_s = a.q_bbl_day * 0.158987 * rho_oil / 86400.0
    Tj = T_src + a.tj_delta
    print(f'TAMOC: {a.q_bbl_day:.0f} bbl/d x {rho_oil:.1f} kg/m3 = {flow_kg_s:.2f} kg/s, GOR {a.gor:.0f} scf/bbl, hole {a.hole} m, '
          f'z0 {abs(a.depth):.1f} m, Tj {Tj - 273.15:.1f} C, {a.bins} bins', flush=True)
    bpm = tamoc.run_tamoc(live_comp, composition, chemdata, chemunits, name, profile_nc, total_flow=flow_kg_s, z0=abs(a.depth),
                          D=a.hole, Tj=Tj, phi_0=-np.pi / 2., theta_0=0., bins=a.bins, tracking=False)

    # 4. plume + particle properties
    plume = Plume(plume_file=name + '_plume.nc')
    plume.get_particle_properties(particles_file=name + '_particles.nc', tracked=False)
    fp = np.asarray(plume.fp_type)
    print(f'TAMOC plume: {plume.nparticles} particle classes ({int((fp == 0).sum())} gas, {int((fp != 0).sum())} oil), '
          f'exit depths (zp) {np.min(plume.zp):.1f}..{np.max(plume.zp):.1f} m, plume width b {float(plume.b):.1f} m', flush=True)

    # 5. far-field elements
    o.seed_plume_elements(a.lon, a.lat, plume, [t0, t0 + timedelta(hours=a.release_hours)], number=a.particles,
                          z_uncertainty=a.z_unc, max_water=0.9)
    return plume


# ------------------------------------------------------------------ model
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.config import CONFIG_LEVEL_BASIC
rng = np.random.default_rng(a.seed)
np.random.seed(a.seed)

o = OpenCiceseOil(loglevel=a.loglevel, weathering_model='cicese', seed=a.seed)
o.required_variables = dict(OpenCiceseOil.required_variables)
o.required_variables['sea_surface_height'] = {'fallback': 0}
o._add_config({'environment:fallback:sea_surface_height': {'type': 'float', 'default': 0., 'min': -10., 'max': 10., 'units': 'm',
                                                          'description': 'fallback sea surface height', 'level': CONFIG_LEVEL_BASIC}})
readers = make_readers_coupled() if a.forcing == 'coupled' else make_readers()
o.add_reader(readers)
o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 30.0)
o.set_config('environment:fallback:sea_water_temperature', 25.0)
o.set_config('environment:fallback:sea_water_salinity', 36.0)
o.set_config('environment:fallback:ocean_mixed_layer_thickness', 20.0)
o.set_config('drift:vertical_mixing', a.mode == '3d')
o.set_config('drift:stokes_drift', True)
o.set_config('general:coastline_action', 'previous')
o.set_config('general:seafloor_action', 'previous')
for k, v in dict(evaporation=True, emulsification=True, dispersion=True, biodegradation=a.bio, photooxidation=a.photo).items():
    o.set_config(f'processes:{k}', v)
o.set_config('processes:dispersion_mode', a.dispersion_mode)
o.set_oiltype('MAYA')

t0 = datetime.strptime(a.release_date, '%Y-%m-%d')
use_plume = (a.mode == '3d' and a.droplets == 'tamoc')
if use_plume:
    tamoc_v1_route(o, readers, t0)                     # seeds the elements itself (seed_plume_elements)
else:
    kw = dict(lon=a.lon, lat=a.lat, number=a.particles, time=t0, m3_per_hour=a.volume_m3, oil_type='MAYA', max_water=0.9)
    if a.mode == '2d':
        kw['z'] = 0.0
    else:
        kw['z'] = -abs(a.depth)
        if a.droplets.startswith('lognormal:'):
            _, d50, sg = a.droplets.split(':')
            kw['diameter'] = np.clip(rng.lognormal(np.log(float(d50) * 1e-3), float(sg), a.particles), 1e-5, 2e-2)
        elif a.droplets != 'default':
            sys.exit('--droplets must be tamoc, default or lognormal:D50_mm:SIGMA')
        if 'diameter' in kw:
            print(f'droplet diameters (mm): min {kw["diameter"].min()*1e3:.2f}  median {np.median(kw["diameter"])*1e3:.2f}  max {kw["diameter"].max()*1e3:.2f}', flush=True)
    o.seed_elements(**kw)

_m = np.atleast_1d(np.asarray(o.elements_scheduled.mass_oil, dtype=float))
seeded = float(np.nansum(_m)) if len(_m) > 1 else float(_m[0]) * o.num_elements_total()   # oil mass of ALL scheduled elements; a scalar means equal mass per element
_g = np.atleast_1d(np.asarray(o.elements_scheduled.mass_gas, dtype=float))
gas_seeded = float(np.nansum(_g)) if len(_g) > 1 else float(_g[0]) * o.num_elements_total()   # bubble elements: released to the atmosphere
tag = f'abkatun2026_{a.mode}_{a.dispersion_mode}_s{a.seed}'
nc = os.path.join(a.out, tag + '.nc')
o.run(duration=timedelta(days=a.days), time_step=timedelta(minutes=a.dt), time_step_output=timedelta(hours=a.out_dt), outfile=nc)

b = o.get_oil_budget()
tot = seeded                              # oil mass only: bubble elements carry their mass in mass_gas and have mass_oil = 0 (checked with seed_elements),
                                          # so `seeded` (sum of mass_oil) already excludes the gas. Subtracting gas_seeded again gave closure 102.50 % in job 22830.
                                          # NOT mass_total[0]: with a continuous release that is only what was seeded at the first step
print(f'\nDONE {nc}')
n_el = o.num_elements_total()
print(f'mode {a.mode} | dispersion_mode {a.dispersion_mode} | droplets {a.droplets if a.mode == "3d" else "-"} | requested N {a.particles} | seeded oil {tot:.0f} kg + gas bubbles {gas_seeded:.0f} kg (released to the atmosphere, not in the oil budget; {100 * gas_seeded / max(tot + gas_seeded, 1):.1f} % of the total)')
print(f'elements in the run: {n_el} -> mean mass per element {tot / max(n_el, 1):.4f} kg (volume x density / N; not a constant)')
print('budget at the end of the run (% of the seeded OIL mass):')
s = 0.0
for k in ('mass_surface', 'mass_submerged', 'mass_evaporated', 'mass_dispersed', 'mass_biodegraded', 'mass_photooxidized', 'mass_stranded', 'mass_dissolved_subsea', 'mass_gas'):
    v = float(np.nan_to_num(b[k][-1])); s += v
    print(f'  {k:24s} {100 * v / tot:7.2f} %')
print(f'  closure (sum of the above) {100 * s / tot:7.2f} %  (should be 100.00; budget mass_total at the end = {100 * float(b["mass_total"][-1]) / tot:.2f} %)')
