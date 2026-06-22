# CICOILv2.0 Reference Manual

**Complete functionality reference for `cicoilv2_deployment/`**

This is the companion reference to [`CICOIL_v2_MANUAL.md`](CICOIL_v2_MANUAL.md)
(overview, deployment, photooxidation/biodegradation calibration). Where the
Manual explains *how to deploy and run* CICOIL v2 and the physics of the two
CICESE weathering additions, this Reference catalogs *every module, class,
method, function, element variable and config key* available in the package,
for use as an API lookup.

**Source:** `OILSPILL/CICOIL_dev/cicoilv2_deployment/`

---

## 1. Package layout

```
cicoilv2_deployment/
├── ciceseoil.py              # CiceseOil (element type) + OpenCiceseOil (model class)
├── fluid_properties.py       # FluidProps: 17-component distillation framework
├── cicoil_estimations.py      # Standalone ADIOS-DB estimation functions
├── tamoc_plume.py             # Plume class — TAMOC BPM output → OpenDrift seeding
├── run_tamoc.py                # run_tamoc() — TAMOC BPM wrapper
├── tamoc_chemical_properties.py
├── readers/
│   ├── reader_NEMO_native_v3.py    # active — CROCO native variable names
│   ├── reader_nemo_optimized.py    # active — generic NEMO-CMEMS names
│   ├── reader_nemo_combined.py     # alternative — xarray-based, see §8.3
│   └── reader_nemo_modified.py     # alternative — netCDF4-based, see §8.3
├── export/
│   └── io_stat_nc.py           # gridded statistical NetCDF export
├── plotting.py                 # publication plotting helpers (8 functions)
├── download_data.py            # GLORYS12/ANFC/ERA5 download + preprocessing
├── preprocess_era5_uv.py        # ERA5 UV → CF-compliant NetCDF (photooxidation forcing)
├── preprocess_real_forcing.py
├── cicoil_common.py             # shared helpers (load_ciceseoil, default_cicese_config)
├── test_cicoil_e2e.py           # 30-test full suite (T1–T8)
├── test_biodegradation.py        # Q10 sanity test
├── test_biodegradation_water.py  # dissolved-phase pathway test
├── CICOIL_v2_MANUAL.md           # deployment + weathering-module physics
└── CICOIL_v2_REFERENCE.md         # this file
```

---

## 2. Element variables (`CiceseOil`)

`CiceseOil` extends OpenDrift's `Oil` element type. Total: **97 variables**
(22 inherited from upstream `Oil`/`OpenOil` + 75 CICOIL-specific below).

### 2.1 Core CICOIL state variables

| Variable | Units | Default | Description |
|---|---|---|---|
| `spillet_area` | m² | 0.114 | 1-side surface area of the spillet |
| `oil_molar_mass` | kg/mol | 0.2 | |
| `spillet_thickness` | m | 0.01 | |
| `oil_density` | kg/m³ | 880 | |
| `DOR` | — | 0 | Dispersant-to-Oil Ratio |
| `IFT` | N/m | 30.0 | Oil–water interfacial tension |
| `max_water` | — | 0 | Maximum water content (emulsification cap) |
| `mass_dissolved_surface` | kg | 0 | |
| `mass_dissolved_subsea` | kg | 0 | Used as the dissolved-pool reference total (§7.3 of Manual) |
| `is_bubble` | — | 0 | Flag: gas-bubble particle |
| `is_dissolved` | — | 0 | Flag: dissolved subsea particle |
| `mass_gas` | kg | 0 | |
| `mass_biodegraded_from_oil` | kg | 0 | NEW — written by `biodegradation_cicese()` |
| `mass_biodegraded_from_water` | kg | 0 | NEW — written by `biodegradation_cicese()` (dissolved pathway) |
| `fraction_dissolved_surface` | % | 0 | |
| `fraction_biodegraded_from_oil` | % | 0 | NEW |
| `fraction_biodegraded_from_water` | % | 0 | NEW |
| `mass_photooxidized` | kg | 0 | NEW — written by `photooxidation_cicese()` |
| `fraction_photooxidized` | % | 0 | NEW |
| `mass_op_dissolved` | kg | 0 | NEW (audit #4) — dissolved oxygenated photoproduct pool, fed by `photooxidation_cicese()` |
| `fraction_op_dissolved` | % | 0 | NEW (audit #4) |
| `mass_op_degraded` | kg | 0 | NEW (audit #4) — decay product of `mass_op_dissolved` via `K_OP_REMOVAL` |
| `fraction_op_degraded` | % | 0 | NEW (audit #4) |

### 2.2 Pseudo-component arrays (17 components × 3 = 51 variables)

For each of the 17 pseudo-components below, three per-particle variables
exist: `<comp>`, `<comp>_evaporated`, `<comp>_dissolved` (all kg, default 0).

```
comp_ali_G1 .. comp_ali_G8     (8 aliphatic groups)
comp_arom_G1 .. comp_arom_G8   (8 aromatic groups)
comp_residue                   (1 non-volatile residue)
```

G1 = lightest (e.g. benzene/light alkanes) → G8 = heaviest
(pyrene+/heavy waxes). These are read/written by `evaporation_cicese()`,
`biodegradation_cicese()`, and `photooxidation_cicese()` (aromatics only),
and tracked in `self.cicese_mass_balance['mass_components']`.

### 2.3 Inherited variables (from upstream `OpenOil`/`Oil`, selected)
`mass_oil`, `mass_evaporated`, `mass_dispersed`, `mass_biodegraded`,
`density`, `viscosity`, `water_fraction`, `age_seconded`, `oil_film_thickness`,
`diameter`, `wind_drift_factor`, `bulltime`, plus position/status/ID — see
OpenDrift `OpenOil` documentation for the full 22-variable inherited set.

---

## 3. Configuration keys

### 3.1 CICOIL/CICESE-specific & OD 1.14.9-added keys

| Key | Type | Default | Description |
|---|---|---|---|
| `seed:m3_per_hour` | float | 1 | Volume of oil released per hour (or total if instantaneous) |
| `seed:droplet_diameter_min_subsea` | float | 0.0005 m | Min droplet diameter, subsea release |
| `seed:droplet_diameter_max_subsea` | float | 0.005 m | Max droplet diameter, subsea release |
| `seed:droplet_size_distribution` | enum | `'uniform'` | `uniform`/`normal`/`lognormal`, subsea release |
| `seed:droplet_diameter_mu` | float | 0.001 m | Mean diameter for normal/lognormal |
| `seed:droplet_diameter_sigma` | float | 0.0005 m | Std dev for normal/lognormal |
| `seed:oil_type` | enum | first in `self.oiltypes` | ADIOS oil name |
| `processes:dispersion` | bool | False | Natural dispersion (entrainment as small droplets) |
| `processes:chemical_dispersion` | enum | No\_treatment | No\_treatment / Surface / SSDI |
| `processes:evaporation` | bool | True | Cox Pv evaporation |
| `processes:emulsification` | bool | True | Bullwinkle-based water uptake |
| `processes:biodegradation` | bool | **False** | NEW — `biodegradation_cicese()`, see Manual §7 |
| `processes:update_oilfilm_thickness` | bool | False | Recompute film thickness each step vs. keep constant |
| `processes:handle_released_gas` | bool | True | Gas-bubble bookkeeping (`handle_gas_particles()`) |
| `processes:subsea_dissolution` | bool | True | Feeds the dissolved-mass pool (`handle_subsea_dissolution()`) |
| `processes:spreading` | bool | True | Nihoul (1984) spreading |
| `processes:surface_dissolution` | bool | False | |
| `processes:photooxidation` | bool | **False** | NEW — `photooxidation_cicese()`, see Manual §6 |
| `wave_entrainment:` `droplet_size_distribution` | enum | Johansen et al. (2015) | Johansen et al. (2015) / Li et al. (2017) |
| `wave_entrainment:entrainment_rate` | enum | Li et al. (2017) | |
| `biodegradation:method` | enum | `'Adcroft'` | Forced to `'Adcroft'`; describes the Q10 law in `biodegradation_cicese()` (not the upstream bulk formula) |

### 3.2 Drift defaults set by `OpenCiceseOil.__init__`

| Key | Value | Notes |
|---|---|---|
| `drift:vertical_advection` | False | |
| `drift:vertical_mixing` | True | **Required** for correct photooxidation surface-fraction (Manual §6.2) |
| `drift:current_uncertainty` | 0.05 | |
| `drift:wind_uncertainty` | 0.5 | |
| `drift:profiles_depth` | 20 | OD 1.14.9 plural key; replaces obsolete `required_profiles_z_range` |

---

## 4. `OpenCiceseOil` — method reference

`class OpenCiceseOil(OpenOil)` in `ciceseoil.py`. Element type:
`ElementType = CiceseOil`. `required_variables` includes the standard ocean
current/wind/wave/ice/temperature/salinity fields plus
`surface_downward_uv_radiation` (fallback `-1.0`, used by
`_compute_uv_irradiance()`).

### 4.1 Construction & oil selection

| Method | Signature | Description |
|---|---|---|
| `__init__` | `(self, weathering_model='noaa', *args, **kwargs)` | `weathering_model` must be `'noaa'` or `'cicese'`. Builds `self.oiltypes` from the ADIOS database (GENERIC oils first, `duplicate_oils` filtered out), registers all config keys (§3), sets drift defaults. |
| `set_oiltype` | `(self, oiltype)` | Sets `self.oiltype` (an `adios_db` `Oil`/`CiceseOil` fluid object) and, if `weathering_model='cicese'`, builds `self.fluid_properties = FluidProps(...)` (§5) for the 17-component framework. |
| `get_max_water_content` | `(self)` | Returns the oil's maximum water (emulsification) fraction, applying `_cicoil_max_water_fraction_override` for the two oils with known ADIOS data gaps (`MARINE GAS OIL 500 ppm S 2017`, `FENJA (PIL) 2015`). |
| `run` | `(self, record_components=False, record_environment=True, *args, **kwargs)` | Thin wrapper over `OpenOil.run()`. By default, `export_variables` is set to exclude the 51 per-component (`comp_*`, `comp_*_evaporated`, `comp_*_dissolved`) arrays, which otherwise dominate output size for large ensembles. `record_components=True` retains them. `record_environment=False` additionally drops the 20 `required_variables` environment fields (current, wind, waves, SST, salinity, etc.) sampled per element -- useful for further shrinking large-ensemble output, since these can be recovered from the forcing data. Both flags are ignored if an explicit `export_variables` kwarg is passed. NetCDF output is zlib-compressed (complevel 6) by OpenDrift core on file close, independent of these flags. |
| `prepare_run` | `(self)` | Pre-run setup: initializes `cicese_mass_balance` per-component tracking arrays, populates `self.max_water_fraction` from ADIOS JSON (NOAA mode). |

### 4.2 Seeding

| Method | Signature | Description |
|---|---|---|
| `seed_elements` | `(self, *args, **kwargs)` | Main seeding entry point. Accepts `lon`, `lat`, `time`, `number`, `z`, `m3_per_hour`, `oil_type` (preferred) or deprecated `oiltype`, `max_water`, `is_bubble`, `is_dissolved`, `diameter`. Computes `mass_oil` per element from `m3_per_hour × duration / number × density`; for subsea `z<0` without explicit `diameter`, draws droplet diameters uniformly between `seed:droplet_diameter_min/max_subsea`. Sets `self.Density`/`self.KinematicViscosity` (NOAA mode) or `density`/`oil_density`/`viscosity` (cicese mode) from `oiltype.density_at_temp(285)` / `kvis_at_temp(285)`. `is_bubble=1` → mass routed to `mass_gas`; `is_dissolved=1` → mass routed to `mass_dissolved_subsea`. |
| `seed_with_dispersant` | `(self, lon, lat, seed_time, appl_method, DOR, dispersant='C9500', time_frac_surf, time_frac_ssdi, effic_surf, effic_ssdi, wind_limit, delay, stock_kg, time_step, appl_rate, number)` | Seeds elements representing a chemically-dispersed slick: computes treated IFT/max-water via `treated_IFT_and_WC()`, samples effective DOR via `dispersant_efficiency_sampled()`, and checks weather windows via `weather_permitting_time()`. |
| `seed_plume_elements` | `(self, lon, lat, plume, seed_time, *args, **kwargs)` | Seeds particles from a TAMOC `Plume` object (§7.1) — i.e. a subsea blowout. Distributes droplets/gas-bubbles/dissolved-phase elements according to the plume's particle-tracking output (`Plume.get_particle_properties()`), sets `is_bubble`/`is_dissolved` flags, and finalizes the environment (with a 28°C temperature fallback if not yet available) before transitioning the model from `Mode.Config` to `Mode.Ready` via `__set_seed_config__()`. |
| `weather_permitting_time` | `(self, lon, lat, treatment_time, wind_limit, time_step)` | Computes the cumulative time window during which wind speed is below `wind_limit` at the given position/time — used to gate dispersant application. |

### 4.3 Weathering processes (dispatch & physics)

| Method | Signature | Description |
|---|---|---|
| `oil_weathering` | `(self)` | Top-level dispatcher; calls `oil_weathering_cicese()` (or the upstream NOAA path) depending on `oil_weathering_model`. |
| `oil_weathering_cicese` | `(self)` | Per-timestep weathering pipeline (CICOIL scheme adapted from NOAA PyGNOME). Calls, in order: `evaporation_cicese()`, `emulsification_cicese()`, `spreading_cicese()`, `dispersion_cicese()`, `handle_subsea_dissolution()`, `handle_gas_particles()`, then conditionally `biodegradation_cicese()` (if `processes:biodegradation`) and `photooxidation_cicese()` (if `processes:photooxidation`). |
| `evaporation_cicese` | `(self, atm_pressure=101325.)` | Cox vapor-pressure evaporation per pseudo-component; uses `diffusion_coeff_in_air()` and `transfer_coeff_air()`. Moves mass from `comp_*` → `comp_*_evaporated` and `mass_oil` → `mass_evaporated`. |
| `spreading_cicese` | `(self, gravity=9.81, drag_coef=0.003, film_limit=0.000005)` | Gravity-viscous spreading (Nihoul, 1984); updates `spillet_area`/`spillet_thickness`. `film_limit` is the minimum film thickness (m). |
| `emulsification_cicese` | `(self)` | Bullwinkle-type water-in-oil emulsification; updates `max_water`, `water_fraction`, `viscosity`, `IFT`. Requires `max_water>0` at seed time for oils lacking ADIOS emulsification data (e.g. MAYA). The effective water-uptake rate is: `k_emul_eff = k_emul * clip(nu_dry/EMUL_STABILITY_NU_REF, EMUL_STABILITY_MIN, EMUL_STABILITY_MAX) * EMUL_RATE_SCALE`, where `nu_dry` [cSt] is the dry-oil (evaporation-corrected) kinematic viscosity computed each step. Class-level constants: `EMUL_STABILITY_NU_REF=1000.0`, `EMUL_STABILITY_MIN=0.3`, `EMUL_STABILITY_MAX=3.0`, `EMUL_RATE_SCALE=0.0142` (calibrated to a 17.2-day tar-ball transition for MAYA crude at 10 m/s wind; achieved 16.92 d). Set all three to 1.0 to recover the stock NOAA Bullwinkle rate. See Manual §8. |
| `dispersion_cicese` | `(self)` | Natural dispersion via wave breaking; calls `oil_wave_entrainment_rate()` and one of the droplet-size-distribution methods (§4.4) to move mass to `mass_dispersed` / subsea droplets. |
| `handle_subsea_dissolution` | `(self)` | Dissolved-phase mass bookkeeping for `is_dissolved==1` particles. On first encounter, transfers `cicese_mass_balance['mass_components']` → `mass_dissolved` (one-time, zeroing components) and sets `mass_dissolved_subsea` as the reference total. See Manual §7.3 for the full lifecycle (this is the entry point that `biodegradation_cicese()`'s dissolved branch consumes). |
| `handle_gas_particles` | `(self)` | For `is_bubble==1` particles: records `mass_gas`, then deactivates the particle (status `released_gas`) — gas bubble mass is tracked but not advected further. |
| `biodegradation_cicese` | `(self)` | **NEW.** Per-component, Q10-scaled (T_REF=20°C, `BIODEG_Q10=3`) first-order decay using `K_BIODEG_ALI`/`K_BIODEG_AROM`/`K_BIODEG_RESIDUE`, with `BIODEG_SURFACE_FACTOR=1/30` for `z==0` particles and a droplet-size scaling factor `(BIODEG_DROPLET_D_REF/d)^BIODEG_DROPLET_EXPONENT` (audit #2; `d=elements.diameter` clamped to `[BIODEG_DROPLET_DIAMETER_MIN, BIODEG_DROPLET_DIAMETER_MAX]`, `diameter==0`→factor=1), plus a dissolved-phase branch (`K_BIODEG_WATER=0.0578 d⁻¹`) for `is_dissolved==1` particles, deactivating once `mass_dissolved < BIODEG_WATER_DEPLETION_FRAC` (0.01) of its initial reference. Writes `mass_biodegraded_from_oil/_water`, `mass_biodegraded`, `fraction_biodegraded_from_oil/_water`. Full physics in Manual §7. |
| `photooxidation_cicese` | `(self)` | **NEW.** UV-driven first-order decay of aromatic groups G1-G8 only, using `K_PHOTO_AROM`, Beer-Lambert depth attenuation (`K_D_UV=0.10 m⁻¹`, cutoff `UV_DEPTH_CUTOFF=50 m`), and `_compute_uv_irradiance()`. Of each step's photooxidized mass, a fraction `OP_DISSOLVED_FRACTION=0.5` routes to the dissolved oxygenated-photoproduct pool `mass_op_dissolved` instead of `mass_photooxidized` (audit #4); `mass_op_dissolved` decays unconditionally (incl. at night) at `K_OP_REMOVAL=0.116 d⁻¹` (Q10-scaled) into `mass_op_degraded`. Writes `mass_photooxidized`, `mass_op_dissolved`, `mass_op_degraded`, `fraction_photooxidized`, `fraction_op_dissolved`, `fraction_op_degraded`. Full physics in Manual §6. |

### 4.4 Wave entrainment & droplet size distributions

| Method | Signature | Description |
|---|---|---|
| `oil_wave_entrainment_rate` | `(self)` | Entrainment rate of oil into the water column due to wave breaking (`wave_entrainment:entrainment_rate` config, currently always Li et al. 2017). |
| `get_wave_breaking_droplet_diameter_johansen2015` | `(self)` | Droplet-size spectrum after wave-breaking entrainment, Johansen et al. (2015) parameterization. |
| `get_wave_breaking_droplet_diameter_liz2017` | `(self)` | Droplet-size spectrum, Li et al. (2017) parameterization. Selected via `wave_entrainment:droplet_size_distribution`. |

### 4.5 Dispersants

| Method | Signature | Description |
|---|---|---|
| `dispersant_efficiency_sampled` | `(self, target_DOR, efficiency, samples)` | Monte-Carlo sampling of effective dispersant-to-oil ratio given a target DOR and an efficiency distribution. |
| `treated_IFT_and_WC` | `(self, dispersant, applied_DORs)` | Computes the post-treatment interfacial tension and max-water-content for a given `dispersant` (one of `OpenCiceseOil.dispersants`: `C9500`, `OSR-52`, `Dasic-NS`, `Average`) and array of applied DORs. |

### 4.6 Physical property helpers

| Method | Signature | Description |
|---|---|---|
| `oil_density` | `(self, indices)` | Per-particle oil density at given element indices, accounting for current composition/temperature. |
| `component_molfrac` | `(self, surface)` | Mole fractions of the 17 pseudo-components for `surface` (boolean mask) elements — used by evaporation. |
| `diffusion_coeff_in_air` | `(self, atm_pressure=101325.)` | Vapour diffusion coefficients in air (Hilal & Karickhoff 2003; Kim & Monroe 2014), per pseudo-component, function of `atm_pressure`. |
| `transfer_coeff_air` | `(self, surface_indices)` | Mass-transfer coefficients at the air-sea interface for `surface_indices` elements (wind-speed dependent), used by `evaporation_cicese()`. |
| `_compute_uv_irradiance` | `(self)` | Surface (z=0) solar UV irradiance (W/m²). Prefers real `surface_downward_uv_radiation` reader data (ERA5), falls back to a synthetic diurnal cosine model. Returns 0 at night. See Manual §6.2. |

### 4.7 Output, budgets & analysis

| Method | Signature | Description |
|---|---|---|
| `get_oil_budget` | `(self)` | Returns a dict with `oil_density`, `mass_dispersed`, `mass_submerged`, `mass_surface`, `mass_stranded`, `mass_dissolved_subsea`, `mass_gas`, `mass_evaporated`, `mass_biodegraded`, `mass_photooxidized`, `mass_total` — each a 1-D time series array (summed over trajectories). `mass_total` should equal the released mass (mass conservation check). Returns `None` for backward simulations. Uses the `_prop()` helper for OD 1.14.9 xarray compatibility. |
| `plot_oil_budget` | `(self, filename=None, ax=None, show_density_viscosity=True, show_wind_and_current=True)` | Stacked-area plot of `get_oil_budget()` categories over time, optionally with density/viscosity and wind/current sub-panels. Saves to `filename` or shows interactively. |
| `plot_density` | `(self, value_array, lon_array, lat_array, title=None, label=None, min_value=None, max_value=None, cmapname=...)` | Plots a gridded density field (e.g. from `get_density_array_framed()`) on a map. |
| `plot_stats_property` | `(self, value_array, lon_array, lat_array, title=None, label=None, min_value=None, max_value=None, cmapname=...)` | Plots a gridded statistical property (e.g. from `get_property_binned()`) on a map. |
| `create_new_map` | `(self, corners=None, buffer=0.1, lscale='auto', fast=False)` | Creates a Cartopy `Figure`/`GeoAxes` for trajectory/density plots, with land features at resolution `lscale`. |
| `get_bins` | `(self, bin_size_deg, lat_min, lat_max, lon_min, lon_max, depths=None)` | Builds a regular lon/lat (and optionally depth) bin grid of size `bin_size_deg` for spatial aggregation. |
| `get_property_binned` | `(self, bins, p_groups, depths, de_sizes, property='particles_no')` | Aggregates a named element `property` (default: particle counts) into the `bins` grid, optionally split by pseudo-component group (`p_groups`), depth bin, and droplet-size bin. |
| `get_unique_visits` | `(self, bins, p_groups, depths, de_sizes)` | Counts unique particle visits per spatial bin (avoids double-counting a particle that revisits a cell). |
| `get_density_array_framed` | `(self, bin_size_deg, lat_min, lat_max, lon_min, lon_max, weight=None)` | Convenience wrapper combining `get_bins()` + `get_property_binned()` to produce a 2-D mass/particle-density array over a fixed lon/lat frame, optionally weighted by `weight` (e.g. `mass_oil`). |

---

## 5. `fluid_properties.FluidProps`

`class FluidProps(object)` in `fluid_properties.py` — the 17-component
pseudo-component ("distillation cut") framework underlying CICOIL's
component-resolved weathering (evaporation, biodegradation, photooxidation).

| Member | Signature | Description |
|---|---|---|
| `_riazi` | `(T, A, B)` (module function) | Riazi correlation helper: generic `Y = exp(A + B*T)`-type fit used by several property estimators below. |
| `__init__` | `(self, oiltype, max_cuts)` | Builds the component framework for `oiltype` (an ADIOS `Oil`/`CiceseOil` object) with up to `max_cuts` distillation cuts (CICOIL uses 17 = 8 aliphatic + 8 aromatic + 1 residue). |
| `_create_TBP` | `(self, max_cuts)` | Builds the True-Boiling-Point (TBP) distillation curve and per-cut boiling points, used to derive molecular weights, densities, and vapor pressures for each pseudo-component. |
| `get_Riazi_model` | `(self)` | Returns the Riazi-correlation-based property set (mol. weight, density, critical properties) for each pseudo-component. |
| `get_live_composition` | `(self, GOR=0)` | Returns the "live" (in-situ, pressure-adjusted) composition accounting for Gas-to-Oil Ratio `GOR` — relevant for subsea/TAMOC releases where dissolved gas affects bulk properties. |
| `setup_treatment` | `(self, treatment, dispersant, DOR, disp_stock, delay, time_frac_surf, time_frac_ssdi, effic_surf, effic_ssdi, wind_limit, equipment_rate)` | Configures a chemical-dispersant treatment scenario (used by `seed_with_dispersant()`/`OpenCiceseOil`); stores treatment parameters for later IFT/max-water adjustments. |
| `verify_cut_fractional_masses` | `(cls, fmass_i, T_i, f_sat_i, f_arom_i, prev_f_sat_i)` (classmethod) | Sanity-checks that per-cut fractional masses (`fmass_i`), boiling points (`T_i`), and saturate/aromatic splits (`f_sat_i`, `f_arom_i`) are physically consistent (monotonic, sum to 1, etc.) given the previous cut's saturate fraction. |

---

## 6. `cicoil_estimations.py` — ADIOS-gap estimation toolbox

Standalone replacement for `adios.util.estimations` (not present in
`adios_db`). 32 functions, all pure/stateless, mostly empirical correlations
from ADIOS2 documentation, Riazi, Fingas, or Bill Lehr's recommendations.
Grouped by purpose:

### 6.1 Density / API gravity / specific gravity
| Function | Signature |
|---|---|
| `density_from_api` | `(api)` |
| `api_from_density` | `(density)` |
| `density_at_temp` | `(ref_density, ref_temp_k, temp_k, k_rho_t=0.0008)` |
| `vol_expansion_coeff` | `(rho_0, t_0, rho_1, t_1)` |
| `specific_gravity` | `(density)` |

### 6.2 Viscosity
| Function | Signature |
|---|---|
| `dvis_to_kvis` | `(dvis, density)` — dynamic → kinematic viscosity |
| `kvis_at_temp` | `(ref_kvis, ref_temp_k, temp_k, k_v2=2416.0)` |

### 6.3 SARA fractions (Saturates/Aromatics/Resins/Asphaltenes)
| Function | Signature |
|---|---|
| `resin_fraction` | `(density, viscosity, f_other=0.0)` |
| `asphaltene_fraction` | `(density, viscosity, f_other=0.0)` |
| `saturates_fraction` | `(density, viscosity, f_other=0.0)` |
| `aromatics_fraction` | `(f_res, f_asph, f_sat)` |
| `_A_coeff` | `(density)` — Fingas empirical coefficient |
| `_B_coeff` | `(density, viscosity)` — Fingas empirical coefficient |

### 6.4 Distillation cuts
| Function | Signature |
|---|---|
| `cut_temps_from_api` | `(api, N=5)` — boiling points for N cuts from API gravity |
| `fmasses_from_cuts` | `(f_evap_i)` — fractional masses from evaporation-curve cuts |
| `fmasses_flat_dist` | `(f_res, f_asph, N=5)` — flat-distribution fallback for N cuts |

### 6.5 Per-cut molecular weight & density (Riazi correlations)
| Function | Signature |
|---|---|
| `saturate_mol_wt` | `(boiling_point)` |
| `aromatic_mol_wt` | `(boiling_point)` |
| `resin_mol_wt` | `(_boiling_point)` |
| `asphaltene_mol_wt` | `(_boiling_point)` |
| `trial_densities` | `(boiling_points, watson_factor)` |
| `saturate_densities` | `(boiling_points)` |
| `aromatic_densities` | `(boiling_points)` |
| `resin_densities` | `(_boiling_points)` |
| `asphaltene_densities` | `(_boiling_points)` |

### 6.6 Hydrocarbon characterization
| Function | Signature |
|---|---|
| `_hydrocarbon_characterization_param` | `(specific_gravity, temp_k)` |
| `refractive_index` | `(hc_char_param)` |
| `_hydrocarbon_grouping_param` | `(mol_wt, specific_gravity, temp_k)` |
| `saturate_mass_fraction` | `(fmass_i, mol_wt, specific_gravity, temp_k)` |

### 6.7 Interfacial tension, pour point, flash point, Bullwinkle
| Function | Signature |
|---|---|
| `oil_water_surface_tension_from_api` | `(api)` |
| `pour_point_from_kvis` | `(ref_kvis, ref_temp_k)` |
| `pour_point_from_sg_mw_kvis` | `(specific_gravity, mol_wt, kvis)` |
| `flash_point_from_bp` | `(temp_k)` |
| `flash_point_from_api` | `(api)` |
| `bullwinkle_fraction_from_asph` | `(f_asph)` — emulsification onset fraction from asphaltene content |
| `bullwinkle_fraction_from_api` | `(api)` — emulsification onset fraction from API gravity |

---

## 7. TAMOC coupling

### 7.1 `tamoc_plume.Plume`
`class Plume(object)` in `tamoc_plume.py` — represents the output of a TAMOC
Bent Plume Model (BPM) run, packaged for OpenDrift subsea seeding.

| Method | Signature | Description |
|---|---|---|
| `__init__` | `(self, chem_data_file=None, plume_file=None)` | Optionally loads chemical-data and plume-output files immediately. |
| `load_chem_data` | `(self, fname=None)` | Loads TAMOC chemical property tables (see `tamoc_chemical_properties.py`). |
| `load_plume` | `(self, fname=None)` | Loads the TAMOC BPM plume-trajectory output (centerline position, radius, velocity, temperature, phase fractions vs. height above release). |
| `get_particle_properties` | `(self, particles_file, tracked=False)` | Reads TAMOC's per-particle output (`particles_file`) and returns arrays of particle positions, sizes, phase (oil droplet / gas bubble / dissolved), and mass — consumed directly by `OpenCiceseOil.seed_plume_elements()`. `tracked=True` retains individual-particle trajectories rather than aggregating. |

### 7.2 `run_tamoc.run_tamoc`
```python
run_tamoc(live_comp, composition, chemdata, chemunits, output_name,
          profile, total_flow, z0, D, Tj, phi_0, theta_0, bins)
```
Wrapper that runs a TAMOC BPM simulation given:
- `live_comp` / `composition` / `chemdata` / `chemunits` — fluid composition and units (from `FluidProps.get_live_composition()` and `tamoc_chemical_properties.py`)
- `profile` — ambient ocean profile (T, S, currents vs. depth)
- `total_flow`, `z0`, `D`, `Tj`, `phi_0`, `theta_0` — release rate, depth, orifice diameter, jet temperature, and discharge angles
- `bins` — number of particle-size bins for the dispersed phase
- `output_name` — base filename for TAMOC output files, later read by `Plume.load_plume()` / `Plume.get_particle_properties()`

### 7.3 `tamoc_chemical_properties.py`
Static property tables (molecular weights, critical properties, Henry's law
constants, etc.) for the chemical components used by TAMOC's equilibrium and
dissolution calculations — feeds `run_tamoc()`.

---

## 8. NEMO/CROCO readers

Both active readers subclass `opendrift.readers.basereader.BaseReader` and
expose ocean fields via a `NEMO_variable_mapping` dict (native var → CF
`standard_name`), consumed automatically by OpenDrift's reader interface.

### 8.1 `reader_NEMO_native_v3.Reader` (active — CROCO native names)
```python
Reader(filename=None, filenameU=None, filenameV=None, name=None,
       gridfile=None, custom_var_mapping=None, _FillValue=None)
```
Variable mapping (CROCO/NEMO native diagnostic names):

| Native | CF standard_name |
|---|---|
| `sossheig` | `sea_surface_height` |
| `vozocrtx` | `x_sea_water_velocity` |
| `vomecrty` | `y_sea_water_velocity` |
| `vovecrtz` | `upward_sea_water_velocity` |
| `votemper` | `sea_water_temperature` |
| `vosaline` | `sea_water_salinity` |
| `utau` | `surface_downward_x_stress` |
| `vtau` | `surface_downward_y_stress` |

`custom_var_mapping` (dict) overrides/extends this mapping.

### 8.2 `reader_nemo_optimized.Reader` (active — generic NEMO/CMEMS names)
```python
Reader(filename=None, filenameU=None, filenameV=None, filenameW=None,
       filename_mask=None, name=None, gridfile=None,
       custom_var_mapping=None, _FillValue=None)
```
Variable mapping (CMEMS-style short names):

| Native | CF standard_name |
|---|---|
| `ssh` | `sea_surface_height` |
| `uoce` | `x_sea_water_velocity` |
| `voce` | `y_sea_water_velocity` |
| `wo` | `upward_sea_water_velocity` |
| `toce` | `sea_water_temperature` |
| `soce` | `sea_water_salinity` |
| `tauuo` | `surface_downward_x_stress` |
| `tauvo` | `surface_downward_y_stress` |
| `tmask` | `land_binary_mask` |

Supports separate U/V/W/mask files (`filenameU/V/W`, `filename_mask`) in
addition to a combined `filename`.

### 8.3 `reader_nemo_combined.py` / `reader_nemo_modified.py` — alternative NEMO readers
Both files were cleaned up (2026-06-14) by removing an earlier
fully-commented-out draft that preceded the working `class Reader` in each
file. The remaining live classes use the same
`ssh`/`uoce`/`voce`/`wo`/`toce`/`soce`/`tauuo`/`tauvo`/`tmask` mapping as
`reader_nemo_optimized.py`, with separate U/V/W/mask file support
(`filenameU/V/W`, `filename_mask`) and an optional `gridfile` for
`nav_lon`/`nav_lat`:

- `reader_nemo_combined.py` — xarray-based (`xr.open_dataset`, loads datasets
  into memory, interpolates via `scipy.ndimage.map_coordinates`). Imported
  and exercised by `test_cicoil_e2e.py` (T5.3).
- `reader_nemo_modified.py` — netCDF4-based (`Dataset`/`MFDataset`, supports
  glob patterns via `MFDataset`), with a `lonlat_to_ij()` helper for
  nearest-grid-point lookup. Not currently used by the test suite or
  production scripts; kept as an alternative netCDF4-based implementation.

`reader_nemo_optimized.py` and `reader_NEMO_native_v3.py` remain the primary
readers for production runs.

---

## 9. `export/io_stat_nc.py` — gridded statistical output

Module-level function (monkey-patched onto an OpenDrift model instance as an
export hook, NOAA/OpenDrift `io_*` convention):

| Function | Signature | Description |
|---|---|---|
| `init` | `(self, model_template, filename, prop_metadata)` | Opens a new NetCDF file (`filename`) for gridded statistical output, writing global attributes (`model_url`, `opendrift_class`, `opendrift_module`, `readers`, `time_step_calculation`) and dimension/variable definitions derived from `prop_metadata` (the gridded properties to be written each timestep, e.g. from `get_density_array_framed()` / `get_property_binned()`). |

Used for producing gridded (lon/lat/depth-binned) summary fields instead of
(or in addition to) per-trajectory NetCDF output — useful for large
ensembles where per-particle output would be too large.

---

## 10. `plotting.py`

8 publication-quality plotting helpers, all standalone functions (operate on
`get_oil_budget()` dicts, NetCDF output files, or an `OpenCiceseOil` instance
— no shared state).

| Function | Signature | Description |
|---|---|---|
| `plot_mass_budget` | `(budget, title=None, timestep_hours=1, unit='tonnes', show_photoox=True, filename=None, figsize=..., dpi=...)` | Stacked-area plot of an oil-budget dict (from `get_oil_budget()` or a saved `.npz`/`.csv`) over time, in `unit` (tonnes/kg/%). `show_photoox` toggles the photooxidized category. |
| `plot_trajectory_map` | `(nc_file, spill_lon=None, spill_lat=None, title=None, time_index=-1, color_by=None, filename=None, figsize=..., dpi=...)` | Cartopy map of particle trajectories from an OpenDrift output `nc_file`, optionally colored by `color_by` (e.g. `'status'`, `'z'`, a `comp_*` variable), marking the spill origin. |
| `plot_photoox_timeseries` | `(budget_with, budget_without, timestep_hours=1, title=None, filename=None, figsize=..., dpi=...)` | Compares `fraction_photooxidized` time series between a run with UV-driven photooxidation (`budget_with`) and a control run without (`budget_without`). |
| `plot_aromatic_budget` | `(o, title=None, filename=None, figsize=..., dpi=...)` | Stacked bar chart of per-aromatic-group (G1–G8) evaporated + photooxidized mass fractions, reading `comp_arom_G*` / `comp_arom_G*_evaporated` / photooxidation totals directly from an `OpenCiceseOil` instance `o`. |
| `plot_f_surface_sensitivity` | `(Hs_range, Tp_range, w_rise=..., n=..., title=None, filename=None, figsize=..., dpi=...)` | Heatmap of the (legacy, now-removed from the model) `f_surface` wave-submergence correction as a function of significant wave height `Hs_range` and peak period `Tp_range` — retained for diagnostic/historical comparison (see Manual §6.2 on why `f_surface` was removed from `photooxidation_cicese()`). |
| `plot_surface_density` | `(nc_file, time_index, spill_lon=None, spill_lat=None, bin_deg=..., weight_var='mass_oil', title=None, filename=None, figsize=..., dpi=...)` | Gridded surface oil-mass (or other `weight_var`) density map at timestep `time_index`, binned at `bin_deg` resolution. |
| `compare_weathering_modes` | `(budgets, labels, timestep_hours=1, title=None, filename=None, figsize=..., dpi=...)` | Side-by-side multi-panel comparison of oil-fate fractions across multiple `budgets` (e.g. with/without biodegradation, with/without photooxidation), each labeled per `labels`. |
| `plot_photoox_rate_constants` | `(filename=None, figsize=..., dpi=...)` | Bar chart of the current `K_PHOTO_AROM` G1–G8 rate constants / half-lives (reads directly from `ciceseoil.OpenCiceseOil.K_PHOTO_AROM`) — useful for documenting the current calibration (Manual §6.3). |

---

## 11. `download_data.py` — forcing acquisition & preprocessing

| Function | Signature | Description |
|---|---|---|
| `download_glorys12` | `(lon_min, lon_max, lat_min, lat_max, date_start, date_end, outfile, depth_min=..., depth_max=...)` | Downloads GLORYS12 reanalysis ocean currents/T/S (Copernicus Marine) for the given bounding box, date range, and depth range — recommended for reproducible hindcast runs (Manual §9 caveat re: ANFC). |
| `download_mercator_anfc` | `(lon_min, lon_max, lat_min, lat_max, date_start, date_end, outfile, depth_min=..., depth_max=..., tmpdir=...)` | Downloads Mercator ANFC near-real-time ocean forecast/analysis fields — use only when reproducibility across download dates is not required. |
| `download_era5_winds` | `(lon_min, lon_max, lat_min, lat_max, year, months, outfile, time_step=...)` | Downloads ERA5 10 m wind components (and related surface fields) via the CDS API. |
| `download_era5_waves` | `(lon_min, lon_max, lat_min, lat_max, year, months, outfile, time_step=...)` | Downloads ERA5 wave fields (Hs, Tp, Stokes drift) for wave-submergence and dispersion calculations. |
| `_preprocess_era5_cf` | `(ds)` | Internal: adds CF `standard_name` attributes and normalizes the time dimension/units on a raw ERA5 `xarray.Dataset` so OpenDrift's `reader_netCDF_CF_generic` recognizes the variables (required — raw ERA5 NetCDF has `standard_name='unknown'`). |
| `preprocess_era5_zip` | `(zip_nc_path, outfile)` | Unpacks a ZIP-disguised ERA5 NetCDF file (old CDS API delivery format) and applies `_preprocess_era5_cf`. |
| `check_domain_coverage` | `(forcing_files, spill_lon, spill_lat, sim_days, buffer_deg=...)` | Verifies that the given `forcing_files` spatially and temporally cover a simulation starting at `(spill_lon, spill_lat)` for `sim_days`, with a `buffer_deg` margin — run before launching a long simulation to catch missing-coverage failures early. |
| `list_croco_files` | `(croco_dir, pattern=...)` | Lists CROCO output NetCDF files in `croco_dir` matching `pattern`, for use with `reader_NEMO_native_v3.Reader`. |

> For ERA5 UV radiation specifically (photooxidation forcing), see
> `preprocess_era5_uv.py` (not part of `download_data.py`) — produces the
> CF-compliant `surface_downward_uv_radiation` files in `data_local/era5_uv_cf/`
> referenced in Manual §6.2.

---

**Version:** CICOILv2.0
**Last updated:** 2026-06-14
