# CICOILv2.0 Manual

**OpenDrift-based oil weathering model (CICOIL fork, OpenDrift 1.14.9 + CICESE extensions)**

This manual merges and supersedes `cicoilv2/hpc_migration.md` (deployment
workflow) and `DEPLOYMENT_MANIFEST.txt` (file inventory / API reference), and
adds the photooxidation, biodegradation, and emulsification-fix modules
developed through 2026-06-15.

**Active development copy:** `OILSPILL/CICOIL_dev/cicoilv2_deployment/`
**Contact:** adomingu@cicese.edu.mx — CICESE, Departamento de Oceanografía Física

> For a complete API/functionality catalog (every class, method, function,
> element variable, and config key — including FluidProps, TAMOC coupling,
> NEMO/CROCO readers, plotting and download utilities), see the companion
> [`CICOIL_v2_REFERENCE.md`](CICOIL_v2_REFERENCE.md).

---

## 1. Overview

CICOILv2.0 = OpenDrift 1.14.9 `OpenOil` + the CICOIL (Kotzakoulakis,
CICESE/SINTEF Ocean) weathering extensions, patched for compatibility with
OpenDrift 1.14.9, plus three CICESE-developed improvements:

- **Photooxidation** (`photooxidation_cicese()`) — UV-driven loss of
  aromatic pseudo-components at/near the surface.
- **Biodegradation** (`biodegradation_cicese()`) — component-resolved,
  temperature-scaled, with surface-slick suppression and a dissolved-phase
  pathway.
- **Emulsification viscosity-stability fix** (`emulsification_cicese()`) —
  viscosity-dependent scaling + empirical rate correction that prevents the
  stock NOAA Bullwinkle model from producing instant (~5 h) mousse formation
  and instead gives a physically calibrated ~17-day tar-ball transition for
  MAYA crude.

Photooxidation and biodegradation are additive (off by default) and conserve
mass exactly (rel. err ~1e-15–1e-16 in all test suites). The emulsification
fix is always active when `processes:emulsification = True`.

### What CICOIL adds over upstream OpenDrift 1.14.9 `OpenOil`
- 17-component weathering (Cox vapor-pressure evaporation, Nihoul spreading,
  per-particle density)
- TAMOC BPM coupling (`seed_plume_elements()`, `Plume`, `run_tamoc()`)
- Dispersant model (surface + SSDI, 4 dispersant types, per-particle
  IFT/max_water)
- 3 particle types: surface droplets, gas bubbles (`is_bubble`), dissolved
  subsea (`is_dissolved`)
- 93 element variables (vs 22 in vanilla OD 1.14.9 `OpenOil`)
- 4 NEMO/CROCO readers (native_v3, combined, modified, optimized)
- Statistical NetCDF output (`export/io_stat_nc.py`)
- **NEW (CICESE):** photooxidation, biodegradation, and emulsification fix
  described in §6–8

---

## 2. Installation & Deployment

### 2.1 Local development environment
- Source tree: `OILSPILL/CICOIL_dev/cicoilv2_deployment/`
- `OILSPILL/CICOIL` symlinks to the OpenDrift install location
  (`~/miniconda3/envs/numerics/lib/python3.11/site-packages/opendrift/`)
- Local conda env: `numerics` (Python 3.11, OpenDrift 1.14.9, TAMOC 4.x)

Run the test suite locally with:
```bash
conda run -n numerics python test_cicoil_e2e.py
conda run -n numerics python test_biodegradation.py
conda run -n numerics python test_biodegradation_water.py
```

### 2.2 HPC deployment (CHAMAN2)

**Target environment:** `cicoil4` (Python 3.10, cloned from `cicoil3`).
**Package manager:** conda (mamba is broken on CHAMAN2 — do not use
`deploy_on_hpc_mamba.sh`, use **`deploy_on_hpc_conda.sh`**).
**Safety:** `cicoil3` is never modified; `cicoil4` is a fresh, dedicated
environment that can be deleted and recreated.

#### Step 1 — Local staging
```bash
cd OILSPILL/CICOIL_dev/cicoilv2_deployment
bash stage_for_hpc.sh   # produces cicoilv2_hpc_deployment.tar.gz + DEPLOYMENT_MANIFEST.txt
```

#### Step 2 — Transfer
```bash
scp cicoilv2_hpc_deployment.tar.gz YOUR_USERNAME@chaman2.cicese.mx:~/cicoilv2_migration/
ssh YOUR_USERNAME@chaman2.cicese.mx
cd ~/cicoilv2_migration
tar -xzf cicoilv2_hpc_deployment.tar.gz
cd cicoilv2_deployment
```

#### Step 3 — Deploy
```bash
bash deploy_on_hpc_conda.sh
```
This script:
1. Creates `cicoil4` (cloned from `cicoil3`, conda)
2. Upgrades OpenDrift 1.9.0 → 1.14.9 in `cicoil4` only
3. Installs all patched/new files into `cicoil4` site-packages
4. Runs the full test suite (`test_cicoil_e2e.py`, `test_biodegradation.py`,
   `test_biodegradation_water.py`)
5. Writes `HPC_DEPLOYMENT_REPORT.txt`

Expected end-of-script summary:
```
=== HPC DEPLOYMENT REPORT ===
✓ OpenDrift version: 1.14.9
✓ Patched files installed
✓ test_cicoil_e2e.py: 30/30 PASSED
✓ test_biodegradation.py: PASSED (Q10 sanity check)
✓ test_biodegradation_water.py: PASSED
✓ Mass conservation: 100.000% (rel. err ~1e-15)
✓ Ready for production — activate with: conda activate cicoil4
```

#### Step 4 — Production run
```bash
cp slurm_template_cicoilv2.sh ~/my_simulation.slurm
nano ~/my_simulation.slurm   # set conda activate cicoil4, edit experiment details
sbatch ~/my_simulation.slurm
squeue -u YOUR_USERNAME
```

#### Step 5 — Verify output
```python
import xarray as xr
ds = xr.open_dataset('output.nc')
print(len(ds.time), ds.sizes['trajectory'])
print(set(ds.status.values.ravel()))
```

---

## 3. Quick Start

```python
import logging; logging.disable(logging.WARNING)
from opendrift.models.openoil.ciceseoil import OpenCiceseOil
from opendrift.readers.reader_netCDF_CF_generic import Reader as CFReader
from datetime import timedelta

o = OpenCiceseOil(weathering_model='cicese')
o.set_config('processes:evaporation', True)
o.set_config('processes:spreading', True)
o.set_config('processes:emulsification', True)
o.set_config('processes:photooxidation', True)
o.set_config('processes:biodegradation', True)

o.add_reader([CFReader('glorys12_ocean.nc'), CFReader('era5_wind_uv_cf.nc')])

o.seed_elements(lon=-92.1, lat=19.7, time='2023-07-06 12:00',
                number=3000, m3_per_hour=100.0,
                oil_type='MAYA', max_water=0.9)

o.run(duration=timedelta(days=30), time_step=timedelta(hours=1),
      time_step_output=timedelta(hours=6), outfile='output.nc')

o.plot_oil_budget(filename='oil_budget.png')
```

**Output size for large ensembles:** by default, `o.run()` drops the 51
per-pseudo-component arrays (`comp_*`, `comp_*_evaporated`, `comp_*_dissolved`)
from `output.nc`, since these dominate file size (17 components x 3 variants
x N particles x N output steps). Set `record_components=True` to keep them
(e.g. for debugging the mass balance). For very large ensembles, pass
`record_environment=False` to also drop the 20 per-element environment fields
(current, wind, waves, SST, salinity, etc.) -- these can be re-derived from
the forcing data and are often the largest remaining contributor to file size.
NetCDF output is zlib-compressed (complevel 6) automatically by OpenDrift on
file close, regardless of these flags.

**TAMOC subsea blowout (plume seeding):**
```python
from opendrift.models.openoil.tamoc_plume import Plume
from opendrift.models.openoil.run_tamoc import run_tamoc
# build Plume from run_tamoc(...) output, then:
o.seed_plume_elements(plume, ...)
```

**CROCO/NEMO forcing:**
```python
from opendrift.readers.reader_nemo_optimized import Reader as CrocoReader
o.add_reader(CrocoReader('croco_avg.nc'))
```

**Raw CROCO sigma-coordinate output (no preprocessing):**
```python
from opendrift.readers.reader_croco_native import Reader as CrocoNative      # installed with readers/*
from opendrift.readers.reader_netCDF_CF_generic import Reader as CF

ocean = CrocoNative('croco_avg.nc', name='croco')    # or croco_his.nc; native sigma, staggered C-grid
wind  = CF('wind_cf.nc', name='wind')                # CROCO files carry no atmosphere: supply wind/waves separately
o.add_reader([ocean, wind])
```
`reader_croco_native` subclasses OpenDrift's `reader_ROMS_native`, derives `mask_u`/`mask_v` from `mask_rho` and de-staggers `u`/`v` onto rho points (`destagger=True`, default; `False` = stock behaviour with a half-cell shift). It reads `croco_avg.nc` / `croco_his.nc` as they come out of CROCO; unlike `reader_NEMO_native_v3` it handles sigma coordinates and the C-grid. Ocean variables only; the first `croco_avg.nc` record is at the middle of the first averaging window, so a run cannot start earlier. Validated for surface transport (24 h backward / 48 h forward, 5000 particles: mean separation 0.2 m / 11 m against converted files, versus 2.1 km for the stock reader); subsurface fields are not yet validated (velocity differences of 0.01-0.04 m/s at 20 m depth against the reference). Tested with OpenDrift 1.14.9.

### Quick verification
```python
from opendrift.models.openoil.ciceseoil import OpenCiceseOil, CiceseOil
o = OpenCiceseOil(weathering_model='cicese')
o.set_oiltype('MAYA')
assert len(CiceseOil.variables) == 93
assert o.get_config('biodegradation:method') == 'Adcroft'
print("CICOIL OK:", len(o.oiltypes), "oils,", len(o.fluid_properties.composition), "components")
```

---

## 4. Configuration Reference (key new/changed config keys)

| Key | Default | Notes |
|---|---|---|
| `processes:evaporation` | True | Cox Pv evaporation |
| `processes:spreading` | True | Nihoul spreading |
| `processes:emulsification` | True | Bullwinkle model (viscosity-stability fix active; see §8); needs `max_water>0` at seed time (see §9 caveat) |
| `processes:dispersion` | True | Natural/wave dispersion |
| `processes:chemical_dispersion` | False | Surface + SSDI dispersant model |
| `processes:subsea_dissolution` | True | Feeds the dissolved-phase pool (§7.3) |
| `processes:surface_dissolution` | False | |
| `processes:update_oilfilm_thickness` | True | |
| `processes:handle_released_gas` | True | |
| `processes:photooxidation` | **False** | NEW — §6 |
| `processes:biodegradation` | **False** | NEW — §7 |
| `biodegradation:method` | `'Adcroft'` | Now describes the per-component Q10 law in `biodegradation_cicese()`, not the upstream bulk formula |
| `seed:oil_type` | — | ADIOS oil name, e.g. `'MAYA'` |
| `seed:m3_per_hour` | — | Required for continuous seeding |
| `seed:droplet_size_distribution` | — | `'normal'`/`'lognormal'`/`'uniform'` for subsea release |
| `seed:droplet_diameter_mu` / `_sigma` | — | Subsea droplet size distribution params |
| `seed:droplet_diameter_min/max_subsea` | — | |
| `drift:profiles_depth` | 20 | Plural key name (OD 1.14.9) |
| `drift:vertical_advection` | False | |
| `drift:vertical_mixing` | True | Required for correct photooxidation surface-fraction behavior (§6.2) |

---

## 5. API Compatibility Fixes Applied (14, vs upstream OD 1.14.9)

1. `adios.computation` → `adios_db.computation` (imports)
2. `adios.util.estimations` → standalone `cicoil_estimations.py`
3. `adios.oil_name_alias` → guarded with `getattr()`
4. `CONFIG_LEVEL_*` now imported from `opendrift.config` (module-level, not class attrs)
5. `max_water_fraction` class attr → `_cicoil_max_water_fraction_override`; `self.max_water_fraction = None` set in `__init__`
6. `self.get_environment()` → `self.env.get_environment()` (returns `(recarray, profiles, missing)`)
7. `self.Density` / `self.KinematicViscosity` initialized in NOAA-mode `seed_elements()`
8. Element ID indexing fixed: `elements.ID - 1` → `elements.ID` (4 occurrences)
9. `__set_seed_config__()` used for `seed:oil_type` (works in both Config and Ready modes — fixes `seed_plume_elements()` mode guard)
10. Mode.Config → Mode.Ready transition fixed in `seed_plume_elements()`
11. `pyproj` imported directly (not via `basereader.pyproj`)
12. `get_oil_budget()` fixed for xarray/masked-array via `_prop()` helper (`.values` + `np.nansum()`)
13. `required_profiles_z_range` removed; replaced by `_set_config_default('drift:profiles_depth', 20)`
14. `seed_plume_elements()` environment-finalization fallback (28°C default temperature before env finalized)

---

## 6. Photooxidation Module (`photooxidation_cicese()`)

Aromatic pseudo-component groups G1–G8 are photooxidized by solar UV when
near the surface. **Status: CLOSED 2026-06-13** — considered sufficient;
no further development planned absent new validation data.

### 6.1 Rate model
```
dM_i = K_PHOTO_AROM[i] * I_UV(z) * dt * M_i      (i = G1..G8, aromatic only)
I_UV(z) = I_UV(0) * exp(-K_D_UV * |z|)            (Beer-Lambert, depth-resolved)
```
- Elements deeper than `UV_DEPTH_CUTOFF = 50 m` are skipped.
- `K_D_UV = 0.10 m⁻¹` ("clear GoM" condition, from `surface_fate_coupled.py`).
- Aliphatic components and residue: unaffected (0.000 kg).
- Only acts at daytime (`I_UV(0) = 0` at night).

### 6.2 UV forcing
- **Preferred**: real ERA5 `surface_downward_uv_radiation` (CF-preprocessed
  via `preprocess_era5_uv.py`, data in `data_local/era5_uv_cf/`, 1972–2026).
  Real ERA5 UV peaks ~120 W/m² (vs the old synthetic ~50 W/m²).
- **Fallback**: synthetic diurnal cosine solar-zenith model (peak 50 W/m²,
  Bay of Campeche, used only if no UV reader is supplied).
- Requires `drift:vertical_mixing = True` — OpenDrift's 3D vertical mixing
  already represents the surface fraction; do **not** also apply a separate
  `f_surface` correction (the old `_compute_f_surface_cicese()` was removed
  as double-counting). With `vertical_mixing=False`, photooxidation is
  overstated by ~11%.

### 6.3 Current calibration (`K_PHOTO_AROM`, m² J⁻¹)
| Group | Compound | k (m² J⁻¹) | t½ |
|---|---|---|---|
| G1 | benzene | 2.9e-8 | ~14 d |
| G2 | toluene | 2.0e-8 | ~20 d |
| G3 | xylenes | 5.0e-8 | ~8 d |
| G4 | naphthalene | 4.260e-8 | ~9.4 d |
| G5 | acenaphthene | 5.145e-8 | ~7.8 d (interpolated G4–G6) |
| G6 | phenanthrene | 6.176e-8 | ~6.5 d |
| G7 | chrysene | 1.581e-8 | ~25.3 d (interpolated G6–G8) |
| G8 | pyrene+ | 5.110e-9 | ~79.8 d |

G4/G6/G8 are full-series least-squares recalibrated (f\*=0.3834) against
`surface_fate_coupled.py`'s 30-day reference (`data_local/surface_fate/surface_fate_base.csv`).
G5/G7 are log-linear interpolations vs. boiling point (no compound-specific
literature found; values fall within Bacosa et al. 2015 ring-class ranges).
A residual shape mismatch remains: the depth-resolved scheme saturates
~day 10–15 vs. the target's gradual 30-day rise (documented limitation, not
fixed).

### 6.4 Reference result & uncertainty
Nohoch Alfa 2023 (N=3000, 30 d, real forcing, MAYA): **total_oxy_frac ≈
5.73%** (~31.9 t for a 557 t spill). OAT uncertainty envelope (K_PHOTO_AROM
±50%, K_D_UV 0.05–0.20, DELTA_M_OX 16–48 g/mol): **2.85–9.60%** (~15.9–53.5 t,
0.5×–1.7× baseline) — conservative non-statistical bound. K_PHOTO_AROM
dominates the uncertainty (factor ~2.4).

**Interpretation**: 5.73% is a **lower-bound** estimate — aromatics (G1–G8)
only, no aliphatic photooxidation pathway implemented (cf. Ward et al. 2018
on radical-mediated aliphatic photooxidation).

### 6.5 New element variables
`mass_photooxidized` (kg), `fraction_photooxidized` (%). `mass_photooxidized`
is a **terminal sink** — not biodegraded further.

### 6.6 Oxygenated photoproduct (OP) tracking (audit #4, 2026-06-14)
Of each step's photooxidized mass `dM_phox` (aromatic groups G1–G8), a
fraction routes to a dissolved oxygenated-photoproduct pool instead of the
oil-phase `mass_photooxidized` sink:
```
dM_op      = OP_DISSOLVED_FRACTION * dM_phox        # -> mass_op_dissolved
dM_residue = (1 - OP_DISSOLVED_FRACTION) * dM_phox  # -> mass_photooxidized (as before)
```
- `OP_DISSOLVED_FRACTION = 0.5` — even 50/50 split, a first-pass estimate
  reflecting OPs (oxy-PAHs, quinones, carboxylic acids) being markedly more
  water-soluble than their parent PAHs (Aeppli et al. 2012); pending
  compound-class-resolved partitioning data.
- `mass_op_dissolved` decays at its own first-order rate
  `K_OP_REMOVAL = 0.116 d⁻¹` (= 2 × `K_BIODEG_WATER_ANCHOR`, t½≈6 d @ 20°C,
  Q10-scaled like `K_BIODEG_*`), into `mass_op_degraded`. This removal step
  is **unconditional** (runs even at night / `I_UV(0)=0`), since it
  represents a water-column process (further
  biodegradation/photodegradation/dilution of OPs) independent of solar UV.
- `fraction_photooxidized`, `fraction_op_dissolved`, `fraction_op_degraded`
  share the same conserved `initial_mass` denominator (now also including
  `mass_op_dissolved` + `mass_op_degraded`); the analogous denominator in
  `biodegradation_cicese()` (`fraction_biodegraded_from_oil`) was updated the
  same way for consistency.
- New element variables: `mass_op_dissolved`, `fraction_op_dissolved`,
  `mass_op_degraded`, `fraction_op_degraded` (all kg/%, default 0).
- Verified by `test_op_tracking.py`: 8-day surface MAYA run, mass
  conservation incl. new pools (rel. err 0), `mass_op_dissolved +
  mass_op_degraded == mass_photooxidized` (exact, by construction of the
  50/50 split), `mass_op_degraded > 0` (K_OP_REMOVAL decay active).

---

## 7. Biodegradation Module (`biodegradation_cicese()`)

Added 2026-06-14, replacing the previously dead/inconsistent pathway (old
code called inherited `biodegradation()` → `biodegradation_adcroft()`, which
decayed bulk `mass_oil` without updating `cicese_mass_balance['mass_components']`
— a mass-balance inconsistency; CICOIL's `mass_biodegraded_from_oil/_water`
tracking variables existed but were always zero).

### 7.1 Rate model
Component-resolved, first-order exponential decay applied per pseudo-component:
```
dM = M * (1 - exp(-K(T) * dt))
K(T) = K(T_REF) * BIODEG_Q10 ** ((T - T_REF) / 10)     # T_REF=20°C, BIODEG_Q10=3
```
- `sea_water_temperature` used per-element (profile-enabled — works for
  subsea droplets).
- Acts on **all active particles** (not depth/light-limited, unlike
  photooxidation).
- Surface particles (`z==0`) use `K * BIODEG_SURFACE_FACTOR` where
  **`BIODEG_SURFACE_FACTOR = 1/30`** — Prince et al. (2003, 2017) found
  floating slicks "almost immune to detectable biodegradation" vs. dispersed
  oil (1–3 week half-life), apparent surface half-life "many months to many
  years". `1/30` maps e.g. a 5 d dispersed-oil half-life → ~5 months surface.

### 7.2 Current calibration (`K_BIODEG_*`, d⁻¹ at T_REF=20°C, dispersed/entrained oil z<0)
| Group | Aliphatic k | t½ | Aromatic k | t½ |
|---|---|---|---|---|
| G1 | 0.1386 | 5 d | 0.0866 | 8 d |
| G2 | 0.1155 | 6 d | 0.0693 | 10 d |
| G3 | 0.0990 | 7 d | 0.0578 | 12 d |
| G4 | 0.0866 | 8 d | 0.0462 | 15 d |
| G5 | 0.0693 | 10 d | 0.0347 | 20 d |
| G6 | 0.0462 | 15 d | 0.0277 | 25 d |
| G7 | 0.0277 | 25 d | 0.0173 | 40 d |
| G8 | 0.0173 | 40 d | 0.0116 | 60 d |

Residue: `K_BIODEG_RESIDUE = 0.00139` (t½ ~500 d, effectively
non-biodegradable on simulation timescales). Surface particles: multiply all
of the above by `BIODEG_SURFACE_FACTOR = 1/30`.

Basis: Adcroft et al. (2010) Q10 law; Bacosa et al. (2015) relative
magnitudes (n-alkanes ~10× faster than PAH photooxidation; biodegradation
comparable-to-faster than photooxidation for heavy 4–5 ring PAHs); mesocosm
aliphatic half-lives 8–19 d. **First-pass, literature-informed — pending
recalibration against data** (same status as pre-recal K_PHOTO_AROM).

### 7.3 Dissolved-phase pathway (`mass_biodegraded_from_water`)
`handle_subsea_dissolution()` now maintains a persistent dissolved-mass pool
instead of immediate deactivation:
- On first encounter (`is_dissolved==1`): one-time transfer of
  `cicese_mass_balance['mass_components']` → `mass_dissolved` (components
  zeroed); `elements.mass_dissolved_subsea` overwritten with this reference
  total.
- Each step: `mass_dissolved` decays at **`K_BIODEG_WATER = 0.0578 d⁻¹`**
  (t½~12 d @ 20°C, Q10-scaled, **no** surface suppression — Adcroft et al.
  2010 dissolved-GoM-plume reference, also consistent with the fast end of
  Prince et al.'s 1–3 week dispersed-oil range).
- Element deactivates (`reason='dissolved_biodegraded'`) once
  `mass_dissolved < BIODEG_WATER_DEPLETION_FRAC (=0.01, i.e. 1%)` of the
  initial reference mass.
- First-pass single rate, not component-resolved.
- Only relevant for TAMOC subsea blowout scenarios
  (`seed_plume_elements()`); not exercised by surface releases
  (`seed_elements()`).

### 7.4 New/updated element variables
`mass_biodegraded_from_oil`, `mass_biodegraded_from_water`, generic
`mass_biodegraded` (for existing reporting/plots), `fraction_biodegraded_from_oil`,
`fraction_biodegraded_from_water`. Updates per-component
`cicese_mass_balance['mass_components']` (keeps component balance consistent
with `mass_oil`).

### 7.5 Reference results — Nohoch Alfa 2023 (N=3000, 30 d, real forcing, MAYA, surface release)
% of total initial oil mass at day 30, both `processes:biodegradation` and
`processes:photooxidation` enabled:

| Pathway | Pre-recal (no surface suppression) | Recalibrated (BIODEG_SURFACE_FACTOR=1/30) |
|---|---|---|
| Remaining | 3.0329% | 4.4324% |
| Evaporated | 73.8170% | 81.5155% |
| Photooxidized | 3.2750% | 4.2874% |
| **Biodegraded (from oil)** | **19.8751%** | **9.7648%** |
| Mass conservation rel. err | 0.00e+00 | 1.36e-16 |

Scaled to spill size: 557 t → 54.4 t biodeg / 23.9 t photoox; 1000 t → 97.6 t
/ 42.9 t.

### 7.6 Droplet-size-resolved biodegradation (audit #2, 2026-06-14)
The per-group rate matrix `K` (after Q10 and surface-suppression scaling) is
further multiplied by a surface-area-to-volume size factor:
```
K_i_eff = K_i * (BIODEG_DROPLET_D_REF / d) ** BIODEG_DROPLET_EXPONENT
```
- `BIODEG_DROPLET_D_REF = 0.001 m` (1 mm) — matches `seed:droplet_diameter_mu`,
  taken as the calibration reference diameter for `K_BIODEG_ALI/AROM/RESIDUE`
  (§7.2).
- `BIODEG_DROPLET_EXPONENT = 1.0` — SA/V ~ 1/d for a sphere (King 1992;
  Brakstad et al. 2015): smaller droplets present proportionally more
  oil-water interface per unit oil volume, degrading faster.
- `d = elements.diameter`, clamped to `BIODEG_DROPLET_DIAMETER_MIN/MAX =
  [1e-5, 1e-2] m` (10 µm – 1 cm), bounding `K_i_eff/K_i` to `[0.1, 100]`.
- `diameter == 0` (the default for most surface "spillet" particles, for
  which `diameter` is not a meaningful droplet size) is treated as `d_ref`,
  i.e. factor = 1 — no change, and no double-counting with
  `BIODEG_SURFACE_FACTOR` (§7.1), which remains the primary surface-vs-
  dispersed discriminator.

Net effect: subsea/dispersed droplets smaller than 1 mm (e.g. fine
wave-entrainment or subsea-blowout droplets) biodegrade faster than the §7.2
baseline; droplets larger than 1 mm biodegrade more slowly. Surface slicks
(`diameter==0`) are unaffected.

Verified by `test_droplet_size_biodeg.py` (90 particles, 4 groups of 30:
`diameter` = 0, 0.1 mm, 1 mm, 1 cm, 10 d @ SST=20°C=`BIODEG_T_REF`): mass
conservation rel. err < 1e-6;
`frac_biodeg(0.1mm)=7.110% > frac_biodeg(1mm)=0.810% > frac_biodeg(1cm)=0.082%`;
`frac_biodeg(diameter=0) == frac_biodeg(diameter=1mm)` (rel. diff 4.7e-8).

---

## 8. Emulsification Module (`emulsification_cicese()`)

### 8.1 Problem with the stock NOAA Bullwinkle rate

The upstream NOAA PyGNOME Bullwinkle water-uptake model computes a single
characteristic timescale `t_char = (6/drop_max) / k_emul_eff` where, at a
typical 10 m/s wind, the unscaled `k_emul ≈ 121 s⁻¹` gives `t_char ≈ 1–5
hours` — i.e. the oil forms a stable mousse within a few hours regardless of
its viscosity or weathering state. For MAYA crude (and most medium-to-light
crudes), field observations and Fingas & Fieldhouse (2004) indicate that
tar-ball/mousse formation takes **days to weeks**, not hours.

### 8.2 Viscosity-stability scaling

The fix introduces a viscosity-dependent stability factor before the NOAA
Bullwinkle water-uptake rate is applied:

```
k_emul_eff = k_emul * clip(nu_dry / EMUL_STABILITY_NU_REF,
                            EMUL_STABILITY_MIN, EMUL_STABILITY_MAX)
           * EMUL_RATE_SCALE
```

where:
- `nu_dry` = dry-oil (no-emulsion), evaporation-corrected kinematic viscosity
  [cSt] (`self._oil_viscosity_dry`, computed each timestep in
  `oil_weathering_cicese()`).
- `EMUL_STABILITY_NU_REF = 1000.0 cSt` — reference viscosity around which
  the stability factor is centred (clip argument = 1.0 at `nu_dry = 1000 cSt`).
- `EMUL_STABILITY_MIN = 0.3` — minimum stability factor (fresh, low-viscosity
  oil emulsifies more slowly).
- `EMUL_STABILITY_MAX = 3.0` — maximum stability factor (very viscous oil
  emulsifies faster, up to 3×).
- `EMUL_RATE_SCALE = 0.0142` — overall multiplicative correction (≈2 orders
  of magnitude below the raw NOAA rate) needed because `EMUL_STABILITY_MIN/MAX`
  alone (a ≤10× range) cannot supply the required ~100× rate reduction.

To recover the original (viscosity-independent) NOAA Bullwinkle behaviour
exactly, set `EMUL_STABILITY_MIN = EMUL_STABILITY_MAX = EMUL_RATE_SCALE = 1.0`.

### 8.3 Empirical calibration

Calibrated against a MAYA crude Abkatun 2026 TAMOC+GLORYS12 realistic run
(`calibrate_emulsification_viscosity_tamoc.py`) using a 17.2-day surface tar-ball
transition target from `TarballFormation_MPB.tex` (Fingas & Fieldhouse 2004
consistent):

| `EMUL_RATE_SCALE` | Tar-ball crossing day (ν ≥ 10,000 cSt) |
|---|---|
| 1.0 (no scaling) | 0.76 d |
| 0.024 | 10.21 d |
| **0.0142** | **16.92 d** ← selected (target 17.2 d, error <1.6%) |
| 0.0072 | never within 25 d (ν = 7,388 cSt @ day 25) |

### 8.4 Impact on case-study mass balances

The emulsification fix is the **dominant change** in the Nohoch 2023 and
Abkatun 2026 re-runs (June 2026):

- **Pre-fix**: every particle reached `max_water=0.9` (the ADIOS cap) and
  viscosity 18–32× above the tar-ball threshold within ~5 hours of release,
  in every case.
- **Post-fix**: water content builds gradually over days, reaching
  `water_fraction` ≈ 0.57–0.60 and viscosity ≈ 5,500–9,300 cSt by the end
  of 22–30 day runs (approaching but not crossing the 10,000 cSt threshold).

Key quantitative effects (see `findings_emulsification_fix_audit24.md`):

| Case | Evaporated (pre/post) | Biodegraded (pre/post) |
|---|---|---|
| Nohoch 30 d | 85.1% / 94.1% | 5.36% / 0.98% |
| Abkatun 22 d (surface seed) | 81.9% / 91.5% | 8.02% / 1.01% |
| Abkatun TAMOC-init (d₀=0.56 mm) | 61.2% / 74.6% | 32.0% / 19.5% |

The mechanism: pre-fix, the rapid emulsification locked most oil mass into a
high-water emulsion, suppressing evaporation and making more oil available
for biodegradation. The fix restores physically realistic evaporation as the
dominant fate process for surface-released MAYA crude in the Bay of Campeche.

### 8.5 Unit test

`test_emulsification_viscosity.py` (4/4 PASS):
- cold oil (5°C, high viscosity → stability > 1.0) emulsifies faster than
  warm oil (30°C, lower viscosity → stability < 1.0)
- with `EMUL_STABILITY_MIN=MAX=1.0`, cold==warm (viscosity-independence
  recovered)
- default stability always < unit stability (scaling suppresses rate for
  fresh MAYA)
- `EMUL_RATE_SCALE=1.0` gives stock NOAA rate (much faster than default)

---

## 9. Testing & Verification

| Test | What it checks | Result |
|---|---|---|
| `test_cicoil_e2e.py` | 30-test full suite (T1–T8): seeding, weathering, TAMOC plume seeding, mass/component balance, biodegradation (T7), photooxidation + OP tracking + regression bound (T8) | 30/30 PASS, mass conservation 100.000%, component balance 0.000% error |
| `test_photoox_vertmix.py` | `vertical_mixing` on/off vs `surface_fate_coupled.py` reference | `vertical_mixing=True`: 3.50% vs 3.47% target (within 1%); `False`: overstates by ~11% |
| `test_biodegradation.py` | Q10 sanity (surface release, 10 d, MAYA, all other weathering off) | Post-recal: 0.273%/0.810%/2.357% biodegraded at 10/20/30°C (vs 7.110%/16.595%/29.301% pre-recal, ~30× smaller as expected); mass conservation exact |
| `test_biodegradation_water.py` | Dissolved-pool decay/deactivation (20 particles, `is_dissolved=1`, 60 d, SST 10/20/30°C) | 10°C: 68.526% biodeg, 0/20 deactivated; 20°C: 96.882%, 0/20; 30°C: 99.043%, 20/20 deactivated; mass conservation ~1e-15 |
| `test_droplet_size_biodeg.py` | Droplet-size-resolved biodegradation (audit #2; 90 particles, 4 diameter groups, 10 d, SST=20°C) | `frac_biodeg(0.1mm)=7.110% > frac_biodeg(1mm)=0.810% > frac_biodeg(1cm)=0.082%`; `diameter=0` == `diameter=1mm` (rel. diff 4.7e-8); mass conservation rel. err 1.22e-16 |
| `test_op_tracking.py` | Oxygenated photoproduct tracking (audit #4; 100 particles, 8 d, SST=29°C, synthetic UV) | mass conservation incl. OP pools rel. err 1.22e-16; `mass_op_dissolved + mass_op_degraded == mass_photooxidized` (OP_DISSOLVED_FRACTION=0.5, exact); `mass_op_degraded > 0` (K_OP_REMOVAL decay active) |
| `test_emulsification_viscosity.py` | Viscosity-stability scaling: cold oil emulsifies faster than warm; unit-stability recovers viscosity independence; default stability suppresses rate vs unity; EMUL_RATE_SCALE=1.0 recovers stock NOAA rate | 4/4 PASS |

| `test_croco_native_reader.py` | `reader_croco_native` on a synthetic CROCO file: derived `mask_u`/`mask_v`, de-staggered `u`/`v` at a rho node equal the analytic value, `destagger=False` shows the stock half-cell offset | 7/7 PASS |

### Known caveats
- **Emulsification**: `max_water=0.9` required at `seed_elements()` for MAYA
  (ADIOS AD02254 has no emulsification data → `emulsification_cicese()` exits
  immediately otherwise). Rate is calibrated for MAYA crude at 10 m/s wind
  (17.2-day tar-ball target); re-calibrate `EMUL_RATE_SCALE` for other oils.
- **ANFC vs GLORYS12**: ANFC ocean forcing gives non-reproducible results
  across download dates (80 km vs 151 km transport for the same nominal
  run). Use GLORYS12 (`--ocean glorys`) for reproducible hindcast
  attribution.
- **Longitude sign**: Bay of Campeche is °W — use `lon=-92.x`, not `+92.x`
  (positive longitude silently places particles in the Indian Ocean).
- **Mass scaling**: CICOIL enforces 185.6 kg/element; `m3_per_element` is
  rejected. Scale outputs by `M_spill / (N_p × 185.6 kg)`.

---

## 10. File Manifest (HPC deployment tarball)

**Core model:**
- `ciceseoil.py` (142 KB) — `OpenCiceseOil` class, 14 API fixes, photooxidation + biodegradation + emulsification fix
- `fluid_properties.py` (19 KB) — 17-component distillation framework (Cox Pv + Mackay)
- `cicoil_estimations.py` (13 KB) — standalone ADIOS-DB estimation workarounds

**TAMOC coupling:**
- `tamoc_plume.py`, `run_tamoc.py`, `tamoc_chemical_properties.py`

**NEMO/CROCO readers** (`readers/`):
- `reader_NEMO_native_v3.py`, `reader_nemo_combined.py`, `reader_nemo_modified.py`, `reader_nemo_optimized.py`
- `reader_croco_native.py` — raw CROCO sigma-coordinate output, de-staggered (subclass of OpenDrift `reader_ROMS_native`)

**Export:**
- `export/io_stat_nc.py` — gridded statistical NetCDF output

**Data** (`data/`, 5 CSVs, ~200 KB):
- `pseudo_chemdata.csv`, `ChemData.csv`, `Aij.csv`, `Bij.csv`, `gas_composition.csv`

**Testing:**
- `test_cicoil_e2e.py` (30-test suite, T1–T8), `test_biodegradation.py`,
  `test_biodegradation_water.py`, `test_droplet_size_biodeg.py`,
  `test_op_tracking.py`, `test_emulsification_viscosity.py`

**Shared helpers:**
- `cicoil_common.py` — `load_ciceseoil()` (importlib boilerplate) and
  `default_cicese_config()` (standard process on/off settings)

**Utilities:**
- `download_data.py` (GLORYS12/ERA5 download helpers)
- `preprocess_era5_uv.py`, `preprocess_real_forcing.py`
- `plotting.py` (10+ plotting routines)

**Deployment:**
- `CICOIL_v2_MANUAL.md` (this file), `stage_for_hpc.sh`, `deploy_on_hpc_conda.sh`,
  `slurm_template_cicoilv2.sh`, `DEPLOYMENT_MANIFEST.txt`

Total: ~500 MB tarball, ~150 MB installed in site-packages.

---

## 11. References

- Adcroft et al. (2010) — Q10=3, T_REF=20°C temperature-scaling law for
  biodegradation; dissolved GoM plume reference rate
- Bacosa et al. (2015), Marine Pollution Bulletin — ring-class
  biodegradation/photooxidation half-life ranges
- Prince et al. (2003, 2017), ACS ES&T "The Rate of Crude Oil Biodegradation
  in the Sea" — surface-slick vs dispersed-oil biodegradation contrast
  (basis for `BIODEG_SURFACE_FACTOR`)
- Delvigne & Sweeney (1988), Oil & Chem. Pollut. 4:281 — wave-submergence
  (historical; superseded by `vertical_mixing`-based surface fraction)
- Ward et al. (2018) — radical-mediated aliphatic photooxidation (discussed,
  not implemented — aromatics-only is a lower bound)
- OpenDrift 1.14.9 documentation
- CICOIL 2023 technical reports (Kotzakoulakis, CICESE/SINTEF Ocean)

**Session logs:**
- `session_2026-06-05_CICOIL_analysis.{html,pdf}` — full code audit
- `session_2026-06-06_CICOIL_install_test.{html,pdf}` — install & testing
- `session_2026-06-13_era5_uv_recalibration.{html,tex}` — photoox UV recalibration
- `session_2026-06-14_CICOILv2_Biodegradation.{tex,pdf}` (12 pp) — biodegradation module
- `BIODEGRADATION.md` — full biodegradation diagnostic + rate tables

---

**Version:** CICOILv2.0
**Last updated:** 2026-06-14
