# CICOILv2 Biodegradation Module — Notes (2026-06-14)

## 1. Background: CICOILv1 → CICOILv2 rename

Before this work, the active development copy was renamed repo-wide from
`CICOILv1` to `CICOILv2` (directories, files, ~118 text references across
source, scripts, manuscripts, and session logs):

- `cicoilv1_deployment/` → `cicoilv2_deployment/` (this directory, active dev copy)
- `cicoilv1/` → `cicoilv2/`
- `abkatun_2026/data_local/{forward,backward,figures}_cicoilv1/` → `..._cicoilv2/`
- Individual files: `slurm_template_cicoilv1.sh`, `cicoilv1_hpc_deployment.tar.gz`,
  `Tarball_formation_CICOILv1_0.pdf`, `backward_full_cicoilv1_test.log`,
  `backward_simulation_cicoilv1.py`, `forward_simulation_cicoilv1.py`, and
  figure PNGs — all renamed to `v2` equivalents.
- All case-variant strings (`CIC-OILv1`, `CICOIL v1`, `CICOILv1`, `cicoilv1`)
  replaced with v2 equivalents via `sed` across 118 text files.

No behavioral changes from the rename itself — `ciceseoil.py` content was
untouched until the biodegradation work below.

## 2. Diagnostic: the pre-existing biodegradation pathway was dead/inconsistent

Before this change, `oil_weathering_cicese()` called the **inherited upstream**
`self.biodegradation()` (from `openoil.py`), which dispatches to
`biodegradation_adcroft()` based on `biodegradation:method`.

Problems identified:

1. **Bulk-only, not component-resolved.** `biodegradation_adcroft()` decays
   the bulk `mass_oil` / generic `mass_biodegraded` element variables using a
   single temperature-dependent timescale
   `tau(T) = 12 * 3^((20-T)/10)` days (Adcroft et al. 2010, GoM dissolved-oil
   plumes — 62 d @ 5°C, 12 d @ 20°C, 7 d @ 25°C).

2. **Does not touch `cicese_mass_balance['mass_components']`** — the
   17-pseudo-component (8 aliphatic G1-G8, 8 aromatic G1-G8, 1 residue) array
   that evaporation and photooxidation read/write. After biodegradation ran,
   `sum(mass_components, axis=1) != mass_oil` for affected elements — a
   mass-balance inconsistency between the bulk and component-resolved
   accounting.

3. **CICOIL's own split-tracking variables were unused scaffolding.**
   `mass_biodegraded_from_oil` / `mass_biodegraded_from_water` (element
   variables) and the corresponding `cicese_mass_balance` arrays existed,
   were initialized to zero in `prepare_run()`, but were **never written to**
   — always zero regardless of `processes:biodegradation`.

Net effect: enabling `processes:biodegradation` previously removed mass from
the bulk oil budget without keeping the per-component composition consistent,
and the dedicated CICOIL biodeg-tracking fields were silently always zero.

## 3. New implementation: `biodegradation_cicese()`

Replaces the call in `oil_weathering_cicese()` (was `self.biodegradation()` →
now `self.biodegradation_cicese()`), structurally mirroring
`photooxidation_cicese()`.

### Rate model

- **Temperature scaling**: reuses the Adcroft et al. (2010) Q10=3,
  T_ref=20°C law: `rate(T) = rate(T_ref) * Q10^((T-T_ref)/10)`, applied
  per-element using local `sea_water_temperature` (profile-enabled —
  works for subsea droplets too).
- **Decay**: first-order exponential, `dM = M * (1 - exp(-K*dt))`, applied
  per pseudo-component, for numerical stability at large `dt` (consistent
  with the upstream Adcroft/half_time formulas).
- **Scope**: oil-phase only (`mass_biodegraded_from_oil`), acts on **all
  active particles** (surface slick + subsea droplets) — biodegradation is
  not light-limited, unlike photooxidation which is depth/UV-gated.

### Rate constants (`K_BIODEG_*` in `ciceseoil.py`)

| Group | Aliphatic k (d⁻¹) | t½ (d) | Aromatic k (d⁻¹) | t½ (d) |
|-------|------------------|--------|------------------|--------|
| G1 | 0.1386 | 5  | 0.0866 | 8  |
| G2 | 0.1155 | 6  | 0.0693 | 10 |
| G3 | 0.0990 | 7  | 0.0578 | 12 |
| G4 | 0.0866 | 8  | 0.0462 | 15 |
| G5 | 0.0693 | 10 | 0.0347 | 20 |
| G6 | 0.0462 | 15 | 0.0277 | 25 |
| G7 | 0.0277 | 25 | 0.0173 | 40 |
| G8 | 0.0173 | 40 | 0.0116 | 60 |

Residue: `K_BIODEG_RESIDUE = 0.00139` (t½ ~500 d, essentially
non-biodegradable on simulation timescales).

Basis: Adcroft et al. (2010) Q10 temperature law; Bacosa et al. (2015, MPB)
relative magnitudes (n-alkane biodegradation ~10x PAH photooxidation rates;
biodegradation comparable-to-faster than photooxidation for heavy 4-5 ring
PAHs); mesocosm aliphatic half-life ranges (8-19 d bulk, C10-C22 most
bioavailable, heavy waxes slower).

**Status (as of the recalibration in §6)**: the table above gives the
**dispersed/entrained-oil** (z<0) rates. Surface-slick particles (z==0) use
these same rates multiplied by `BIODEG_SURFACE_FACTOR=1/30` (see §6) —
Prince et al. (2003, 2017) report floating slicks are "almost immune to
detectable biodegradation" relative to dispersed oil.

### Config updates

- `processes:biodegradation` description rewritten to describe the new
  component-resolved CICESE scheme.
- `biodegradation:method` description rewritten — `'Adcroft'` (default) now
  refers to the Q10 temperature-scaling law applied per-component inside
  `biodegradation_cicese()`, not the old upstream bulk dispatch. Enum/default
  unchanged (no behavioral dispatch needed — `biodegradation_cicese()` is the
  only path now).

`mass_biodegraded_from_water` (dissolved-phase pathway) is implemented as of
§7 below — see that section for details (it was originally flagged here as
out-of-scope future work).

## 4. Testing

- **`test_biodegradation.py`** (new): MAYA oil, 50 particles, 10-day sim,
  calm constant forcing, all other weathering off. Checks:
  - Mass conservation (`mass_components + mass_evaporated + mass_dissolved +
    mass_photooxidized + mass_biodegraded_from_oil + mass_biodegraded_from_water
    == initial`): rel. err ~0 (1e-16 to 0.0).
  - Q10 sanity check at SST = 10/20/30°C: 7.110% / 16.595% / 29.301%
    biodegraded over 10 days — monotonically increasing with temperature. PASS.
  - `mass_biodegraded_from_water == 0` for all cases (expected). PASS.
- **`test_cicoil_e2e.py`** (existing 20-test suite): 20/20 pass, no
  regression (biodegradation defaults to `False`).

## 5. Nohoch Alfa 2023 real-forcing test (in progress)

`run_nohoch30d_biodeg_realforcing_3000p.py` — adapted from
`run_nohoch30d_photoox_realforcing_3000p_calib2.py`: same real forcing
(ERA5 UV depth-resolved + Mercator currents/SST/salinity + ERA5 wind/waves,
N=3000, 30 days, release Jul6 2023 at LON=-92.0/LAT=19.4, MAYA oil), but with
`processes:biodegradation=True` **and** `processes:photooxidation=True`
together (combined weathering effect on the validated scenario).

Reports, per pseudo-component (17 groups) and in aggregate:
- `mass_biodegraded_from_oil` as fraction of total initial oil mass
- `mass_photooxidized` (aromatics G1-G8) as fraction of total initial oil mass
- Full mass balance breakdown (remaining/evaporated/dissolved/photoox/biodeg)
  and mass-conservation relative error

Output: `data_local/nohoch2023_biodeg_estimate/cicoil_nohoch30d_biodeg_realforcing_3000p.nc`
and `.npz` summary.

### Results, pre-recalibration (2026-06-14, N=3000, 30 days, real forcing)

These were the first-pass results, using dispersed-oil-calibrated
`K_BIODEG_*` applied uniformly to ALL particles (no surface suppression).
**Superseded by the recalibrated results in §6** — kept here for comparison.

Mass balance, day 30 (fraction of total initial oil mass):

| Pathway | % of total oil |
|---|---|
| Remaining (mass_components) | 3.0329% |
| Evaporated | 73.8170% |
| Dissolved | 0.0000% |
| Photooxidized (aromatics G1-G8) | 3.2750% |
| **Biodegraded (from oil, all 17 components)** | **19.8751%** |
| Biodegraded (from water) | 0.0000% (future work) |
| **Mass conservation rel. err** | **0.00e+00 (exact)** |

Per-component biodegraded mass (% of total initial oil mass):

| Group | Aliphatic | Aromatic |
|---|---|---|
| G1 | 0.00000% | 0.00000% |
| G2 | 0.31971% | 0.08840% |
| G3 | 0.61077% | 0.18119% |
| G4 | 1.18437% | 0.51069% |
| G5 | 3.47104% | 1.94491% |
| G6 | 1.37270% | 0.90555% |
| G7 | 2.31819% | 1.78017% |
| G8 | 0.91881% | 0.76291% |
| Residue | 3.50573% | — |

Scaled to spill sizes:

| Total spill mass | Biodegraded | Photooxidized |
|---|---|---|
| 557 t (3000 x 185.6 kg, CICOILv2 default) | 110.7 t | 18.2 t |
| 1000 t (illustrative) | 198.8 t | 32.8 t |

### Interpretation (pre-recalibration)

- Biodegradation removes ~6x more mass than photooxidation over 30 days
  (19.9% vs 3.3% of total oil) — consistent with Bacosa et al. (2015)'s
  finding that biodegradation rates generally exceed photooxidation rates,
  especially for the lighter/more bioavailable pseudo-components.
- G1 (lightest aliphatic and aromatic/benzene) shows **0% biodegraded** —
  these components evaporate essentially completely (within hours to ~1 day)
  before the slower biodegradation timescales (t½ 5-8 d) can act on
  meaningful remaining mass; evaporation dominates for G1 at Campeche SSTs
  (~29°C).
- G5 (both aliphatic and aromatic) shows the largest biodegraded fraction —
  reflects both a moderate intrinsic rate (t½ 10/20 d) and a large remaining
  pool after evaporation removes the lighter groups first.
- Residue (3.5% of total mass biodegraded) is non-trivial despite the very
  slow base rate (t½ ~500 d) because of the warm Campeche SST (~29°C):
  Q10=3 scaling gives `rate(29C) = rate(20C) * 3^0.9 ≈ 2.7x`, i.e.
  effective t½ ≈ 185 d, giving ~10.6% of the *residue pool* biodegraded over
  30 days (≈3.5% of total oil mass; residue is **64.2%** of total initial
  mass for MAYA — the prior estimate of "~1/3" here was incorrect).
- These first-pass `K_BIODEG_*` constants applied uniformly (no surface
  suppression) over-predict biodegradation for a surface slick — see §6.

## 6. Recalibration (2026-06-14): surface-slick suppression (`BIODEG_SURFACE_FACTOR`)

### Motivation

The pre-recalibration result (§5) applies the dispersed-oil-calibrated
`K_BIODEG_*` (basis: mesocosm/dispersed-oil half-lives, §3) to **all**
particles, including surface-slick particles (z==0) — which dominate the
Nohoch Alfa scenario (a surface release). This over-predicts biodegradation:
floating slicks have far less surface-area-to-volume for microbial
colonization than dispersed droplets.

### Literature basis: Prince et al. (2003, 2017)

**Prince et al., Environ. Sci. Technol., "The Rate of Crude Oil
Biodegradation in the Sea"**: dispersed oil at environmentally relevant
(sub-ppm) concentrations has an apparent half-life of **1-3 weeks**, whereas
floating surface slicks are "almost immune to detectable biodegradation" with
apparent half-lives of **"many months to many years"** — roughly a 1-2 order
of magnitude difference, driven by the surface-area-to-volume ratio available
for microbial attack.

### Implementation: two-regime model

`ciceseoil.py`, `biodegradation_cicese()`:

```python
BIODEG_SURFACE_FACTOR = 1.0 / 30.0

is_surface = self.elements.z[active] == 0
surface_mult = np.where(is_surface, self.BIODEG_SURFACE_FACTOR, 1.0)
K = K * surface_mult[:, np.newaxis]
```

- `K_BIODEG_ALI`/`K_BIODEG_AROM`/`K_BIODEG_RESIDUE` (§3 table, t½ 5-60 d /
  ~500 d) now apply to **dispersed/entrained oil** (z<0).
- Surface-slick particles (z==0) get these rates × 1/30 → effective t½
  **~5 months to ~5 years**, spanning Prince et al.'s reported range.
- `1/30` was chosen as a round number mapping the dispersed-oil half-life
  range onto the low end of "months to years" (5 d × 30 = 150 d ≈ 5 months;
  60 d × 30 = 1800 d ≈ 5 years).

### Verification

- `test_biodegradation.py` (surface release, z=0, 10 days): fraction
  biodegraded dropped from 7.110%/16.595%/29.301% (10/20/30°C,
  pre-recalibration) to **0.273%/0.810%/2.357%** (≈30x smaller, consistent
  with `BIODEG_SURFACE_FACTOR=1/30`). Mass conservation rel. err = 0 for all.
  PASS.
- `test_cicoil_e2e.py`: 20/20 PASS, no regression.

### Recalibrated Nohoch Alfa 2023 results (N=3000, 30 days, real forcing)

Mass balance, day 30 (fraction of total initial oil mass):

| Pathway | Pre-recalibration | **Recalibrated** |
|---|---|---|
| Remaining (mass_components) | 3.0329% | **4.4324%** |
| Evaporated | 73.8170% | **81.5155%** |
| Dissolved | 0.0000% | **0.0000%** |
| Photooxidized (aromatics G1-G8) | 3.2750% | **4.2874%** |
| **Biodegraded (from oil, all 17 components)** | **19.8751%** | **9.7648%** |
| Biodegraded (from water) | 0.0000% | **0.0000%** (dissolved pool not seeded in this scenario; see §7) |
| **Mass conservation rel. err** | 0.00e+00 | **1.36e-16** |

Per-component biodegraded mass (% of total initial oil mass), recalibrated:

| Group | Aliphatic | Aromatic |
|---|---|---|
| G1 | 0.00000% | 0.00000% |
| G2 | 0.05603% | 0.01634% |
| G3 | 0.15794% | 0.04203% |
| G4 | 0.38270% | 0.14859% |
| G5 | 1.58088% | 0.81048% |
| G6 | 0.79434% | 0.44804% |
| G7 | 1.34292% | 0.98857% |
| G8 | 0.53214% | 0.43591% |
| Residue | 2.02785% | — |

Scaled to spill sizes:

| Total spill mass | Biodegraded | Photooxidized |
|---|---|---|
| 557 t (3000 x 185.6 kg, CICOILv2 default) | 54.4 t | 23.9 t |
| 1000 t (illustrative) | 97.6 t | 42.9 t |

### Interpretation (recalibrated)

- Total biodegraded mass drops from 19.9% to **9.8%** of total oil — about
  **half** of the pre-recalibration estimate, reflecting that the Nohoch Alfa
  scenario is a surface release where most particles spend most of the
  simulation at z=0.
- Photooxidized mass correspondingly **rises** from 3.3% to **4.3%** — see §8
  (comparison vs. photooxidation-only) for the mechanism: less
  biodegradation-driven depletion of the aromatic pool leaves more mass
  available for photooxidation.
- The relative pattern across pseudo-components is unchanged (G5/G7 still
  dominate, G1 still 0%) — only the overall magnitude is rescaled (not
  uniformly, since particles with more subsurface/dispersed time during the
  30 days retain a larger dispersed-rate contribution).
- Mass conservation remains exact to numerical precision (1.36e-16).

## 7. Dissolved-phase biodegradation (`mass_biodegraded_from_water`)

### Motivation

Previously (§3/§5) `mass_biodegraded_from_water` was flagged as out-of-scope:
`is_dissolved==1` particles (seeded by `seed_plume_elements()` for TAMOC
far-field dissolved plume mass) had their `mass_dissolved_subsea` set once at
seeding and were **deactivated the same step** by
`handle_subsea_dissolution()` — called *before* `biodegradation_cicese()` in
`oil_weathering_cicese()` — leaving no persistent dissolved-mass pool for a
from-water decay to act on.

### Implementation: persistent dissolved pool + lifecycle

**`handle_subsea_dissolution()`** (rewritten): instead of immediate
deactivation,
1. On the first encounter (`mass_dissolved` sum == 0), transfers each
   element's `cicese_mass_balance['mass_components']` into
   `cicese_mass_balance['mass_dissolved']` (zeroing `mass_components`), and
   overwrites `elements.mass_dissolved_subsea` with the transferred total —
   making it a self-consistent "initial dissolved mass" reference.
2. On later steps, checks whether `mass_dissolved` has decayed below
   `BIODEG_WATER_DEPLETION_FRAC` (1%) of that reference, and if so,
   deactivates the element (`reason='dissolved_biodegraded'`).

**`biodegradation_cicese()`** (extended): for `is_dissolved==1` elements,
decays `cicese_mass_balance['mass_dissolved']` at rate `K_BIODEG_WATER`
(Q10-scaled, no surface suppression — dissolved mass is never at the
surface), moving mass into `mass_biodegraded_from_water` and
`mass_biodegraded`, and updates `fraction_biodegraded_from_water` (relative
to `elements.mass_dissolved_subsea`).

### Rate constant

```python
K_BIODEG_WATER_ANCHOR = 0.0578  # d^-1, t1/2 ~ 12 d @ 20C (Adcroft et al. 2010)
K_BIODEG_WATER = (
    K_BIODEG_WATER_ANCHOR
    / np.mean(np.concatenate([K_BIODEG_ALI, K_BIODEG_AROM]))
    * np.concatenate([K_BIODEG_ALI, K_BIODEG_AROM, [K_BIODEG_RESIDUE]])
).astype(np.float32)  # (17,)
BIODEG_WATER_DEPLETION_FRAC = 0.01
```

Adcroft et al. (2010) report a 12-day apparent half-life for dissolved GoM
oil-plume hydrocarbons at T_ref=20°C. Rather than apply that single rate
uniformly, `K_BIODEG_WATER` is now **component-resolved** (17 values,
matching the G1-G8 aliphatic, G1-G8 aromatic, and residue layout of
`cicese_mass_balance`'s pseudo-components): the dispersed-oil
`K_BIODEG_ALI`/`K_BIODEG_AROM`/`K_BIODEG_RESIDUE` table from §3 is rescaled by
a single factor (≈0.972) so that the mean rate across the 16 hydrocarbon
groups equals the Adcroft anchor (0.0578 d⁻¹). This preserves the
lighter/more-soluble-degrades-faster ordering already established for
dispersed oil while keeping the ensemble-mean rate consistent with the
Adcroft observation. Residue is rescaled by the same factor and remains
effectively non-biodegradable in the dissolved phase too (t½ > 500 d even at
30°C after Q10 scaling).

### Verification: `test_biodegradation_water.py` (new)

MAYA oil, 20 particles seeded directly with `is_dissolved=1`,
`mass_dissolved_subsea=1.0` kg, 60-day sim, calm conditions, SST = 10/20/30°C:

| SST | Fraction biodegraded (water) after 60 d | Active / deactivated |
|---|---|---|
| 10°C | 24.278% | 20 / 0 |
| 20°C | 36.493% | 20 / 0 |
| 30°C | 49.124% | 20 / 0 |

- **Mass conservation**: `mass_dissolved + mass_biodegraded_from_water +
  mass_components` == initial `sum(mass_components)` to within 3.39e-16
  rel. err, for all SSTs. PASS.
- **One-time transfer**: `mass_components` → 0 after the first step for all
  is_dissolved particles. PASS.
- **Q10 sanity check**: 24.3% < 36.5% < 49.1% — monotonically increasing
  with temperature. PASS.
- **Long-lived residue tail**: with the component-resolved rate table, the
  dissolved residue component has a very long half-life (mirroring
  `K_BIODEG_RESIDUE`'s role for dispersed oil in §3), so the total dissolved
  pool never drops below the 1% `BIODEG_WATER_DEPLETION_FRAC` threshold within
  60 days at any tested temperature — `n_deactivated == 0` everywhere is
  expected and physically correct (analogous to dispersed-oil residue
  persisting indefinitely). The test instead checks that the 30°C
  biodegraded fraction stays well below 100% (residue persists). PASS.
- `test_cicoil_e2e.py`: 20/20 PASS, no regression (the plume-seeding test
  T4.1-T4.3 exercises `is_dissolved==1` particles through the new lifecycle).

### Status

Not yet exercised in the Nohoch Alfa scenario (§6): that is a surface
release via `seed_elements()`, so no `is_dissolved==1` particles are ever
created (`mass_biodegraded_from_water` correctly stays at 0%, as shown in
§6's table). The pathway is relevant for **subsea blowout** scenarios seeded
via `seed_plume_elements()` (TAMOC far-field dissolved plume mass), which is
the configuration verified by `test_cicoil_e2e.py` T4.1-T4.3 and
`test_biodegradation_water.py` above. `K_BIODEG_WATER` is now
component-resolved (rescaled from the §3 `K_BIODEG_*` table) but the rescale
factor itself is still a first-pass estimate — same caveat status as
`K_BIODEG_*` in §3.

## 8. Audit improvements #2 and #4 (2026-06-14)

Following the CICOILv2 framework audit (13 ideas logged on the Research
Ideas Board, IDEA-004..013), improvements #2 and #4 were implemented.

### #2 — Droplet-size-resolved biodegradation

`biodegradation_cicese()` now applies a surface-area-to-volume scaling
factor to `K` (the per-group, Q10-scaled, surface-suppression-adjusted rate
matrix), using `self.elements.diameter`:

```python
K_i_eff = K_i * (BIODEG_DROPLET_D_REF / d) ** BIODEG_DROPLET_EXPONENT
```

- `BIODEG_DROPLET_D_REF = 0.001` m (1 mm) — matches
  `seed:droplet_diameter_mu`, taken as the calibration reference diameter for
  `K_BIODEG_ALI/AROM/RESIDUE`.
- `BIODEG_DROPLET_EXPONENT = 1.0` — from SA/V ~ 1/d for a sphere (King 1992;
  Brakstad et al. 2015): smaller droplets present proportionally more
  oil-water interface per unit oil volume, degrading faster.
- `d` is clamped to `[BIODEG_DROPLET_DIAMETER_MIN, BIODEG_DROPLET_DIAMETER_MAX]
  = [1e-5, 1e-2]` m (10 µm – 1 cm), bounding the factor to `[0.1, 100]` —
  avoids unphysical rates for sub-resolution wave-entrained droplets (down to
  ~1 µm) or unusually large subsea droplets.
- `diameter == 0` (the default for most surface "spillet" particles, for
  which `diameter` is not a meaningful droplet size) is treated as `d_ref`,
  i.e. factor = 1 — no change for these particles, and no double-counting
  with the existing surface suppression (`BIODEG_SURFACE_FACTOR`, §6), which
  remains the primary surface-vs-dispersed discriminator.

Net effect: subsea/dispersed droplets smaller than 1 mm (e.g. the
wave-entrainment spectrum down to ~1 µm, or fine subsea-blowout droplets)
biodegrade faster than the §3/§6 baseline; droplets larger than 1 mm
biodegrade more slowly. Surface slicks (`diameter==0`) are unaffected by this
change.

### #4 — Oxygenated photoproduct (OP) tracking

Mirrors the `mass_dissolved`/`K_BIODEG_WATER` pattern of §7. Two new pools
are added to `cicese_mass_balance` (and per-particle element variables
`mass_op_dissolved`, `mass_op_degraded`, `fraction_op_dissolved`,
`fraction_op_degraded`), both 17-component arrays restricted in practice to
the aromatic indices 8–15.

`photooxidation_cicese()` is restructured:

1. **OP removal (unconditional, runs even at night)**: the dissolved OP pool
   decays at `K_OP_REMOVAL = 0.116 d⁻¹` (= 2 × `K_BIODEG_WATER_ANCHOR`,
   t½ ≈ 6 d @ 20°C), Q10-scaled like `K_BIODEG_*`:
   `mass_op_dissolved → mass_op_degraded`. This represents further
   biodegradation/photodegradation/dilution of OPs in the water column and is
   independent of solar UV.
2. **Photooxidation (daytime only, `I_UV0 > 0`)**: of the mass `dM_phox`
   removed from oil-phase aromatics each step, a fraction
   `OP_DISSOLVED_FRACTION = 0.5` (50/50, first-pass estimate) routes to
   `mass_op_dissolved` instead of `mass_photooxidized`; the remainder stays
   in `mass_photooxidized` (oil-phase terminal sink) as before. This reflects
   OPs (oxy-PAHs, quinones, carboxylic acids) being markedly more
   water-soluble than their parent PAHs (Aeppli et al. 2012).
3. **Fractions**: `fraction_photooxidized`, `fraction_op_dissolved`, and
   `fraction_op_degraded` are all computed relative to the same conserved
   `initial_mass` denominator, which now also includes `mass_op_dissolved`
   and `mass_op_degraded`. The analogous `initial_mass` formula in
   `biodegradation_cicese()` (for `fraction_biodegraded_from_oil`) was updated
   the same way, for consistency.

### Verification

`test_cicoil_e2e.py` T3.5–T3.7 (mass conservation) and T4.1–T4.3 (dissolved
lifecycle) re-run 20/20 PASS — biodegradation/photooxidation are off by
default in those tests, so the new pools stay at zero and don't change the
totals. New targeted tests for droplet-size scaling and OP routing/removal
are in `test_biodegradation.py` and `test_photooxidation_op.py`.
