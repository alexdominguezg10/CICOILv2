# Calibrated rate constants for CICOILv2 (MAYA crude, Bay of Campeche).
#
# All 17 pseudo-components in order:
#   [1..8]  = aliphatic groups  ali_G1 (light, <C10) .. ali_G8 (heavy waxes C36+)
#   [9..16] = aromatic groups   arom_G1 (benzene-like) .. arom_G8 (pyrene+)
#   [17]    = residue           (asphaltenes + resins)
#
# NTuple is a fixed-size, stack-allocated tuple from Base.
# On GPU (future), these can be converted to SVector for register storage.

const N_COMPONENTS = 17
const N_ALI        = 8
const N_AROM       = 8

# ── Biodegradation ────────────────────────────────────────────────────────────
# Rate constants (d⁻¹ at T_REF=20°C) for DISPERSED/ENTRAINED oil droplets
# (z<0).  Half-lives from Prince et al. (2003, 2017) + Bacosa et al. (2015).
# Surface slicks (z==0) are suppressed by BIODEG_SURFACE_FACTOR.

const K_BIODEG_ALI = NTuple{8, Float32}((
    0.1386f0,   # G1  t½ ~5 d   (light, <C10 — most bioavailable)
    0.1155f0,   # G2  t½ ~6 d
    0.0990f0,   # G3  t½ ~7 d
    0.0866f0,   # G4  t½ ~8 d
    0.0693f0,   # G5  t½ ~10 d
    0.0462f0,   # G6  t½ ~15 d
    0.0277f0,   # G7  t½ ~25 d  (waxes)
    0.0173f0,   # G8  t½ ~40 d  (heavy waxes C36+)
))

const K_BIODEG_AROM = NTuple{8, Float32}((
    0.0866f0,   # G1  t½ ~8 d   benzene-like (BTEX)
    0.0693f0,   # G2  t½ ~10 d  toluene-like
    0.0578f0,   # G3  t½ ~12 d  xylenes
    0.0462f0,   # G4  t½ ~15 d  naphthalene (2-ring)
    0.0347f0,   # G5  t½ ~20 d  acenaphthene (2–3 ring)
    0.0277f0,   # G6  t½ ~25 d  phenanthrene (3-ring)
    0.0173f0,   # G7  t½ ~40 d  chrysene (4-ring)
    0.0116f0,   # G8  t½ ~60 d  pyrene+ (4–5 ring)
))

const K_BIODEG_RESIDUE = 0.00139f0  # t½ ~500 d — essentially non-biodegradable

# Full 17-component vector in canonical order [ali, arom, residue]
const K_BIODEG_ALL = NTuple{17, Float32}((
    K_BIODEG_ALI..., K_BIODEG_AROM..., K_BIODEG_RESIDUE
))

const BIODEG_T_REF            = 20.0f0   # °C reference temperature (Adcroft 2010)
const BIODEG_Q10              = 3.0f0    # Q10 temperature coefficient
const BIODEG_SURFACE_FACTOR   = 1.0f0 / 30.0f0   # suppression for surface slicks (z==0)
const BIODEG_DROPLET_D_REF    = 0.001f0  # m — reference droplet diameter
const BIODEG_DROPLET_EXPONENT = 1.0f0    # SA/V ~ 1/d scaling (King 1992)
const BIODEG_DROPLET_D_MIN    = 1.0f-5   # m (10 µm — clamp lower bound)
const BIODEG_DROPLET_D_MAX    = 1.0f-2   # m (1 cm — clamp upper bound)

# Dissolved-phase (is_dissolved particles, TAMOC far-field plume)
# Rescaled so mean rate across 16 HC groups = 0.0578 d⁻¹ (Adcroft 12-day t½)
const K_BIODEG_WATER_ANCHOR = 0.0578f0
const _k_biodeg_mean = sum(Float32[K_BIODEG_ALI..., K_BIODEG_AROM...]) / 16.0f0
const K_BIODEG_WATER = NTuple{17, Float32}(
    Tuple((K_BIODEG_WATER_ANCHOR / _k_biodeg_mean) .* Float32[K_BIODEG_ALL...])
)
const BIODEG_WATER_DEPLETION_FRAC = 0.01f0

# ── Photooxidation ────────────────────────────────────────────────────────────
# Rate constants (m² J⁻¹) for the 8 aromatic groups ONLY.
# Calibrated for floating crude oil film (Bacosa et al. 2015, GoM July UV).
# Full-series recalibration 2026-06-13: f*=0.3834 applied to G4/G6/G8 anchors;
# G5/G7 are log-linear interpolations between their bracketing recalibrated anchors.
# G1/G2/G3 from surface_fate_coupled.py (CICESE), unchanged.

const K_PHOTO_AROM = NTuple{8, Float32}((
    2.9f-8,    # G1  benzene       t½ ~14 d @ GoM July UV
    2.0f-8,    # G2  toluene       t½ ~20 d
    5.0f-8,    # G3  xylenes       t½ ~8 d
    4.260f-8,  # G4  naphthalene   t½ ~9.4 d  (recal f*=0.3834)
    5.145f-8,  # G5  acenaphthene  t½ ~7.8 d  (interpolated G4–G6)
    6.176f-8,  # G6  phenanthrene  t½ ~6.5 d  (recal f*=0.3834)
    1.581f-8,  # G7  chrysene      t½ ~25.3 d (interpolated G6–G8)
    5.110f-9,  # G8  pyrene+       t½ ~79.8 d (recal f*=0.3834)
))

const K_D_UV          = 0.10f0    # m⁻¹  Beer-Lambert UV attenuation (clear GoM)
const UV_DEPTH_CUTOFF = 50.0f0    # m    skip particles deeper than this
const OP_DISSOLVED_FRACTION = 0.5f0   # fraction of photoox mass → dissolved OP pool
const K_OP_REMOVAL    = 0.116f0   # d⁻¹  OP pool decay (t½ ~6 d @ 20°C, Q10-scaled)

# ── Emulsification ────────────────────────────────────────────────────────────
# Fingas & Fieldhouse (2004) viscosity-stability scaling + EMUL_RATE_SCALE.
# EMUL_RATE_SCALE=0.0142 calibrated so mean surface viscosity crosses
# NU_TARBALL=10,000 cSt at day ~17 (Norte 2 window, Abkatún 2026).

const EMUL_STABILITY_NU_REF = 1000.0f0   # cSt — reference viscosity for stability factor
const EMUL_STABILITY_MIN    = 0.3f0      # floor — very fluid fresh oil
const EMUL_STABILITY_MAX    = 3.0f0      # ceiling
const EMUL_RATE_SCALE       = 0.0142f0
const EMUL_MAX_WATER_DEFAULT = 0.9f0     # mousse saturation cap
const NU_TARBALL            = 10_000.0f0 # cSt — viscosity threshold for tar balls

# ── Evaporation — Cox vapour-pressure coefficients ───────────────────────────
# Cox coefficients are computed analytically from boiling temperature Tb (K):
#   A0 =  1.77838e-3 * Tb + 2.186529
#   A1 = -1.436577e-7 * Tb - 1.908404e-3
#   A2 = -2.734678e-9 * Tb + 2.858042e-6
# These are set in OilProps (computed from Tb loaded from pseudo_chemdata.csv).

const R_GAS   = 8.31446f0    # J mol⁻¹ K⁻¹
const P_ATM   = 101325.0f0   # Pa
const MAX_AGE_EVAP_DAYS = 120.0f0  # particles older than this skip evaporation
