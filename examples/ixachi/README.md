# CICOILv2 examples for Ixachi (CICESE)

Scripts used for the CICOILv2 hand-over to the oil-spill group (2026-11-23). They are written for the shared
installation on Ixachi under `/LUSTRE/adomingu/CICOIL_SHARED` (environment, preprocessed forcing); the paths are at the top of
each `.slurm` file. All of them go through `sbatch` (partition `oc_cpus`): no compute on the login node.

| File | What it does |
|---|---|
| `abkatun2026.py`, `run_abkatun2026.slurm` | Abkatun 2026. `MODE=2d`: surface release. `MODE=3d`: subsea release at 19 m with TAMOC droplet sizes (v1 route: `vertical_profile.create_profile` -> `run_tamoc` -> `Plume(tracked=False)` -> `seed_plume_elements`). Prints the oil budget as % of the seeded oil; gas bubbles are reported separately (they carry `mass_gas`, `mass_oil = 0`). The 3D TAMOC case needs forcing with vertical levels (coupled CROCO-WRF-WW3 files; the GLORYS/ANFC files have one level). |
| `run_nohoch2023.slurm` | Nohoch 2023 (S0 polygon release, Maya, GLORYS12 + ERA5) through the S0 adapter staged in the shared directory. |
| `ejemplo_derrame_simple.{py,slurm}` | Minimal 24 h example with constant forcing (not a realistic case). |
| `readers_check/check_readers.py`, `hpc/10_check_readers.slurm` | Reads one file per reader (CROCO native, WW3, WRF via CF files, NEMO) and reports PASS / FAIL / SKIP. |
| `readers_check/analyze_budget.py`, `analyze_nan.py`, `hpc/13_*.slurm`, `hpc/14_*.slurm` | Read-only analysis of an output file: where the oil mass goes; which elements keep status 0 with NaN mass. |
| `readers_check/diag_croco_native.py`, `inspect_forcing.py` | Diagnostics of the native CROCO reader and of a forcing file. |

Requires CICOILv2 >= v2.1.0 (gas-bubble mask fix, `oil_density` ID fix, `vertical_profile.py`). Verified on Ixachi 2026-11-23:
Abkatun 2026 3D with v2.1.0 closes the oil budget at 100.00 % (job 22831); reader checks 9 PASS, NEMO skipped (no file).
`abkatun2026.py` still applies a run-time gas-bubble mask patch, which is harmless on v2.1.0 and keeps it working on b032140.
