#!/bin/bash
#SBATCH -J cicoil_readers
#SBATCH -p oc_cpus
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH -o /LUSTRE/adomingu/CICOIL_SHARED/logs/check_readers_%j.out
#SBATCH -e /LUSTRE/adomingu/CICOIL_SHARED/logs/check_readers_%j.err
# Reader smoke test in the shared env (CROCO native, CROCO/WRF/WW3 via their CF files, GLORYS/ERA5/ANFC, NEMO).
# Reads big files lazily but this is netCDF work: compute node only.  Exit code 1 if any reader FAILS; SKIP is not a pass.
# Needs 00_build_shared_env.slurm done.  Optional: NEMO_FILE=/path/to/nemo.nc  (otherwise the NEMO check is reported SKIP)
#   sbatch /LUSTRE/adomingu/CICOIL_SHARED/hpc/10_check_readers.slurm
#   sbatch --export=ALL,NEMO_FILE=/path/nemo.nc /LUSTRE/adomingu/CICOIL_SHARED/hpc/10_check_readers.slurm
set -uo pipefail
SH=/LUSTRE/adomingu/CICOIL_SHARED
ENV=${ENV:-$SH/env/cicoilv2}
source "${CONDA_ROOT:-/data/apps/python/miniforge2025}/etc/profile.d/conda.sh"
conda activate "$ENV" || { echo "FAIL: cannot activate $ENV (run 00_build_shared_env.slurm first)"; exit 2; }
export ABK_DIR=${ABK_DIR:-$SH/examples/abkatun2026/inputs}
[ -d "$ABK_DIR" ] || export ABK_DIR=/LUSTRE/adomingu/Abkatun_ixachi/forcing
export COUPLED_DIR=${COUPLED_DIR:-/LUSTRE/adomingu/CICOIL_COMPARISON/forcing/S1_feb2026}
export CROCO_NATIVE=${CROCO_NATIVE:-/LUSTRE/adomingu/TRIPLE_COUPLED_GOM/FORECAST/runs/ir/coupled_forecast_2026-10-07_2026-10-12_forecastB_r2/croco_avg.nc}
export NEMO_FILE=${NEMO_FILE:-}
echo "user $USER host $(hostname) env $ENV start $(date)"
echo "ABK_DIR=$ABK_DIR"; echo "COUPLED_DIR=$COUPLED_DIR"; echo "CROCO_NATIVE=$CROCO_NATIVE"; echo "NEMO_FILE=${NEMO_FILE:-(not set)}"
python $SH/readers_check/check_readers.py --json $SH/logs/check_readers_${SLURM_JOB_ID:-manual}.json
rc=$?
echo "check_readers exit $rc  end $(date)"
exit $rc
