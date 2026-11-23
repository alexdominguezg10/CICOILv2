#!/bin/bash
#SBATCH -J analyze_nan
#SBATCH -p oc_cpus
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -t 00:10:00
#SBATCH -o /LUSTRE/adomingu/CICOIL_SHARED/logs/analyze_nan_%j.out
#SBATCH -e /LUSTRE/adomingu/CICOIL_SHARED/logs/analyze_nan_%j.err
# Read-only: which elements keep status 0 but lose their mass?   usage: sbatch --export=ALL,NC=/path/file.nc 13_analyze_nan.slurm
set -uo pipefail
SH=/LUSTRE/adomingu/CICOIL_SHARED
source "${CONDA_ROOT:-/data/apps/python/miniforge2025}/etc/profile.d/conda.sh"
conda activate ${ENV:-$SH/env/cicoilv2}
python $SH/readers_check/analyze_nan.py "${NC:?set NC=/path/file.nc}"
echo "exit $?"
