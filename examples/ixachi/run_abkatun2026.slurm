#!/bin/bash
#SBATCH -J abk2026
#SBATCH -p oc_cpus
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o %x_%j.out
# Abkatun 2026 example, 2D (surface) or 3D (subsea seed + TAMOC droplet sizes).  Submit from YOUR own directory:
#   mkdir -p /LUSTRE/$USER/CICOIL_practica && cd /LUSTRE/$USER/CICOIL_practica
#   sbatch --export=ALL,MODE=2d /LUSTRE/adomingu/CICOIL_SHARED/ejemplo/run_abkatun2026.slurm
#   sbatch --export=ALL,MODE=3d /LUSTRE/adomingu/CICOIL_SHARED/ejemplo/run_abkatun2026.slurm
# Optional (comma-free values only):  N=2000  DAYS=10  SEED=1  DMODE=legacy|transport|budget  DROPLETS=tamoc|default  DEPTH=19  FORCING=coupled|glorys (default coupled; 3d tamoc needs coupled)  OUTDT=6 (hours)  LOGLEVEL=30 (20 = info)
# Smoke test:  --export=ALL,MODE=3d,N=200,DAYS=2
set -euo pipefail
SH=/LUSTRE/adomingu/CICOIL_SHARED
MODE=${MODE:?set MODE=2d or MODE=3d with --export=ALL,MODE=2d}
source "${CONDA_ROOT:-/data/apps/python/miniforge2025}/etc/profile.d/conda.sh"
conda activate ${ENV:-$SH/env/cicoilv2}
echo "user $USER host $(hostname) mode $MODE start $(date)"
python $SH/ejemplo/abkatun2026.py --mode $MODE --forcing ${FORCING:-coupled} --forcing-dir $SH/examples/abkatun2026/inputs --out ./output_abkatun2026 \
  -n ${N:-2000} --days ${DAYS:-10} --seed ${SEED:-1} --dispersion-mode ${DMODE:-legacy} --droplets ${DROPLETS:-tamoc} --depth ${DEPTH:-19} \
  --out-dt ${OUTDT:-6} --loglevel ${LOGLEVEL:-30}
echo "end $(date)"
