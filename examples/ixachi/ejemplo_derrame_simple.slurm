#!/bin/bash
#SBATCH -J cicoil_ejemplo
#SBATCH -p oc_cpus
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o %x_%j.out
# Ejemplo minimo de CICOILv2 (ver ejemplo_derrame_simple.py).  Se envia desde TU directorio de trabajo:
#   mkdir -p /LUSTRE/$USER/CICOIL_ejemplo && cd /LUSTRE/$USER/CICOIL_ejemplo
#   cp /LUSTRE/adomingu/CICOIL_SHARED/ejemplo/ejemplo_derrame_simple.* .
#   sbatch ejemplo_derrame_simple.slurm
# Para otro modo de dispersion:  sbatch --export=ALL,MODE=budget ejemplo_derrame_simple.slurm
set -euo pipefail
ENV=${ENV:-/LUSTRE/adomingu/CICOIL_SHARED/env/cicoilv2}
source "${CONDA_ROOT:-/data/apps/python/miniforge2025}/etc/profile.d/conda.sh"
conda activate "$ENV"
echo "user $USER host $(hostname) start $(date)"
python ejemplo_derrame_simple.py --out salida_ejemplo_${MODE:-legacy}.nc --dispersion-mode ${MODE:-legacy}
echo "end $(date)"
