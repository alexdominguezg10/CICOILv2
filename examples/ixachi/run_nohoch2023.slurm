#!/bin/bash
#SBATCH -J nohoch2023
#SBATCH -p oc_cpus
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o %x_%j.out
# Nohoch 2023 example (S0: 6 Jul 2023 polygon release from Jiao et al. 2025, Maya, Mercator GLORYS12 + ERA5), CICOILv2 b032140.
# Runs the validated S0 adapter from the COPY staged by 20_stage_examples.slurm; the shared forcing is only read:
# a private data dir is built here (symlinks to mercator/ and era5/, real dirs for era5_cf/ and opendrift/, which the script writes).
# Submit from YOUR own directory:
#   mkdir -p /LUSTRE/$USER/CICOIL_practica && cd /LUSTRE/$USER/CICOIL_practica
#   sbatch /LUSTRE/adomingu/CICOIL_SHARED/ejemplo/run_nohoch2023.slurm
# Optional: N=300 DAYS=3 SEED=1 DMODE=legacy|transport|budget   (defaults are a smoke test; the full case is N=3000 DAYS=30)
set -euo pipefail
SH=/LUSTRE/adomingu/CICOIL_SHARED
EX=$SH/examples/nohoch2023
source "${CONDA_ROOT:-/data/apps/python/miniforge2025}/etc/profile.d/conda.sh"
conda activate ${ENV:-$SH/env/cicoilv2}
unset DIFF_K_PATCH
export CICOIL_COMMIT=b032140
export NOHOCH_DATA=$PWD/nohoch_data
mkdir -p $NOHOCH_DATA/era5_cf $NOHOCH_DATA/opendrift
ln -sfn $EX/forcing/mercator $NOHOCH_DATA/mercator
ln -sfn $EX/forcing/era5 $NOHOCH_DATA/era5
cp -n $EX/forcing/jiao2025_contours.geojson $NOHOCH_DATA/ 2>/dev/null || true
echo "user $USER host $(hostname) start $(date)  NOHOCH_DATA=$NOHOCH_DATA"
python $EX/adapters/run_nohoch_v2fix.py --dispersion-mode ${DMODE:-legacy} --seed ${SEED:-1} --n ${N:-300} --days ${DAYS:-3} \
  --out $PWD/output_nohoch2023/${DMODE:-legacy}_s${SEED:-1}
echo "end $(date)"
