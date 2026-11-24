#!/bin/bash

# Directivas  (equivalente a cicoil.slurm de CICOIL v1, para Ixachi)
#SBATCH -J cicoilv2-job
#SBATCH -p oc_cpus
#SBATCH -o cicoilv2-job.out
#SBATCH -e cicoilv2-job.err
#SBATCH -n 1
#SBATCH --cpus-per-task=4
# Abre Jupyter Notebook en un nodo de cómputo. Copie este archivo a su directorio de trabajo y cambie el nombre del ambiente (línea conda activate).
# Envíe con:  sbatch cicoilv2.slurm   y use el comando ssh -N -L ... que aparece en cicoilv2-job.out

source /data/apps/python/miniforge2025/etc/profile.d/conda.sh
conda activate cicoilv2          # <- nombre de SU ambiente

# Configurar Jupyter
export XDG_RUNTIME_DIR="/tmp/jupyter-runtime-${USER}"
mkdir -p $XDG_RUNTIME_DIR

cd $SLURM_SUBMIT_DIR

JUPYTERPORT=9999
echo "ssh -N -L $JUPYTERPORT:$(hostname):$JUPYTERPORT $USER@logixachi.cicese.mx"
jupyter notebook --no-browser --ip=$(hostname) --port=$JUPYTERPORT --NotebookApp.token=''
