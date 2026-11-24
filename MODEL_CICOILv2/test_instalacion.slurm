#!/bin/bash
#SBATCH -J cicoilv2_install
#SBATCH -p oc_cpus
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -t 01:30:00
#SBATCH -o cicoilv2_install_%j.out
# Prueba de instalación desde cero (pasos 2, 4 y 5 de INSTRUCCIONES_INSTALACION_CICOILv2.txt) más una prueba de edición de ciceseoil.py.
# Antes, en el login node (necesita git y red):
#   mkdir -p $DEST && cd $DEST && git clone https://github.com/alexdominguezg10/CICOILv2.git   (o copie el repositorio a $DEST/CICOILv2)
#   bash $DEST/CICOILv2/MODEL_CICOILv2/instalar_cicoilv2.sh preparar
# Luego:  sbatch --export=ALL,DEST=$DEST,ENVPREFIX=$DEST/env_cicoilv2 $DEST/CICOILv2/MODEL_CICOILv2/test_instalacion.slurm
# Para una instalación normal use  ENVNAME=cicoilv2  en lugar de ENVPREFIX.
set -eo pipefail
DEST=${DEST:?defina DEST con --export=ALL,DEST=...}
REPO=$DEST/CICOILv2
source "${CONDA_ROOT:-/data/apps/python/miniforge2025}/etc/profile.d/conda.sh"
echo "user $USER host $(hostname) dest $DEST start $(date)"
[ -d $REPO/MODEL_CICOILv2 ] && [ -d $DEST/tamoc/.git ] || { echo "FAIL: faltan $REPO o $DEST/tamoc (corra 'preparar' en el login node)"; exit 2; }

echo "== 2. crear el ambiente"
if [ -n "${ENVPREFIX:-}" ]; then
  [ -e "$ENVPREFIX" ] && { echo "FAIL: $ENVPREFIX ya existe"; exit 4; }
  conda env create -y -f $REPO/MODEL_CICOILv2/entorno_cicoilv2.yml -p "$ENVPREFIX"; conda activate "$ENVPREFIX"
else
  conda env create -y -f $REPO/MODEL_CICOILv2/entorno_cicoilv2.yml -n ${ENVNAME:-cicoilv2}; conda activate ${ENVNAME:-cicoilv2}
fi
python --version; python -c "import numpy; print('numpy', numpy.__version__)"; which gfortran || echo "WARN: gfortran no encontrado"

echo "== 4. instalar"
TAMOC_DIR=$DEST/tamoc bash $REPO/MODEL_CICOILv2/instalar_cicoilv2.sh instalar

echo "== 5. pruebas"
cd $REPO
for t in test_deactivate_masks test_oil_density_ids test_vertical_profile; do
  echo "--- $t"; python $t.py 2>&1 | grep -v "^_____\|^setting oil_type\|^Selected Oil" | tail -4
done

echo "== prueba de edición: una línea agregada a ciceseoil.py debe verse al importar"
cp $REPO/ciceseoil.py $REPO/ciceseoil.py.bak_test
printf '\nCICOIL_EDIT_TEST = 42\n' >> $REPO/ciceseoil.py      # ciceseoil.py no termina en salto de línea
python -c "import opendrift.models.openoil.ciceseoil as m; assert m.CICOIL_EDIT_TEST == 42; print('edición de ciceseoil.py usada: OK')"
mv $REPO/ciceseoil.py.bak_test $REPO/ciceseoil.py        # deshacer la prueba
echo "DONE $(date)"
