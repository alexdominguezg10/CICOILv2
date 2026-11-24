#!/bin/bash
# Prueba de instalación desde cero en una laptop o cualquier máquina Linux (sin Slurm): pasos 2, 4 y 5 de las instrucciones + prueba de edición.
# Antes:  mkdir -p $DEST && cd $DEST && git clone https://github.com/alexdominguezg10/CICOILv2.git
#         bash CICOILv2/MODEL_CICOILv2/instalar_cicoilv2.sh preparar
# Luego:  DEST=$DEST ENVPREFIX=$DEST/env_cicoilv2 bash $DEST/CICOILv2/MODEL_CICOILv2/test_instalacion.sh
# (ENVNAME=cicoilv2 en lugar de ENVPREFIX para un ambiente con nombre.)  Necesita conda (miniconda/miniforge), gcc y gfortran.
set -eo pipefail
DEST=${DEST:?defina DEST}
REPO=$DEST/CICOILv2
CONDA_BASE=${CONDA_BASE:-$(conda info --base)}
source "$CONDA_BASE/etc/profile.d/conda.sh"
echo "host $(hostname) dest $DEST start $(date)"
[ -d $REPO/MODEL_CICOILv2 ] && [ -d $DEST/tamoc/.git ] || { echo "FAIL: faltan $REPO o $DEST/tamoc (corra 'preparar')"; exit 2; }
echo "== 2. crear el ambiente"
if [ -n "${ENVPREFIX:-}" ]; then
  [ -e "$ENVPREFIX" ] && { echo "FAIL: $ENVPREFIX ya existe"; exit 4; }
  conda env create -y -f $REPO/MODEL_CICOILv2/entorno_cicoilv2.yml -p "$ENVPREFIX"; conda activate "$ENVPREFIX"
else
  conda env create -y -f $REPO/MODEL_CICOILv2/entorno_cicoilv2.yml -n ${ENVNAME:-cicoilv2}; conda activate ${ENVNAME:-cicoilv2}
fi
python --version; python -c "import numpy; print('numpy', numpy.__version__)"; which gfortran gcc || echo "WARN: falta gfortran o gcc"
echo "== 4. instalar"
TAMOC_DIR=$DEST/tamoc bash $REPO/MODEL_CICOILv2/instalar_cicoilv2.sh instalar
echo "== 5. pruebas"
cd $REPO
for t in test_deactivate_masks test_oil_density_ids test_vertical_profile; do
  echo "--- $t"; python $t.py 2>&1 | grep -v "^_____\|^setting oil_type\|^Selected Oil" | tail -4
done
echo "== prueba de edición"
cp $REPO/ciceseoil.py $REPO/ciceseoil.py.bak_test
printf '\nCICOIL_EDIT_TEST = 42\n' >> $REPO/ciceseoil.py      # ciceseoil.py no termina en salto de línea
python -c "import opendrift.models.openoil.ciceseoil as m; assert m.CICOIL_EDIT_TEST == 42; print('edición de ciceseoil.py usada: OK')"
mv $REPO/ciceseoil.py.bak_test $REPO/ciceseoil.py        # deshacer la prueba
echo "DONE $(date)"
