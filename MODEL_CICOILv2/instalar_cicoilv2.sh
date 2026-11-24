#!/bin/bash
# Instalación de CICOILv2 para el grupo de modelación de derrames (esquema de CICOIL v1: ambiente conda + TAMOC editable +
# archivos de CICOIL en su propia copia, que usted puede modificar).
#
#   bash MODEL_CICOILv2/instalar_cicoilv2.sh preparar    # necesita git y red (login node de Ixachi: solo clonar, nada de cómputo)
#   bash MODEL_CICOILv2/instalar_cicoilv2.sh instalar    # con el ambiente conda ACTIVO; en un nodo de cómputo (compila TAMOC)
#
# preparar: clona TAMOC en ../tamoc (junto a este repositorio), lo deja en el commit f689dc7 (= V4.1.2 de CHANGES.txt) y aplica el parche
#           tamoc_patches/tamoc_f689dc7_cicoil.patch (numpy 2.x en f2py, un escalar en single_bubble_model.py, compuestos PAH en ChemData.csv).
# instalar: pip install -e de TAMOC y enlaza (ln -s) los archivos de CICOIL de ESTE repositorio dentro del OpenDrift 1.14.9 instalado,
#           de modo que editar ciceseoil.py (etc.) aquí se usa de inmediato y aparece en "git diff".
# Variables: TAMOC_DIR (por defecto <repositorio>/../tamoc), TAMOC_URL.
set -euo pipefail
MODE=${1:-}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/.." && pwd)
TAMOC_DIR=${TAMOC_DIR:-$(dirname "$REPO")/tamoc}
TAMOC_URL=${TAMOC_URL:-https://github.com/socolofs/tamoc.git}
TAMOC_COMMIT=$(cat "$HERE/tamoc_patches/TAMOC_UPSTREAM_COMMIT.txt")
PATCH="$HERE/tamoc_patches/tamoc_${TAMOC_COMMIT}_cicoil.patch"

case "$MODE" in
preparar)
  command -v git >/dev/null || { echo "FAIL: no hay git en este nodo (use el login node)"; exit 2; }
  [ -e "$TAMOC_DIR" ] && { echo "FAIL: $TAMOC_DIR ya existe (bórrelo o use otro TAMOC_DIR)"; exit 4; }
  git clone -q "$TAMOC_URL" "$TAMOC_DIR"
  cd "$TAMOC_DIR"
  git checkout -q -b cicoil-4.1.2 "$TAMOC_COMMIT"
  git apply --check "$PATCH"
  git apply --whitespace=nowarn "$PATCH"
  git config user.name >/dev/null 2>&1 || git config user.name "CICOIL user"
  git config user.email >/dev/null 2>&1 || git config user.email "cicoil@localhost"
  git add -A && git commit -q -m "TAMOC $TAMOC_COMMIT + parche CICOIL (numpy 2.x, PAH en ChemData.csv)"
  git log --oneline -2
  echo "TAMOC listo en $TAMOC_DIR (rama cicoil-4.1.2). Siguiente paso: activar el ambiente y correr  instalar"
  ;;
instalar)
  [ -n "${CONDA_PREFIX:-}" ] || { echo "FAIL: active primero el ambiente conda (conda activate cicoilv2)"; exit 2; }
  [ -d "$TAMOC_DIR/.git" ] || { echo "FAIL: falta $TAMOC_DIR; corra primero:  bash $0 preparar"; exit 2; }
  python - <<'PY'
import opendrift
assert opendrift.__version__ == '1.14.9', f'se necesita OpenDrift 1.14.9, hay {opendrift.__version__}'
print('OpenDrift', opendrift.__version__)
PY
  echo "== TAMOC editable"
  which gfortran >/dev/null || echo "WARN: no se encuentra gfortran (TAMOC lo necesita; conda install -c conda-forge gfortran)"
  ( cd "$TAMOC_DIR" && pip install --no-build-isolation --no-deps -e . )
  echo "== enlazar los archivos de CICOILv2 dentro de OpenDrift"
  OD=$(python -c "import opendrift,os;print(os.path.dirname(opendrift.__file__))")
  link() { local src="$1" dst="$2"; [ -e "$src" ] || { echo "FAIL: falta $src"; exit 3; }; ln -sfn "$src" "$dst"; }
  mkdir -p "$OD/models/openoil/data"
  for f in ciceseoil slick_field fluid_properties cicoil_estimations tamoc_plume run_tamoc tamoc_chemical_properties vertical_profile; do
    link "$REPO/$f.py" "$OD/models/openoil/$f.py"
  done
  for f in "$REPO"/data/*.csv; do link "$f" "$OD/models/openoil/data/$(basename "$f")"; done
  for f in "$REPO"/readers/*.py; do link "$f" "$OD/readers/$(basename "$f")"; done
  link "$REPO/export/io_stat_nc.py" "$OD/export/io_stat_nc.py"
  echo "== verificar"
  REPO="$REPO" TAMOC_DIR="$TAMOC_DIR" python - <<'PY'
import os, opendrift, tamoc
import opendrift.models.openoil.ciceseoil as m
import opendrift.models.openoil.vertical_profile as v
repo = os.path.realpath(os.environ['REPO']); tdir = os.path.realpath(os.environ['TAMOC_DIR'])
for name, mod, base in (('ciceseoil', m, repo), ('vertical_profile', v, repo), ('tamoc', tamoc, tdir)):
    p = os.path.realpath(mod.__file__); print(f'{name:17s} {p}'); assert p.startswith(base), f'{name} no viene de {base}'
import tamoc.dbm_f
print('opendrift', opendrift.__version__, '| tamoc', tamoc.__version__, '| dbm_f (Fortran) OK')
PY
  echo "LISTO. Probar:  cd $REPO && python test_deactivate_masks.py && python test_oil_density_ids.py && python test_vertical_profile.py"
  ;;
*)
  echo "uso: bash $0 preparar | instalar"; exit 1 ;;
esac
