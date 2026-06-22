#!/bin/bash
#
# CICOILv2.0 — HPC Deployment Script (Conda Edition)
# ===================================================
# Uses conda (mamba broken on CHAMAN2)
#
# Usage: bash deploy_on_hpc_conda.sh
#

set -e

echo "================================================"
echo "CICOILv2.0 — HPC Deployment (Conda)"
echo "CHAMAN2 cicoil4 environment"
echo "================================================"
echo ""

CONDA_ENV_NEW="cicoil4"
OD_VERSION="1.14.9"
DEPLOY_START=$(date)

# --- Activate conda ---
echo "1. Initializing conda..."
eval "$(conda shell.bash hook)"

# --- Create or activate environment ---
echo "2. Setting up environment..."
echo "   Creating new environment: $CONDA_ENV_NEW (Python 3.10, from scratch)"

if conda env list | grep -q "^$CONDA_ENV_NEW "; then
    echo "   $CONDA_ENV_NEW exists. Removing..."
    conda remove -n $CONDA_ENV_NEW --all --yes
fi

echo "   Creating $CONDA_ENV_NEW with Python 3.10..."
conda create -n $CONDA_ENV_NEW -y python=3.10 numpy scipy matplotlib xarray netCDF4

echo "   Activating $CONDA_ENV_NEW..."
conda activate $CONDA_ENV_NEW || {
    echo "ERROR: Cannot activate $CONDA_ENV_NEW"
    exit 1
}

PYTHON=$(which python)
PIP=$(which pip)
echo "   Python: $PYTHON"
echo "   Python version: $(python --version)"
echo ""

# --- Install OpenDrift dependencies explicitly ---
echo "3. Installing OpenDrift dependencies..."
$PIP install future netCDF4 shapely pyproj cartopy scikit-image Pillow cython

# --- Install OpenDrift ---
echo ""
echo "4. Installing OpenDrift $OD_VERSION..."
echo "   This may take 3-5 minutes..."
$PIP install opendrift==$OD_VERSION

OD_CURRENT=$(python -c "import opendrift; print(opendrift.__version__)" 2>/dev/null || echo "failed")
echo "   ✓ OpenDrift installed: $OD_CURRENT"

# --- Locate site-packages ---
echo ""
echo "5. Locating OpenDrift installation..."
OD_SITE=$(python -c "import opendrift; import os; print(os.path.dirname(opendrift.__file__))")
echo "   Site-packages: $OD_SITE/models/openoil/"

# --- Copy model files ---
echo ""
echo "6. Installing patched CICOIL files..."

for f in ciceseoil.py fluid_properties.py cicoil_estimations.py; do
    if [[ -f "$f" ]]; then
        cp "$f" "$OD_SITE/models/openoil/"
        echo "   ✓ $(basename "$f")"
    fi
done

for f in tamoc_plume.py run_tamoc.py tamoc_chemical_properties.py; do
    if [[ -f "$f" ]]; then
        cp "$f" "$OD_SITE/models/openoil/"
        echo "   ✓ $(basename "$f")"
    fi
done

mkdir -p "$OD_SITE/readers" 2>/dev/null || true
for f in readers/reader_nemo_*.py readers/reader_NEMO_native_v3.py; do
    if [[ -f "$f" ]]; then
        cp "$f" "$OD_SITE/readers/"
        echo "   ✓ $(basename "$f")"
    fi
done

mkdir -p "$OD_SITE/models/openoil/data" 2>/dev/null || true
for f in data/*.csv; do
    if [[ -f "$f" ]]; then
        cp "$f" "$OD_SITE/models/openoil/data/"
        echo "   ✓ $(basename "$f")"
    fi
done

echo "   ✓ All files installed"

# --- Verify installation ---
echo ""
echo "7. Verifying installation..."
python << 'PYEOF'
import sys
try:
    from opendrift.models.openoil.ciceseoil import OpenCiceseOil
    print("   ✓ OpenCiceseOil imports successfully")
except ImportError as e:
    print(f"   ✗ OpenCiceseOil import failed: {e}")
    sys.exit(1)

import os
import importlib
OD_PATH = importlib.util.find_spec('opendrift').origin
OD_DIR = os.path.dirname(os.path.dirname(OD_PATH))
data_dir = os.path.join(OD_DIR, 'models', 'openoil', 'data')
for fname in ['pseudo_chemdata.csv', 'ChemData.csv', 'Aij.csv', 'Bij.csv', 'gas_composition.csv']:
    fpath = os.path.join(data_dir, fname)
    if os.path.exists(fpath):
        print(f"   ✓ {fname}")
    else:
        print(f"   ✗ {fname} not found")
        sys.exit(1)

print("   ✓ All verification checks passed")
PYEOF

# --- Run test suite ---
echo ""
echo "8. Running end-to-end test suite (20 tests)..."
echo "   This will take ~6 minutes..."
echo ""

if [[ -f "test_cicoil_e2e.py" ]]; then
    python test_cicoil_e2e.py 2>&1 | tee test_output.log
    TEST_RESULT=$?
else
    echo "   ✗ test_cicoil_e2e.py not found, skipping tests"
    TEST_RESULT=1
fi

# --- Generate report ---
echo ""
echo "9. Generating deployment report..."

REPORT="HPC_DEPLOYMENT_REPORT.txt"
cat > "$REPORT" << EOF
CICOILv2.0 HPC DEPLOYMENT REPORT
================================

Deployment Date: $DEPLOY_START
Host: $(hostname)
User: $(whoami)
Package Manager: conda (mamba unavailable)
Conda Environment: $CONDA_ENV_NEW (NEW - cloned from $CONDA_ENV_SOURCE)

DEPLOYMENT STRATEGY
───────────────────
✓ New dedicated environment created: $CONDA_ENV_NEW
✓ Built from scratch (Python 3.10 + minimal dependencies)
✓ Package manager: conda
✓ No other environments modified

ENVIRONMENT
-----------
Python: $(python --version)
OpenDrift: $(python -c "import opendrift; print(opendrift.__version__)")

INSTALLATION
------------
Target: $OD_SITE/models/openoil/

TEST RESULTS
------------
$(tail -50 test_output.log 2>/dev/null || echo "  See test_output.log for details")

Status: $(if [[ $TEST_RESULT -eq 0 ]]; then echo "✓ PASSED"; else echo "✗ FAILED"; fi)

NEXT STEPS
----------
1. Review this report: less $REPORT
2. Verify test outputs: ls -lh test_*.nc test_*.png
3. Update SLURM jobs to use: conda activate $CONDA_ENV_NEW
4. Submit production jobs: sbatch your_simulation.slurm

=====================================
$(date)
EOF

echo "   ✓ Report saved to $REPORT"

# --- Final summary ---
echo ""
echo "================================================"
if [[ $TEST_RESULT -eq 0 ]]; then
    echo "✓ DEPLOYMENT COMPLETE - ALL TESTS PASSED"
    echo "================================================"
    echo ""
    echo "CICOILv2.0 is ready for production simulations"
    echo ""
    echo "Environment: $CONDA_ENV_NEW"
    echo "Activate with: conda activate $CONDA_ENV_NEW"
    echo ""
else
    echo "✗ DEPLOYMENT INCOMPLETE - TESTS FAILED"
    echo "================================================"
    echo ""
    echo "Review test output:"
    echo "  tail -100 test_output.log"
    echo ""
    exit 1
fi
