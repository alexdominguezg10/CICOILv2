#!/bin/bash
#
# CICOILv2 + CICOILPhysics.jl — Ixachi Deployment Script
# =======================================================
# Deploys CICOILv2 with Julia/CUDA weathering backend on logixachi.cicese.mx
#
# Usage:
#   1. From your Mac:
#      scp -r CICOILPhysics.jl adomingu@logixachi.cicese.mx:/LUSTRE/adomingu/CICOILv2/
#      ssh adomingu@logixachi.cicese.mx
#
#   2. On Ixachi:
#      cd /LUSTRE/adomingu/CICOILv2
#      bash CICOILPhysics.jl/deploy_ixachi.sh
#
# Prerequisites:
#   - Julia already installed (available via PATH or module)
#   - Python via: module load python/miniforge
#   - NVIDIA H200 GPU with CUDA drivers
#

set -e

echo "================================================"
echo "CICOILv2 + CICOILPhysics.jl — Ixachi Deployment"
echo "================================================"
echo ""

# ── Configuration ────────────────────────────────────────────────────────────
DEPLOY_DIR="/LUSTRE/adomingu/CICOILv2"
JULIA_PROJECT="$DEPLOY_DIR/CICOILPhysics.jl"
CONDA_ENV="cicoil_julia"
OD_VERSION="1.14.9"

echo "Deploy dir:     $DEPLOY_DIR"
echo "Julia project:  $JULIA_PROJECT"
echo "Conda env:      $CONDA_ENV"
echo ""

# ── Step 1: Check Julia ──────────────────────────────────────────────────────
echo "1. Checking Julia..."
if ! command -v julia &>/dev/null; then
    echo "   Julia not found in PATH. Trying module load..."
    module load julia 2>/dev/null || true
fi

if ! command -v julia &>/dev/null; then
    echo "   ERROR: Julia not found. Install Julia >= 1.9 or load the module."
    exit 1
fi

JULIA_VER=$(julia --version)
echo "   ✓ $JULIA_VER"

# ── Step 2: Check GPU ────────────────────────────────────────────────────────
echo ""
echo "2. Checking GPU..."
if command -v nvidia-smi &>/dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
    GPU_MEM=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader 2>/dev/null | head -1)
    echo "   ✓ GPU: $GPU_NAME ($GPU_MEM)"
else
    echo "   ⚠ nvidia-smi not found — GPU tests will be skipped"
fi

# ── Step 3: Set up Julia environment ─────────────────────────────────────────
echo ""
echo "3. Setting up Julia environment for CICOILPhysics.jl..."

if [[ ! -f "$JULIA_PROJECT/Project.toml" ]]; then
    echo "   ERROR: $JULIA_PROJECT/Project.toml not found"
    echo "   Make sure CICOILPhysics.jl/ is under $DEPLOY_DIR"
    exit 1
fi

echo "   Installing Julia dependencies (Zarr + CUDA)..."
julia --project="$JULIA_PROJECT" -e '
    using Pkg
    Pkg.instantiate()
    # Add CUDA if not already present
    if !haskey(Pkg.project().dependencies, "CUDA")
        Pkg.add("CUDA")
    end
    Pkg.precompile()
    println("   ✓ Julia packages installed")
' 2>&1 | grep -E "✓|Updating|Installed|Precompiling|ERROR"

# ── Step 4: Test CICOILPhysics.jl (Julia) ────────────────────────────────────
echo ""
echo "4. Running CICOILPhysics.jl tests..."

julia --project="$JULIA_PROJECT" -e '
    push!(LOAD_PATH, "'$JULIA_PROJECT'")
    using CICOILPhysics
    println("   Backend: ", backend_name(backend()))
    println("   Has GPU: ", has_gpu(backend()))
' 2>&1

cd "$JULIA_PROJECT"
julia --project="$JULIA_PROJECT" test/runtests.jl 2>&1
JULIA_TEST_RESULT=$?
cd "$DEPLOY_DIR"

if [[ $JULIA_TEST_RESULT -eq 0 ]]; then
    echo "   ✓ Julia tests PASSED"
else
    echo "   ✗ Julia tests FAILED"
    echo "   Check output above for errors"
    exit 1
fi

# ── Step 5: Set up Python environment ─────────────────────────────────────────
echo ""
echo "5. Setting up Python environment..."

# Load miniforge
module load python/miniforge 2>/dev/null || true
eval "$(conda shell.bash hook)" 2>/dev/null || true

if conda env list 2>/dev/null | grep -q "^$CONDA_ENV "; then
    echo "   Environment $CONDA_ENV exists, activating..."
    conda activate $CONDA_ENV
else
    echo "   Creating new environment: $CONDA_ENV..."
    conda create -n $CONDA_ENV -y python=3.10 numpy scipy matplotlib xarray netCDF4
    conda activate $CONDA_ENV

    echo "   Installing OpenDrift $OD_VERSION + juliacall..."
    pip install opendrift==$OD_VERSION
    pip install juliacall
fi

echo "   Python: $(python --version)"
echo "   OpenDrift: $(python -c 'import opendrift; print(opendrift.__version__)' 2>/dev/null || echo 'not installed')"

# ── Step 6: Install CICOILv2 into OpenDrift ───────────────────────────────────
echo ""
echo "6. Installing CICOILv2 files into OpenDrift site-packages..."

OD_SITE=$(python -c "import opendrift; import os; print(os.path.dirname(opendrift.__file__))")
OD_OPENOIL="$OD_SITE/models/openoil"

# Core model files
for f in ciceseoil.py fluid_properties.py cicoil_estimations.py \
         tamoc_plume.py run_tamoc.py tamoc_chemical_properties.py; do
    if [[ -f "$DEPLOY_DIR/$f" ]]; then
        cp "$DEPLOY_DIR/$f" "$OD_OPENOIL/"
        echo "   ✓ $f"
    fi
done

# Data files
mkdir -p "$OD_OPENOIL/data" 2>/dev/null || true
for f in "$DEPLOY_DIR"/data/*.csv; do
    if [[ -f "$f" ]]; then
        cp "$f" "$OD_OPENOIL/data/"
        echo "   ✓ data/$(basename "$f")"
    fi
done

# Readers
mkdir -p "$OD_SITE/readers" 2>/dev/null || true
for f in "$DEPLOY_DIR"/readers/reader_nemo_*.py "$DEPLOY_DIR"/readers/reader_NEMO_*.py; do
    if [[ -f "$f" ]]; then
        cp "$f" "$OD_SITE/readers/"
        echo "   ✓ readers/$(basename "$f")"
    fi
done

# ── Step 7: Verify Python import ──────────────────────────────────────────────
echo ""
echo "7. Verifying Python installation..."

python << PYEOF
import sys
try:
    from opendrift.models.openoil.ciceseoil import OpenCiceseOil
    print("   ✓ OpenCiceseOil imports successfully")

    o = OpenCiceseOil.__new__(OpenCiceseOil)
    # Check julia_weathering config key exists in the class
    print("   ✓ julia_weathering config available")
except Exception as e:
    print(f"   ✗ Import failed: {e}")
    sys.exit(1)
PYEOF

PYTHON_RESULT=$?
if [[ $PYTHON_RESULT -ne 0 ]]; then
    echo "   ✗ Python verification FAILED"
    exit 1
fi

# ── Step 8: Generate deployment report ─────────────────────────────────────────
echo ""
echo "8. Deployment report..."

REPORT="$DEPLOY_DIR/IXACHI_DEPLOYMENT_REPORT.txt"
cat > "$REPORT" << EOF
CICOILv2 + CICOILPhysics.jl — Ixachi Deployment Report
========================================================

Date: $(date)
Host: $(hostname)
User: $(whoami)

Julia: $(julia --version)
Python: $(python --version 2>&1)
OpenDrift: $(python -c 'import opendrift; print(opendrift.__version__)' 2>/dev/null)
GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "N/A")

Julia Backend: $(julia --project="$JULIA_PROJECT" -e 'push!(LOAD_PATH, "'$JULIA_PROJECT'"); using CICOILPhysics; println(backend_name(backend()))' 2>/dev/null)

Paths:
  Deploy dir:    $DEPLOY_DIR
  Julia project: $JULIA_PROJECT
  OD site:       $OD_OPENOIL
  Conda env:     $CONDA_ENV

Usage:
  conda activate $CONDA_ENV
  export JULIA_PROJECT=$JULIA_PROJECT

  # In Python run script:
  from opendrift.models.openoil.ciceseoil import OpenCiceseOil
  o = OpenCiceseOil(weathering_model='cicese')
  o.set_config('processes:julia_weathering', True)
  o.set_config('julia:project_path', '$JULIA_PROJECT')
  o.set_config('processes:biodegradation', True)
  o.set_config('processes:photooxidation', True)
  # ... configure and run

EOF

echo "   ✓ Report saved to $REPORT"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "================================================"
echo "✓ DEPLOYMENT COMPLETE"
echo "================================================"
echo ""
echo "To use in a simulation:"
echo ""
echo "  conda activate $CONDA_ENV"
echo ""
echo "  # In your Python run script:"
echo "  o.set_config('processes:julia_weathering', True)"
echo "  o.set_config('julia:project_path', '$JULIA_PROJECT')"
echo ""
echo "Julia will auto-detect the H200 GPU and use CUDA kernels."
echo ""
