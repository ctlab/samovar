#!/bin/bash
# Conda-friendly installer for the Python SamovaR pipeline.
# R / samovaR is optional (generative abundance models on the r-package branch).
set -euo pipefail

echo "Installing SamovaR (Python pipeline)..."

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p build bin

if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/python" ]; then
    PYTHON_PATH="${PYTHON_PATH:-${CONDA_PREFIX}/bin/python}"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_PATH="${PYTHON_PATH:-python3}"
else
    PYTHON_PATH="${PYTHON_PATH:-python}"
fi

if [ -f build/config.json ]; then
    CFG_PY=$(grep -o '"python_path": *"[^"]*"' build/config.json | sed 's/"python_path": *"\(.*\)"/\1/' || true)
    if [ -n "${CFG_PY:-}" ] && command -v "$CFG_PY" >/dev/null 2>&1; then
        PYTHON_PATH="$CFG_PY"
    fi
fi

if ! command -v "$PYTHON_PATH" >/dev/null 2>&1; then
    echo "Python 3 is not installed. Create a conda env first, e.g.:"
    echo "  conda env create -f environment.yml && conda activate samovar"
    exit 1
fi

PYTHON_VERSION=$("$PYTHON_PATH" --version)
echo "Using $PYTHON_VERSION ($PYTHON_PATH)"
if [ -n "${CONDA_PREFIX:-}" ]; then
    echo "Conda prefix: $CONDA_PREFIX"
fi

echo "Installing Python package (editable)..."
"$PYTHON_PATH" -m pip install --upgrade pip
"$PYTHON_PATH" -m pip install -e ".[dev]"

cat > build/config.json << EOF
{
    "python_path": "$(command -v "$PYTHON_PATH")",
    "r_path": "$(command -v R 2>/dev/null || echo "")",
    "r_lib_path": ""
}
EOF

chmod +x bin/* workflow/database_prep/samovar_build_database.sh workflow/database_prep/samovar_build_database.py workflow/compare_annotations.py 2>/dev/null || true

echo "Building C++ annotation combiner..."
if command -v g++ >/dev/null 2>&1; then
    make -C src/cpp
    chmod +x bin/samovar_combine_annotations 2>/dev/null || true
else
    echo "Warning: g++ not found; annotation merge will try to compile on first use."
fi

if [ "${SAMOVAR_INSTALL_R:-0}" != "0" ] && command -v R >/dev/null 2>&1; then
    echo "SAMOVAR_INSTALL_R=1: installing optional R package..."
    R --quiet -e "if (!require('remotes')) install.packages('remotes', repos='https://cloud.r-project.org/')"
    R --quiet -e "library(remotes); remotes::install_local('.', upgrade='never')"
else
    echo "Skipping R package (not required for samovar prepare/exec). Set SAMOVAR_INSTALL_R=1 to install it."
fi

if [ -z "${CI:-}" ] && [ -f "${HOME}/.bashrc" ]; then
    if ! grep -q "export PATH=\$PATH:${ROOT}/bin" "${HOME}/.bashrc"; then
        echo "export PATH=\$PATH:${ROOT}/bin" >> "${HOME}/.bashrc"
        echo "Added ${ROOT}/bin to PATH in ~/.bashrc"
    fi
fi

echo "Installation completed successfully!"
echo "Config: $ROOT/build/config.json"
