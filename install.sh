#!/bin/bash
# Conda-friendly installer for the Python SamovaR pipeline.
# R / samovaR is optional (generative abundance models on the r-package branch).
#
# Environment:
#   SAMOVAR_OFFLINE=1          pip --offline (air-gapped; optional SAMOVAR_WHEELHOUSE)
#   SAMOVAR_INSTALL_DEV=1      also install pytest/flake8 extras
#   SAMOVAR_INSTALL_R=1        optional R package
#   NCBI_EMAIL / ENTREZ_EMAIL / SAMOVAR_EMAIL
#   CI=true                    non-interactive; NCBI_EMAIL defaults to test@samovar.com
set -euo pipefail

echo "Installing SamovaR (Python pipeline)..."

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p build bin
XDG_CONFIG_HOME="${XDG_CONFIG_HOME:-$HOME/.config}"
USER_CFG_DIR="${XDG_CONFIG_HOME}/samovar"
mkdir -p "$USER_CFG_DIR"

if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/python" ]; then
    PYTHON_PATH="${PYTHON_PATH:-${CONDA_PREFIX}/bin/python}"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_PATH="${PYTHON_PATH:-python3}"
else
    PYTHON_PATH="${PYTHON_PATH:-python}"
fi

if ! command -v "$PYTHON_PATH" >/dev/null 2>&1; then
    echo "Python 3 is not installed. Create a conda env first, e.g.:"
    echo "  conda env create -f environment.yml && conda activate samovar"
    exit 1
fi
PYTHON_PATH="$(command -v "$PYTHON_PATH" || true)"
PYTHON_PATH="${PYTHON_PATH:-python3}"

PYTHON_VERSION=$("$PYTHON_PATH" --version)
echo "Using $PYTHON_VERSION ($PYTHON_PATH)"
if [ -n "${CONDA_PREFIX:-}" ]; then
    echo "Conda prefix: $CONDA_PREFIX"
fi

# NCBI Entrez email (genome fetch)
if [ -n "${CI:-}" ]; then
    NCBI_EMAIL="${NCBI_EMAIL:-test@samovar.com}"
fi
NCBI_EMAIL="${NCBI_EMAIL:-${ENTREZ_EMAIL:-${SAMOVAR_EMAIL:-}}}"
if [ -z "$NCBI_EMAIL" ]; then
    if [ -t 0 ] && [ -z "${CI:-}" ]; then
        printf "NCBI Entrez email (required for genome downloads): "
        read -r NCBI_EMAIL
    fi
fi
if [ -z "$NCBI_EMAIL" ]; then
    NCBI_EMAIL="anonymous@example.com"
    echo "Warning: no NCBI email set; using $NCBI_EMAIL. Re-run with NCBI_EMAIL=you@institution.edu"
fi
export NCBI_EMAIL
echo "NCBI email: $NCBI_EMAIL"

PIP_OPTS=()
if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
    echo "Offline mode (SAMOVAR_OFFLINE=1): skipping pip upgrade / PyPI"
    PIP_OPTS+=(--offline)
    if [ -n "${SAMOVAR_WHEELHOUSE:-}" ]; then
        PIP_OPTS+=(--no-index --find-links "$SAMOVAR_WHEELHOUSE")
    fi
else
    echo "Installing Python package (editable)..."
    "$PYTHON_PATH" -m pip install --upgrade pip
fi

if [ "${SAMOVAR_INSTALL_DEV:-0}" != "0" ] || [ -n "${CI:-}" ]; then
    echo "Installing package + dev extras..."
    "$PYTHON_PATH" -m pip install -e ".[dev]" "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}" \
        || "$PYTHON_PATH" -m pip install -e ".[dev]"
else
    echo "Installing package (no dev extras; set SAMOVAR_INSTALL_DEV=1 for pytest/flake8)..."
    "$PYTHON_PATH" -m pip install -e . "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}" \
        || "$PYTHON_PATH" -m pip install -e .
fi

# ISS CLI (InSilicoSeq)
if ! command -v iss >/dev/null 2>&1; then
    echo "iss not on PATH; installing insilicoseq..."
    if [ "${SAMOVAR_OFFLINE:-0}" = "0" ]; then
        "$PYTHON_PATH" -m pip install "insilicoseq>2.0.0" || true
    fi
fi
if command -v iss >/dev/null 2>&1; then
    echo "iss: $(command -v iss)"
else
    echo "Warning: iss CLI still not found. Install insilicoseq or add it to PATH."
fi

chmod +x bin/* workflow/database_prep/samovar_build_database.sh workflow/database_prep/samovar_build_database.py workflow/compare_annotations.py 2>/dev/null || true

echo "Building C++ annotation combiner..."
CXX_BIN="${CXX:-}"
if [ -z "$CXX_BIN" ] || ! command -v "$CXX_BIN" >/dev/null 2>&1; then
    for cand in g++ c++ clang++; do
        if command -v "$cand" >/dev/null 2>&1; then
            CXX_BIN="$(command -v "$cand")"
            break
        fi
    done
fi
if [ -n "${CXX_BIN:-}" ]; then
    echo "C++ compiler: $CXX_BIN"
    make -C src/cpp CXX="$CXX_BIN"
    chmod +x bin/samovar_combine_annotations 2>/dev/null || true
else
    echo "Warning: no C++ compiler (g++/c++/clang++); annotation merge will try to compile on first use."
fi

if [ "${SAMOVAR_INSTALL_R:-0}" != "0" ] && command -v R >/dev/null 2>&1; then
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping R package install in offline mode."
    else
        echo "SAMOVAR_INSTALL_R=1: installing optional R package..."
        R --quiet -e "if (!require('remotes')) install.packages('remotes', repos='https://cloud.r-project.org/')"
        R --quiet -e "library(remotes); remotes::install_local('.', upgrade='never')"
    fi
else
    echo "Skipping R package (not required for samovar prepare/exec). Set SAMOVAR_INSTALL_R=1 to install it."
fi

# Write ~/.config/samovar/config.json (and build/config.json copy)
export SAMOVAR_ROOT="$ROOT"
export PYTHON_PATH
"$PYTHON_PATH" - <<'PY'
import os, shutil
from samovar.paths import discover_tools, write_config, PACKAGE_VERSION

root = os.environ["SAMOVAR_ROOT"]
python_path = os.environ["PYTHON_PATH"]
tools = discover_tools()
iss = shutil.which("iss") or tools.get("iss", "")
write_config({
    "version": PACKAGE_VERSION,
    "root": root,
    "python_path": python_path,
    "r_path": shutil.which("R") or "",
    "r_lib_path": "",
    "iss_path": iss,
    "ncbi_email": os.environ.get("NCBI_EMAIL", ""),
    "test_genomes": os.path.join(root, "data", "test_genomes"),
    "tools": tools,
})
print("Wrote user + repo config")
PY

echo "export NCBI_EMAIL=\"$NCBI_EMAIL\"" > "$USER_CFG_DIR/env"
echo "export SAMOVAR_ROOT=\"$ROOT\"" >> "$USER_CFG_DIR/env"

if [ -z "${CI:-}" ] && [ -f "${HOME}/.bashrc" ]; then
    if ! grep -q "export PATH=\$PATH:${ROOT}/bin" "${HOME}/.bashrc"; then
        echo "export PATH=\$PATH:${ROOT}/bin" >> "${HOME}/.bashrc"
        echo "Added ${ROOT}/bin to PATH in ~/.bashrc"
    fi
fi

echo "Running install smoke test..."
SMOKE_FAIL=0
"$PYTHON_PATH" -c "from samovar.paths import smoke_test; p=smoke_test();
print('\\n'.join('  WARN: '+x for x in p) or '  samovar.paths.smoke_test: ok');
raise SystemExit(1 if p else 0)" || SMOKE_FAIL=1
"$PYTHON_PATH" -c "import samovar; print('  import samovar: ok')" || SMOKE_FAIL=1
if command -v iss >/dev/null 2>&1; then
    echo "  iss: $(command -v iss)"
else
    echo "  iss: MISSING"; SMOKE_FAIL=1
fi
if command -v snakemake >/dev/null 2>&1; then
    echo "  snakemake: $(snakemake --version 2>/dev/null | head -1)"
else
    echo "  snakemake: MISSING (install via conda environment.yml)"; SMOKE_FAIL=1
fi
if [ ! -f "$ROOT/workflow/annotators/Snakefile" ]; then
    echo "  workflow: MISSING"; SMOKE_FAIL=1
fi
if [ "$SMOKE_FAIL" != "0" ]; then
    echo "Smoke test reported missing pieces (see above). Core Python package may still work."
else
    echo "Smoke test passed."
fi

echo "Installation completed successfully!"
echo "Config: $USER_CFG_DIR/config.json"
echo "Also:   $ROOT/build/config.json"
echo "Use NCBI_EMAIL=$NCBI_EMAIL (export it in new shells, or: source $USER_CFG_DIR/env)"
