#!/bin/bash
# Conda-friendly installer for the Python SamovaR pipeline.
# Optional R package lives on the GitHub ``r-package`` branch (not in this tree).
#
# Usage:
#   ./install.sh                 Python pipeline only
#   ./install.sh R-package       install samovaR from GitHub branch r-package
#   ./install.sh OPAL            optional CAMI OPAL (https://github.com/CAMI-challenge/OPAL)
#   ./install.sh MultiQC         optional MultiQC HTML reports
#   ./install.sh CAMISIM         optional CAMISIM community/read simulator
#   ./install.sh NanoSim         optional NanoSim in a separate conda env (ONT / hybrid)
#   ./install.sh ART             optional ART Illumina simulator in a separate conda env
#   ./install.sh seqtk           optional seqtk (FASTQ subsample / rarefaction)
#   ./install.sh CAMISIM NanoSim ART   several optionals without reinstalling the core
#
# Environment:
#   SAMOVAR_OFFLINE=1          pip --offline (air-gapped; optional SAMOVAR_WHEELHOUSE)
#   SAMOVAR_INSTALL_DEV=1      also install pytest/flake8 extras
#   SAMOVAR_INSTALL_R=1        also install optional R package (same as ./install.sh R-package)
#   SAMOVAR_INSTALL_OPAL=1     also install optional CAMI OPAL (same as ./install.sh OPAL)
#   SAMOVAR_INSTALL_MULTIQC=1  also install optional MultiQC (same as ./install.sh MultiQC)
#   SAMOVAR_INSTALL_CAMISIM=1  also clone optional CAMISIM (same as ./install.sh CAMISIM)
#   SAMOVAR_INSTALL_NANOSIM=1  also sidecar-conda NanoSim (same as ./install.sh NanoSim)
#   SAMOVAR_INSTALL_ART=1      also sidecar-conda ART (same as ./install.sh ART)
#   SAMOVAR_INSTALL_SEQTK=1    also seqtk (same as ./install.sh seqtk)
#   SAMOVAR_CONDA              conda/mamba/micromamba executable for sidecar envs
#   SAMOVAR_R_REPO             default ctlab/samovar
#   SAMOVAR_R_BRANCH           default r-package
#   SAMOVAR_UPDATE_SHELL=1     default: add bin/ to PATH this session and ~/.bashrc
#   SAMOVAR_UPDATE_SHELL=0     shared/read-only: skip PATH and ~/.bashrc edits
#   SAMOVAR_CONFIG=/path/to/config.json
#                              write main config to this file (dir holds env/)
#   SAMOVAR_CONFIG_DIR=/dir    write $dir/config.json instead of ~/.config/samovar
#   After install, the absolute config path is stored in build/config_path
#   NCBI_EMAIL / ENTREZ_EMAIL / SAMOVAR_EMAIL
#   CI=true                    non-interactive; NCBI_EMAIL defaults to test@samovar.com
#                              and shell/PATH edits are skipped
set -euo pipefail

echo "Installing SamovaR (Python pipeline)..."

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p bin
if ! mkdir -p build 2>/dev/null; then
    echo "Warning: cannot write $ROOT/build (read-only checkout); config pointer may be skipped"
fi

# Main config location (default: $XDG_CONFIG_HOME/samovar/config.json).
XDG_CONFIG_HOME="${XDG_CONFIG_HOME:-$HOME/.config}"
if [ -n "${SAMOVAR_CONFIG:-}" ]; then
    # Absolute file path for config.json
    case "${SAMOVAR_CONFIG}" in
        /*) ;;
        ~*) SAMOVAR_CONFIG="${SAMOVAR_CONFIG/#\~/$HOME}" ;;
        *) SAMOVAR_CONFIG="$(pwd)/${SAMOVAR_CONFIG}" ;;
    esac
    USER_CFG_DIR="$(dirname "$SAMOVAR_CONFIG")"
elif [ -n "${SAMOVAR_CONFIG_DIR:-}" ]; then
    case "${SAMOVAR_CONFIG_DIR}" in
        /*) USER_CFG_DIR="${SAMOVAR_CONFIG_DIR}" ;;
        ~*) USER_CFG_DIR="${SAMOVAR_CONFIG_DIR/#\~/$HOME}" ;;
        *) USER_CFG_DIR="$(pwd)/${SAMOVAR_CONFIG_DIR}" ;;
    esac
    SAMOVAR_CONFIG="${USER_CFG_DIR}/config.json"
else
    USER_CFG_DIR="${XDG_CONFIG_HOME}/samovar"
    SAMOVAR_CONFIG="${USER_CFG_DIR}/config.json"
fi
export SAMOVAR_CONFIG
export USER_CFG_DIR
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
PY_BIN="$(cd "$(dirname "$PYTHON_PATH")" && pwd)"
case ":$PATH:" in
    *":$PY_BIN:"*) ;;
    *) export PATH="$PY_BIN:$PATH" ;;
esac

install_samovar_r_package() {
    local repo="${SAMOVAR_R_REPO:-ctlab/samovar}"
    local branch="${SAMOVAR_R_BRANCH:-r-package}"
    local script_dest="$USER_CFG_DIR/annotation_regenerate.R"

    if ! command -v R >/dev/null 2>&1; then
        echo "R is not on PATH; cannot install the optional samovaR package."
        return 1
    fi

    local info
    info="$(R --vanilla -s -e 'if (requireNamespace("samovaR", quietly=TRUE)) { cat("INSTALLED", as.character(packageVersion("samovaR")), find.package("samovaR")) } else cat("MISSING")' 2>/dev/null || true)"
    if echo "$info" | grep -q '^INSTALLED'; then
        local ver lib
        ver="$(echo "$info" | awk '{print $2}')"
        lib="$(echo "$info" | awk '{print $3}')"
        echo "Warning: samovaR is already installed (version ${ver:-unknown}${lib:+ at $lib})."
        echo "Skipping GitHub install. Uninstall the package (not its dependencies) to reinstall from $repo@$branch."
    else
        if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
            echo "Skipping samovaR install in offline mode."
            return 1
        fi
        echo "Installing samovaR from https://github.com/${repo}/tree/${branch} ..."
        R --vanilla -s -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org')"
        if ! R --vanilla -s -e "library(remotes); remotes::install_github('${repo}', ref='${branch}', upgrade='never', dependencies=NA)"; then
            echo "GitHub install failed; trying a local checkout of origin/${branch}..."
            local tmp
            tmp="$(mktemp -d)"
            if git -C "$ROOT" fetch origin "$branch" 2>/dev/null; then
                git -C "$ROOT" archive --format=tar "origin/${branch}" | tar -x -C "$tmp"
                R --vanilla -s -e "library(remotes); remotes::install_local('${tmp}', upgrade='never', dependencies=NA)"
            else
                echo "Could not install samovaR from ${repo}@${branch}"
                rm -rf "$tmp"
                return 1
            fi
            rm -rf "$tmp"
        fi
        info="$(R --vanilla -s -e 'if (requireNamespace("samovaR", quietly=TRUE)) { cat("INSTALLED", as.character(packageVersion("samovaR"))) } else cat("MISSING")' 2>/dev/null || true)"
        echo "samovaR: $info"
    fi

    mkdir -p "$USER_CFG_DIR"
    export SAMOVAR_R_SCRIPT="$script_dest"
    "$PYTHON_PATH" - <<'PY'
import os, shutil, sys
from pathlib import Path

dest = Path(os.environ["SAMOVAR_R_SCRIPT"])
try:
    from samovar.table2iss import R_REGENERATE_DRIVER
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(R_REGENERATE_DRIVER, encoding="utf-8")
    print("R regenerator script:", dest)
except Exception as exc:
    print("Warning: could not write R driver:", exc)
    sys.exit(0)
try:
    from samovar.paths import load_config, write_config
except ImportError:
    print(
        "Warning: samovar Python package not importable; set SAMOVAR_R_REGENERATE="
        + str(dest)
    )
    sys.exit(0)
cfg = load_config()
cfg["annotation_regenerate_r"] = str(dest)
cfg["r_path"] = shutil.which("R") or cfg.get("r_path") or ""
write_config(cfg)
print("Updated config annotation_regenerate_r")
PY
    return 0
}

install_opal() {
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping OPAL install in offline mode."
        return 1
    fi
    echo "Installing optional CAMI OPAL (https://github.com/CAMI-challenge/OPAL) ..."
    # cami-opal pins exact numpy/pandas/scipy; those pins fight the SamovaR env
    # and can hang pip. Install the package only; scientific stack is already here.
    if ! "$PYTHON_PATH" -m pip install "cami-opal" --no-deps "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}"; then
        echo "Warning: pip install cami-opal failed. Profiling HTML from OPAL will be skipped."
        echo "Install later with: ./install.sh OPAL"
        return 1
    fi
    "$PYTHON_PATH" -m pip install "bokeh>=3" "jinja2>=3" "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}" || true
    hash -r 2>/dev/null || true
    export PY_BIN
    "$PYTHON_PATH" - <<'PY'
import os, shutil, stat, sys
from pathlib import Path

try:
    from samovar.paths import discover_opal, load_config, write_config
except ImportError:
    print("Warning: samovar is not importable yet; re-run ./install.sh to record opal_path")
    raise SystemExit(0)

found = discover_opal()
py_bin = Path(os.environ.get("PY_BIN") or Path(sys.executable).parent)
for extra in (py_bin / "opal.py", Path.home() / ".local" / "bin" / "opal.py"):
    if extra.is_file() and (not found):
        found = str(extra.resolve())
        break
if not found:
    print("Warning: cami-opal installed but opal.py was not found on PATH or next to Python.")
    print("Add it with: opal_path in ~/.config/samovar/config.json")
    raise SystemExit(1)
path = Path(found)
if not os.access(path, os.X_OK):
    try:
        path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    except OSError:
        pass
cfg = load_config()
cfg["opal_path"] = str(path)
tools = dict(cfg.get("tools") or {})
tools["opal.py"] = str(path)
cfg["tools"] = tools
write_config(cfg)
print("OPAL:", path)
print("Updated config opal_path:", cfg.get("opal_path"))
PY
}

install_multiqc() {
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping MultiQC install in offline mode."
        return 1
    fi
    echo "Installing optional MultiQC (https://seqera.io/multiqc/) ..."
    if ! "$PYTHON_PATH" -m pip install "multiqc>=1.21" "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}"; then
        echo "Warning: pip install multiqc failed. HTML reports from MultiQC will be skipped."
        echo "Install later with: ./install.sh MultiQC"
        return 1
    fi
    hash -r 2>/dev/null || true
    export PY_BIN
    "$PYTHON_PATH" - <<'PY'
import os, shutil, sys
from pathlib import Path

try:
    from samovar.paths import discover_multiqc, load_config, write_config
except ImportError:
    print("Warning: samovar is not importable yet; re-run ./install.sh to record multiqc_path")
    raise SystemExit(0)

found = discover_multiqc() or shutil.which("multiqc")
py_bin = Path(os.environ.get("PY_BIN") or Path(sys.executable).parent)
for extra in (py_bin / "multiqc", Path.home() / ".local" / "bin" / "multiqc"):
    if extra.is_file() and (not found):
        found = str(extra.resolve())
        break
if not found:
    print("Warning: multiqc installed but the CLI was not found on PATH.")
    print("Add it with: multiqc_path in ~/.config/samovar/config.json")
    raise SystemExit(1)
cfg = load_config()
cfg["multiqc_path"] = str(found)
tools = dict(cfg.get("tools") or {})
tools["multiqc"] = str(found)
cfg["tools"] = tools
write_config(cfg)
print("MultiQC:", found)
print("Updated config multiqc_path:", cfg.get("multiqc_path"))
PY
}

install_camisim() {
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping CAMISIM clone in offline mode."
        return 1
    fi
    echo "Installing optional CAMISIM (https://github.com/CAMI-challenge/CAMISIM) ..."
    if ! command -v git >/dev/null 2>&1; then
        echo "Warning: git is required to clone CAMISIM."
        return 1
    fi
    export PYTHON_PATH
    "$PYTHON_PATH" - <<'PY'
import os, shutil, subprocess, sys
from pathlib import Path

try:
    from samovar.paths import load_config, write_config
except ImportError:
    print("Warning: samovar is not importable yet; re-run ./install.sh CAMISIM after the Python package is installed")
    raise SystemExit(0)

dest = Path(
    os.environ.get("SAMOVAR_CAMISIM")
    or (Path(os.environ["SAMOVAR_ROOT"]) / ".cache" / "CAMISIM")
)
if (dest / "main.nf").is_file() or (dest / "metagenomesimulation.py").is_file():
    print("CAMISIM already present:", dest)
else:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not any(dest.iterdir()):
        dest.rmdir()
    print("Cloning CAMISIM into", dest)
    subprocess.check_call(
        ["git", "clone", "--depth", "1", "https://github.com/CAMI-challenge/CAMISIM", str(dest)]
    )
cfg = load_config()
cfg["camisim_path"] = str(dest)
nxt = shutil.which("nextflow") or ""
if nxt:
    cfg["nextflow_path"] = nxt
tools = dict(cfg.get("tools") or {})
tools["camisim"] = str(dest)
if nxt:
    tools["nextflow"] = nxt
cfg["tools"] = tools
write_config(cfg)
print("CAMISIM:", dest)
print("Nextflow:", nxt or "not found (install nextflow to run CAMISIM 2.0 read simulation)")
print("Updated config camisim_path")
if not nxt:
    print("Note: `samovar generate --camisim-mode table` still works (community design + ISS).")
    print("Illumina/ONT/wgsim/hybrid need Nextflow + CAMISIM simulators (ART / NanoSim / wgsim).")
    print("ONT/hybrid: ./install.sh NanoSim   (separate conda env; do not mix with SamovaR Python)")
    print("Illumina:   ./install.sh ART       (optional sidecar; Nextflow can also create ART via conda)")
PY
}

install_nanosim() {
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping NanoSim sidecar env in offline mode."
        return 1
    fi
    echo "Installing optional NanoSim in a separate conda env (Python 3.8; not the SamovaR env) ..."
    if ! "$PYTHON_PATH" -m samovar.sidecar nanosim; then
        echo "Warning: NanoSim sidecar install failed."
        echo "Install later with: ./install.sh NanoSim"
        echo "Or: conda create -p \"\$SAMOVAR_ROOT/.cache/samovar/envs/nanosim\" -c conda-forge -c bioconda python=3.8 nanosim=3.2 htseq gffutils numpy=1.23.5 scikit-learn=0.23.2"
        echo "Then set tool_envs.nanosim and nanosim_path in ~/.config/samovar/config.json"
        return 1
    fi
}

install_art() {
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping ART sidecar env in offline mode."
        return 1
    fi
    echo "Installing optional ART (art_illumina) in a separate conda env ..."
    if ! "$PYTHON_PATH" -m samovar.sidecar art; then
        echo "Warning: ART sidecar install failed."
        echo "Install later with: ./install.sh ART"
        echo "Or: conda create -p \"\$SAMOVAR_ROOT/.cache/samovar/envs/art\" -c bioconda -c conda-forge art samtools"
        echo "Then set tool_envs.art and art_path in ~/.config/samovar/config.json"
        return 1
    fi
}

install_seqtk() {
    if command -v seqtk >/dev/null 2>&1; then
        echo "seqtk: $(command -v seqtk)"
        return 0
    fi
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        echo "Skipping seqtk install in offline mode. conda install -c bioconda seqtk"
        return 1
    fi
    local conda_bin="${SAMOVAR_CONDA:-}"
    if [ -z "$conda_bin" ]; then
        conda_bin="$(command -v mamba || command -v conda || true)"
    fi
    if [ -z "$conda_bin" ]; then
        echo "seqtk is not on PATH and conda/mamba was not found."
        echo "Install with: conda install -c bioconda seqtk"
        return 1
    fi
    echo "Installing optional seqtk (https://github.com/lh3/seqtk) via bioconda ..."
    if ! "$conda_bin" install -y -c bioconda -c conda-forge seqtk; then
        echo "Warning: seqtk install failed. Rarefaction with seqtk sample will be unavailable."
        return 1
    fi
    hash -r 2>/dev/null || true
    if command -v seqtk >/dev/null 2>&1; then
        echo "seqtk: $(command -v seqtk)"
        return 0
    fi
    echo "seqtk installed but not on PATH in this shell."
    return 1
}

print_install_status() {
    "$PYTHON_PATH" - <<'PY' || true
try:
    from samovar.paths import format_install_status
    print("")
    print(format_install_status())
except Exception as exc:
    print("Could not print tool status:", exc)
PY
}

normalize_optional_arg() {
    local raw="${1:-}"
    local lower
    lower="$(printf '%s' "$raw" | tr '[:upper:]' '[:lower:]')"
    case "$lower" in
        r-package|rpackage|r) echo "r-package" ;;
        opal) echo "opal" ;;
        multiqc) echo "multiqc" ;;
        camisim) echo "camisim" ;;
        nanosim|nanosim3) echo "nanosim" ;;
        art|art_illumina) echo "art" ;;
        nextflow) echo "nextflow" ;;
        seqtk) echo "seqtk" ;;
        *) echo "" ;;
    esac
}

run_optional_named() {
    local name="$1"
    case "$name" in
        r-package) install_samovar_r_package ;;
        opal)
            PIP_OPTS=()
            if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
                PIP_OPTS+=(--offline --no-build-isolation)
                if [ -n "${SAMOVAR_WHEELHOUSE:-}" ]; then
                    PIP_OPTS+=(--no-index --find-links "$SAMOVAR_WHEELHOUSE")
                fi
            fi
            install_opal
            ;;
        multiqc)
            PIP_OPTS=()
            if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
                PIP_OPTS+=(--offline --no-build-isolation)
                if [ -n "${SAMOVAR_WHEELHOUSE:-}" ]; then
                    PIP_OPTS+=(--no-index --find-links "$SAMOVAR_WHEELHOUSE")
                fi
            fi
            install_multiqc
            ;;
        camisim) install_camisim ;;
        nanosim) install_nanosim ;;
        art) install_art ;;
        seqtk) install_seqtk ;;
        nextflow)
            echo "Nextflow is not bundled. Install with: conda install -c bioconda nextflow"
            echo "or https://www.nextflow.io/docs/latest/install.html"
            echo "Then set nextflow_path in ~/.config/samovar/config.json"
            if command -v nextflow >/dev/null 2>&1; then
                echo "Found: $(command -v nextflow)"
            else
                return 1
            fi
            ;;
        *) return 1 ;;
    esac
}

OPTIONAL_ONLY_ARGS=()
if [ "$#" -gt 0 ]; then
    ALL_OPTIONAL=1
    for arg in "$@"; do
        mapped="$(normalize_optional_arg "$arg")"
        if [ -z "$mapped" ]; then
            ALL_OPTIONAL=0
            break
        fi
        OPTIONAL_ONLY_ARGS+=("$mapped")
    done
    if [ "$ALL_OPTIONAL" = "1" ] && [ "${#OPTIONAL_ONLY_ARGS[@]}" -gt 0 ]; then
        FAIL=0
        for name in "${OPTIONAL_ONLY_ARGS[@]}"; do
            echo "Optional install: $name"
            run_optional_named "$name" || FAIL=1
        done
        print_install_status
        exit "$FAIL"
    fi
fi

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
    PIP_OPTS+=(--offline --no-build-isolation)
    if [ -n "${SAMOVAR_WHEELHOUSE:-}" ]; then
        PIP_OPTS+=(--no-index --find-links "$SAMOVAR_WHEELHOUSE")
    else
        echo "Warning: SAMOVAR_WHEELHOUSE is unset; offline pip needs local wheels."
    fi
else
    echo "Installing Python package (editable)..."
    "$PYTHON_PATH" -m pip install --upgrade pip
fi

if [ "${SAMOVAR_INSTALL_DEV:-0}" != "0" ] || [ -n "${CI:-}" ]; then
    echo "Installing package + dev extras..."
    "$PYTHON_PATH" -m pip install -e ".[dev]" "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}"
else
    echo "Installing package (no dev extras; set SAMOVAR_INSTALL_DEV=1 for pytest/flake8)..."
    "$PYTHON_PATH" -m pip install -e . "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}"
fi

# ISS / snakemake: try pip, then hard-exit if still missing
ensure_cli() {
    local name="$1"
    local pip_spec="$2"
    local hint="$3"
    if command -v "$name" >/dev/null 2>&1; then
        echo "$name: $(command -v "$name")"
        return 0
    fi
    if [ -x "$PY_BIN/$name" ]; then
        echo "$name: $PY_BIN/$name"
        return 0
    fi
    echo "$name not on PATH; installing ${pip_spec}..."
    "$PYTHON_PATH" -m pip install "$pip_spec" "${PIP_OPTS[@]+"${PIP_OPTS[@]}"}" || true
    hash -r 2>/dev/null || true
    if command -v "$name" >/dev/null 2>&1; then
        echo "$name: $(command -v "$name")"
        return 0
    fi
    if [ -x "$PY_BIN/$name" ]; then
        echo "$name: $PY_BIN/$name"
        return 0
    fi
    echo "ERROR: $name is required but was not found after install attempt."
    echo "$hint"
    return 1
}

ensure_cli iss "insilicoseq>2.0.0" \
    "Install with: $PYTHON_PATH -m pip install 'insilicoseq>2.0.0' (conda: insilicoseq)" || exit 1
ensure_cli snakemake "snakemake>=7.0" \
    "Install with: $PYTHON_PATH -m pip install 'snakemake>=7.0' or conda install -c bioconda snakemake-minimal" || exit 1

chmod +x bin/* workflow/database_prep/samovar_build_database.sh workflow/database_prep/samovar_build_database.py workflow/compare_annotations.py workflow/annotation_regenerate.py workflow/combine_annotation_tables.py workflow/remap_taxids.py workflow/ML.py 2>/dev/null || true

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

if [ "${SAMOVAR_INSTALL_R:-0}" != "0" ]; then
    install_samovar_r_package || echo "Warning: optional R package install did not complete."
else
    echo "Skipping R package (not required). Use ./install.sh R-package to install samovaR from GitHub branch r-package."
fi

if [ "${SAMOVAR_INSTALL_OPAL:-0}" != "0" ]; then
    install_opal || echo "Warning: optional OPAL install did not complete."
else
    echo "Skipping OPAL (not required). Use ./install.sh OPAL to install CAMI OPAL (cami-opal)."
fi

if [ "${SAMOVAR_INSTALL_MULTIQC:-0}" != "0" ]; then
    install_multiqc || echo "Warning: optional MultiQC install did not complete."
else
    echo "Skipping MultiQC (not required). Use ./install.sh MultiQC to install MultiQC."
fi

if [ "${SAMOVAR_INSTALL_CAMISIM:-0}" != "0" ]; then
    install_camisim || echo "Warning: optional CAMISIM install did not complete."
else
    echo "Skipping CAMISIM (not required). Use ./install.sh CAMISIM to clone https://github.com/CAMI-challenge/CAMISIM"
fi

if [ "${SAMOVAR_INSTALL_NANOSIM:-0}" != "0" ]; then
    install_nanosim || echo "Warning: optional NanoSim sidecar install did not complete."
else
    echo "Skipping NanoSim (not required). Use ./install.sh NanoSim for a separate conda env (CAMISIM ONT/hybrid)."
fi

if [ "${SAMOVAR_INSTALL_ART:-0}" != "0" ]; then
    install_art || echo "Warning: optional ART sidecar install did not complete."
else
    echo "Skipping ART (not required). Use ./install.sh ART for a separate conda env (CAMISIM Illumina)."
fi

if [ "${SAMOVAR_INSTALL_SEQTK:-0}" != "0" ]; then
    install_seqtk || echo "Warning: optional seqtk install did not complete."
else
    echo "Skipping seqtk (not required). Use ./install.sh seqtk for FASTQ subsample / rarefaction."
fi

# Write main config.json (location: $SAMOVAR_CONFIG) and build/config_path
export USER_CFG_DIR
export SAMOVAR_CONFIG
export SAMOVAR_ROOT="$ROOT"
export SAMOVAR_ROOT="$ROOT"
export PYTHON_PATH
"$PYTHON_PATH" - <<'PY'
import os, shutil
from pathlib import Path
from samovar.main_config import build_install_config
from samovar.paths import (
    PACKAGE_VERSION,
    collect_runtime_path_dirs,
    discover_multiqc,
    discover_opal,
    discover_tools,
    load_config,
    write_config,
)
try:
    from samovar.camisim import discover_camisim
except Exception:
    def discover_camisim():
        return None
try:
    from samovar.paths import discover_art, discover_nanosim, cxx_compiler
except Exception:
    def discover_art():
        return None
    def discover_nanosim():
        return None
    def cxx_compiler():
        return None

root = os.environ["SAMOVAR_ROOT"]
python_path = os.environ["PYTHON_PATH"]
existing = load_config()
tools = dict(discover_tools())
# Preserve previously configured tool paths (any schema).
from samovar.main_config import iter_tools, tool_path as _tool_path
for name, spec in iter_tools(existing).items():
    path = _tool_path(spec, name)
    if path:
        tools.setdefault(name, path)
for name, path in (existing.get("tools") or {}).items():
    if isinstance(path, str) and path.strip():
        tools.setdefault(name, path)

_home = Path.home().resolve()

def _not_under_home(path: Path) -> bool:
    try:
        return not path.expanduser().resolve().is_relative_to(_home)
    except (OSError, ValueError):
        return True

default_genomes = ""
if os.environ.get("SAMOVAR_GENOMES", "").strip():
    default_genomes = os.environ["SAMOVAR_GENOMES"].strip()
elif os.environ.get("XDG_CACHE_HOME", "").strip():
    default_genomes = str(Path(os.environ["XDG_CACHE_HOME"]) / "samovar" / "genomes")

def _usable_dir(raw, *, allow_home: bool = False):
    text = str(raw or "").strip()
    if not text:
        return ""
    path = Path(text).expanduser()
    if not allow_home and not _not_under_home(path):
        return ""
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError:
        return ""
    return str(path)

from samovar.main_config import first_dir, processed_genome_dirs, raw_genome_dirs
genomes = _usable_dir(first_dir(raw_genome_dirs(existing))) or _usable_dir(default_genomes) or ""
processed = (
    _usable_dir(first_dir(processed_genome_dirs(existing)))
    or _usable_dir(first_dir(raw_genome_dirs(existing)))
    or genomes
)
kept_libs = []
from samovar.main_config import extra_genome_dirs as _extra_dirs
for item in _extra_dirs(existing) or existing.get("genome_dirs") or []:
    text = str(item or "").strip()
    if not text:
        continue
    path = Path(text).expanduser()
    try:
        ok = path.is_dir()
    except OSError:
        ok = False
    if ok and _not_under_home(path):
        kept_libs.append(str(path))

path_extra = existing.get("path") or existing.get("extra_path") or []
from samovar.main_config import compiler_python_libs
path_extra = path_extra or compiler_python_libs(existing)

sidecars = []
from samovar.main_config import tool_env_prefix as _tep
for name in ("nanosim", "art"):
    if _tep((iter_tools(existing) or {}).get(name), name):
        sidecars.append(name)
if discover_nanosim():
    sidecars.append("nanosim")
if discover_art():
    sidecars.append("art")
sidecars = list(dict.fromkeys(sidecars))

discovered = dict(tools)
if discover_opal():
    discovered.setdefault("opal.py", discover_opal())
if discover_multiqc():
    discovered.setdefault("multiqc", discover_multiqc())
cam = discover_camisim()
if cam:
    discovered.setdefault("camisim", cam)
nano = discover_nanosim()
if nano:
    discovered.setdefault("nanosim", nano)
    discovered.setdefault("simulator.py", nano)
art = discover_art()
if art:
    discovered.setdefault("art_illumina", art)

payload = build_install_config(
    root=root,
    python_path=python_path,
    version=PACKAGE_VERSION,
    existing=existing,
    discovered_tools=discovered,
    ncbi_email=os.environ.get("NCBI_EMAIL", ""),
    genomes_default=genomes,
    processed_default=processed,
    extra_genome_dirs=kept_libs,
    extra_path=path_extra if isinstance(path_extra, list) else [path_extra] if path_extra else [],
    bash=shutil.which("bash") or "",
    cxx=(cxx_compiler() or shutil.which("g++") or ""),
    r_path=shutil.which("R") or "",
    conda_prefix=os.environ.get("CONDA_PREFIX") or "",
    conda_sidecars=sidecars,
)
write_config(payload)
from samovar.paths import user_config_path, write_install_config_pointer
cfg_path = user_config_path()
pointer = write_install_config_pointer(cfg_path, root=Path(root))
user_cfg = Path(os.environ["USER_CFG_DIR"])
user_cfg.mkdir(parents=True, exist_ok=True)
email = os.environ.get("NCBI_EMAIL", "")
dirs = collect_runtime_path_dirs()
path_export = ":".join(dirs + ["$PATH"])
conda_lib = ""
prefix = os.environ.get("CONDA_PREFIX") or ""
if prefix and Path(prefix, "lib").is_dir():
    conda_lib = f'export LD_LIBRARY_PATH="{prefix}/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"\n'
env_path = user_cfg / "env"
# Do not clobber a runtime NCBI_EMAIL; the saved install value is only the default.
env_path.write_text(
    f'export NCBI_EMAIL="${{NCBI_EMAIL:-{email}}}"\n'
    f'export SAMOVAR_ROOT="{root}"\n'
    f'export SAMOVAR_CONFIG="{cfg_path}"\n'
    f'export PYTHON_PATH="{python_path}"\n'
    f'export PATH="{path_export}"\n'
    f"{conda_lib}",
    encoding="utf-8",
)
print(f"Wrote config: {cfg_path}")
if pointer:
    print(f"Wrote pointer: {pointer}")
print(f"Wrote {env_path} PATH with {len(dirs)} tool bin dir(s)")
from samovar.main_config import extra_genome_dirs as _libs, first_dir as _first, processed_genome_dirs as _pg, raw_genome_dirs as _rg, iter_tools as _it, tool_path as _tp
genomes = _first(_rg(payload))
processed = _first(_pg(payload)) or genomes
lib_dirs = _libs(payload)
print("")
print("Genome cache (NCBI / user assemblies, reused by `samovar prepare`):")
print(f"  Downloaded genomes:  {genomes or '(per-run: $out_dir/.cache/samovar/genomes)'}")
print(f"  Processed genomes:   {processed or '(same as downloaded / per-run outdir)'}")
print(f"  Extra libraries:     {', '.join(lib_dirs) if lib_dirs else '(none yet; add folders under genomes.raw in config.json)'}")
print("")
print("IMPORTANT: bulky genome caches must NOT live under $HOME (quota).")
print("  Each `samovar exec` pins XDG_CACHE_HOME/SAMOVAR_GENOMES to $out_dir/.cache.")
print("  Optional shared scratch library: SAMOVAR_GENOMES=/scratch/... ./install.sh")
print("IMPORTANT: data/test_genomes contains TRUNCATED assemblies for ISS and CI only.")
print("They are never used as NCBI substitutes. Real metagenomes reuse only genomes")
print("previously downloaded by SamovaR or listed in genomes.raw.")
print("  samovar prepare --reuse-genomes      # default: symlink from cache if found")
print("  samovar prepare --no-reuse-genomes   # download into the genomes cache")
print("  samovar prepare --genome-dirs DIR[:DIR]")
print("  samovar prepare --test-genomes       # allow truncated stubs (ISS/CI only)")
print("Prefer symlinks into $out_dir/genomes; copy with a warning if linking fails.")
_tools = _it(payload)
if _tp(_tools.get("opal.py") or _tools.get("opal"), "opal.py"):
    print(f"OPAL (cami-opal): {_tp(_tools.get('opal.py') or _tools.get('opal'), 'opal.py')}")
else:
    print("OPAL: not found. ./install.sh OPAL to install CAMI OPAL and record tools.opal.py.")
if _tp(_tools.get("multiqc"), "multiqc"):
    print(f"MultiQC: {_tp(_tools.get('multiqc'), 'multiqc')}")
else:
    print("MultiQC: not found. ./install.sh MultiQC to install MultiQC and record tools.multiqc.")
if _tp(_tools.get("camisim"), "camisim"):
    print(f"CAMISIM: {_tp(_tools.get('camisim'), 'camisim')}")
else:
    print("CAMISIM: not found. ./install.sh CAMISIM")
if _tp(_tools.get("nanosim") or _tools.get("simulator.py"), "nanosim"):
    print(f"NanoSim: {_tp(_tools.get('nanosim') or _tools.get('simulator.py'), 'nanosim')}")
else:
    print("NanoSim: not found. ./install.sh NanoSim (separate conda env for CAMISIM ONT)")
if _tp(_tools.get("art_illumina") or _tools.get("art"), "art_illumina"):
    print(f"ART: {_tp(_tools.get('art_illumina') or _tools.get('art'), 'art_illumina')}")
else:
    print("ART: not found. ./install.sh ART (separate conda env) or let CAMISIM Nextflow create it")
PY

BIN_DIR="$ROOT/bin"
UPDATE_SHELL="${SAMOVAR_UPDATE_SHELL:-1}"
if [ -n "${CI:-}" ] || [ "${SAMOVAR_UPDATE_SHELL:-}" = "0" ]; then
    UPDATE_SHELL=0
fi

EXISTING_SAMOVAR="$(command -v samovar 2>/dev/null || true)"
if [ -n "$EXISTING_SAMOVAR" ]; then
    EXISTING_DIR="$(cd "$(dirname "$EXISTING_SAMOVAR")" && pwd)"
    if [ "$EXISTING_DIR" != "$BIN_DIR" ]; then
        echo "Warning: another samovar is already on PATH: $EXISTING_SAMOVAR"
        echo "This install is $BIN_DIR/samovar."
        if [ "$UPDATE_SHELL" != "0" ]; then
            echo "Switching to SAMOVAR_UPDATE_SHELL=0 (no ~/.bashrc or PATH edits)."
            UPDATE_SHELL=0
        fi
    fi
fi

if [ "$UPDATE_SHELL" != "0" ]; then
    case ":$PATH:" in
        *":$BIN_DIR:"*) ;;
        *)
            export PATH="$BIN_DIR:$PATH"
            echo "Added $BIN_DIR to PATH for this session"
            ;;
    esac
    if [ -f "${HOME}/.bashrc" ]; then
        if ! grep -Fq "$BIN_DIR" "${HOME}/.bashrc"; then
            echo "export PATH=\$PATH:${BIN_DIR}" >> "${HOME}/.bashrc"
            echo "Added ${BIN_DIR} to PATH in ~/.bashrc"
        fi
    fi
else
    echo "Shell/PATH edits skipped (SAMOVAR_UPDATE_SHELL=0 or another samovar on PATH)."
    echo "Load this install with: source $USER_CFG_DIR/env"
fi

if [ -n "${CONDA_PREFIX:-}" ] && [ -d "${CONDA_PREFIX}/lib" ]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

echo "Running install smoke test..."
SMOKE_FAIL=0
"$PYTHON_PATH" -c "from samovar.paths import smoke_test; p=smoke_test();
print('\\n'.join('  WARN: '+x for x in p) or '  samovar.paths.smoke_test: ok');
raise SystemExit(1 if p else 0)" || SMOKE_FAIL=1
"$PYTHON_PATH" -c "import samovar; print('  import samovar: ok')" || SMOKE_FAIL=1
"$PYTHON_PATH" -c "from samovar.viz_annotation import require_viz_backend; print('  viz backend:', require_viz_backend())" || {
    echo "ERROR: cnsplots (preferred) or altair is required for pipeline plots."
    echo "Re-run ./install.sh or: $PYTHON_PATH -m pip install 'cnsplots>=0.6.0' altair"
    exit 1
}
if command -v iss >/dev/null 2>&1 || [ -x "$PY_BIN/iss" ]; then
    echo "  iss: $(command -v iss 2>/dev/null || echo "$PY_BIN/iss")"
else
    echo "ERROR: iss CLI is required (tried pip install insilicoseq)."
    exit 1
fi
if command -v snakemake >/dev/null 2>&1 || [ -x "$PY_BIN/snakemake" ]; then
    echo "  snakemake: $(snakemake --version 2>/dev/null | head -1 || echo found)"
else
    echo "ERROR: snakemake is required (tried pip install snakemake)."
    echo "Or: conda install -c bioconda snakemake-minimal"
    exit 1
fi
if [ ! -f "$ROOT/workflow/annotators/Snakefile" ]; then
    echo "  workflow: MISSING"; SMOKE_FAIL=1
fi
"$PYTHON_PATH" -c "from samovar.opal import opal_executable; p=opal_executable(); print('  opal.py:', p or 'not found (optional; ./install.sh OPAL)')" || true
"$PYTHON_PATH" -c "from samovar.paths import discover_multiqc; p=discover_multiqc(); print('  multiqc:', p or 'not found (optional; ./install.sh MultiQC)')" || true
"$PYTHON_PATH" -c "from samovar.camisim import discover_camisim; p=discover_camisim(); print('  CAMISIM:', p or 'not found (optional; ./install.sh CAMISIM)')" || true
"$PYTHON_PATH" -c "from samovar.paths import discover_nanosim; p=discover_nanosim(); print('  NanoSim:', p or 'not found (optional; ./install.sh NanoSim)')" || true
"$PYTHON_PATH" -c "from samovar.paths import discover_art; p=discover_art(); print('  ART:', p or 'not found (optional; ./install.sh ART)')" || true
if [ "$SMOKE_FAIL" != "0" ]; then
    echo "Smoke test reported missing pieces (see above). Core Python package may still work."
else
    echo "Smoke test passed."
fi

print_install_status

echo "Installation completed successfully!"
echo "Config: $SAMOVAR_CONFIG"
if [ -f "$ROOT/build/config_path" ]; then
    echo "Pointer: $ROOT/build/config_path -> $(tr -d '[:space:]' < "$ROOT/build/config_path")"
fi
if [ -f "$ROOT/build/config.json" ] && [ "$(realpath -m "$ROOT/build/config.json" 2>/dev/null || echo "$ROOT/build/config.json")" != "$(realpath -m "$SAMOVAR_CONFIG" 2>/dev/null || echo "$SAMOVAR_CONFIG")" ]; then
    echo "Mirror:  $ROOT/build/config.json"
fi
echo "Use NCBI_EMAIL=$NCBI_EMAIL (export it in new shells, or: source $USER_CFG_DIR/env)"
echo "Custom config location: SAMOVAR_CONFIG=/path/config.json ./install.sh"
echo "Shared/HPC install without PATH edits: SAMOVAR_UPDATE_SHELL=0 ./install.sh"
