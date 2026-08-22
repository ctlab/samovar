#!/bin/bash
# Conda-friendly installer for the Python SamovaR pipeline.
# Optional R package lives on the GitHub ``r-package`` branch (not in this tree).
#
# Usage:
#   ./install.sh                 Python pipeline only
#   ./install.sh R-package       install samovaR from GitHub branch r-package
#   ./install.sh OPAL            optional CAMI OPAL (https://github.com/CAMI-challenge/OPAL)
#
# Environment:
#   SAMOVAR_OFFLINE=1          pip --offline (air-gapped; optional SAMOVAR_WHEELHOUSE)
#   SAMOVAR_INSTALL_DEV=1      also install pytest/flake8 extras
#   SAMOVAR_INSTALL_R=1        also install optional R package (same as ./install.sh R-package)
#   SAMOVAR_INSTALL_OPAL=1     also install optional CAMI OPAL (same as ./install.sh OPAL)
#   SAMOVAR_R_REPO             default ctlab/samovar
#   SAMOVAR_R_BRANCH           default r-package
#   SAMOVAR_UPDATE_SHELL=1     default: add bin/ to PATH this session and ~/.bashrc
#   SAMOVAR_UPDATE_SHELL=0     shared/read-only: skip PATH and ~/.bashrc edits
#   NCBI_EMAIL / ENTREZ_EMAIL / SAMOVAR_EMAIL
#   CI=true                    non-interactive; NCBI_EMAIL defaults to test@samovar.com
#                              and shell/PATH edits are skipped
set -euo pipefail

echo "Installing SamovaR (Python pipeline)..."

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p bin
if ! mkdir -p build 2>/dev/null; then
    echo "Warning: cannot write $ROOT/build (read-only checkout); config will live under XDG only"
fi
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

if [ "${1:-}" = "R-package" ] || [ "${1:-}" = "r-package" ]; then
    echo "R-package mode (Python tree is unchanged)."
    install_samovar_r_package
    exit $?
fi

if [ "${1:-}" = "OPAL" ] || [ "${1:-}" = "opal" ]; then
    echo "OPAL mode (optional CAMI profiler assessment)."
    PIP_OPTS=()
    if [ "${SAMOVAR_OFFLINE:-0}" != "0" ]; then
        PIP_OPTS+=(--offline --no-build-isolation)
        if [ -n "${SAMOVAR_WHEELHOUSE:-}" ]; then
            PIP_OPTS+=(--no-index --find-links "$SAMOVAR_WHEELHOUSE")
        fi
    fi
    install_opal
    exit $?
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

# Write ~/.config/samovar/config.json (and build/config.json copy)
export USER_CFG_DIR
export SAMOVAR_ROOT="$ROOT"
export PYTHON_PATH
"$PYTHON_PATH" - <<'PY'
import os, shutil
from pathlib import Path
from samovar.paths import (
    PACKAGE_VERSION,
    collect_runtime_path_dirs,
    discover_opal,
    discover_tools,
    load_config,
    write_config,
)

root = os.environ["SAMOVAR_ROOT"]
python_path = os.environ["PYTHON_PATH"]
existing = load_config()
tools = dict(discover_tools())
for name, path in (existing.get("tools") or {}).items():
    if str(path or "").strip():
        tools[name] = path
path_extra = existing.get("path") or existing.get("extra_path") or []
tool_envs = existing.get("tool_envs") or {}
payload = {
    "version": PACKAGE_VERSION,
    "root": root,
    "python_path": python_path,
    "r_path": shutil.which("R") or "",
    "r_lib_path": existing.get("r_lib_path", ""),
    "iss_path": shutil.which("iss") or tools.get("iss", ""),
    "opal_path": existing.get("opal_path") or discover_opal() or shutil.which("opal.py") or shutil.which("opal") or "",
    "ncbi_email": os.environ.get("NCBI_EMAIL", ""),
    "test_genomes": os.path.join(root, "data", "test_genomes"),
    "genomes": existing.get("genomes") or str(Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "samovar" / "genomes"),
    "processed_genomes": existing.get("processed_genomes")
        or existing.get("genomes")
        or str(Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "samovar" / "genomes"),
    "genome_dirs": existing.get("genome_dirs") or [],
    "tools": tools,
}
if path_extra:
    payload["path"] = path_extra
if tool_envs:
    payload["tool_envs"] = tool_envs
if existing.get("annotation_regenerate_r"):
    payload["annotation_regenerate_r"] = existing["annotation_regenerate_r"]
if payload.get("opal_path"):
    tools.setdefault("opal.py", payload["opal_path"])
    payload["tools"] = tools
write_config(payload)
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
env_path.write_text(
    f'export NCBI_EMAIL="{email}"\n'
    f'export SAMOVAR_ROOT="{root}"\n'
    f'export PYTHON_PATH="{python_path}"\n'
    f'export PATH="{path_export}"\n'
    f"{conda_lib}",
    encoding="utf-8",
)
print("Wrote user + repo config")
print(f"Wrote {env_path} PATH with {len(dirs)} tool bin dir(s)")
genomes = payload.get("genomes") or ""
processed = payload.get("processed_genomes") or genomes
lib_dirs = payload.get("genome_dirs") or []
print("")
print("Genome cache (NCBI / user assemblies, reused by `samovar prepare`):")
print(f"  Downloaded genomes:  {genomes}")
print(f"  Processed genomes:   {processed}")
print(f"  Extra libraries:     {', '.join(lib_dirs) if lib_dirs else '(none yet; add dirs to genome_dirs in config.json)'}")
print("")
print("IMPORTANT: data/test_genomes contains TRUNCATED assemblies for ISS and CI only.")
print("They are never used as NCBI substitutes. Real metagenomes reuse only genomes")
print("previously downloaded by SamovaR or listed in genome_dirs.")
print("  samovar prepare --reuse-genomes      # default: symlink from cache if found")
print("  samovar prepare --no-reuse-genomes   # download into the genomes cache")
print("  samovar prepare --genome-dirs DIR[:DIR]")
print("  samovar prepare --test-genomes       # allow truncated stubs (ISS/CI only)")
print("Prefer symlinks into $out_dir/genomes; copy with a warning if linking fails.")
if payload.get("opal_path"):
    print(f"OPAL (cami-opal): {payload['opal_path']}")
else:
    print("OPAL: not found. ./install.sh OPAL to install CAMI OPAL and record opal_path.")
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
if [ "$SMOKE_FAIL" != "0" ]; then
    echo "Smoke test reported missing pieces (see above). Core Python package may still work."
else
    echo "Smoke test passed."
fi

echo "Installation completed successfully!"
echo "Config: $USER_CFG_DIR/config.json"
if [ -f "$ROOT/build/config.json" ]; then
    echo "Also:   $ROOT/build/config.json"
fi
echo "Use NCBI_EMAIL=$NCBI_EMAIL (export it in new shells, or: source $USER_CFG_DIR/env)"
echo "Shared/HPC install without PATH edits: SAMOVAR_UPDATE_SHELL=0 ./install.sh"
