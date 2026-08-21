#!/usr/bin/env bash
# Shared helpers for example pipelines.
# Source from examples/*/pipeline.sh after setting SCRIPT_DIR.
#
# Overridable environment:
#   SAMOVAR              repo root (otherwise inferred from this file)
#   CONDA_PREFIX         conda env; its bin/ is prepended when set
#   NCBI_EMAIL / ENTREZ_EMAIL / SAMOVAR_EMAIL
#   SAMOVAR_KAIJU_DB, SAMOVAR_KRAKEN2_DB, SAMOVAR_KRAKEN_DB, SAMOVAR_KRAKENUNIQ_DB
#     optional local indexes (skip download / rebuild)
#   KAIJU_INDEX_URL, KRAKEN2_INDEX_URL, MINIKRAKEN_URL
#     public archives used when a local index is not provided

if [[ -z "${SAMOVAR:-}" || ! -d "${SAMOVAR}/src/samovar" ]]; then
  SAMOVAR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

samovar_setup_env() {
  local conda_bin=""
  if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/bin" ]]; then
    conda_bin="${CONDA_PREFIX}/bin"
  fi
  export PATH="${SAMOVAR}/bin${conda_bin:+:${conda_bin}}:${PATH}"
  export PYTHONPATH="${SAMOVAR}/src${PYTHONPATH:+:${PYTHONPATH}}"
  export NCBI_EMAIL="${NCBI_EMAIL:-${ENTREZ_EMAIL:-${SAMOVAR_EMAIL:-anonymous@example.com}}}"
}

# wget + extract. Skip if dest_dir already contains marker (filename).
samovar_fetch_archive() {
  local url="$1"
  local dest_dir="$2"
  local marker="${3:-}"
  mkdir -p "$dest_dir"
  if [[ -n "$marker" ]] && find -L "$dest_dir" -name "$marker" | grep -q .; then
    echo "Already present under ${dest_dir}: ${marker}"
    return 0
  fi
  local tmp
  tmp="$(mktemp "${dest_dir}/.download.XXXXXX.tgz")"
  wget -O "$tmp" "$url"
  tar -xzf "$tmp" -C "$dest_dir"
  rm -f "$tmp"
  if [[ -n "$marker" && ! -e "${dest_dir}/${marker}" ]]; then
    local nested
    nested="$(find "$dest_dir" -name "$marker" -print -quit)"
    if [[ -n "$nested" ]]; then
      local nested_dir
      nested_dir="$(dirname "$nested")"
      if [[ "$nested_dir" != "$dest_dir" ]]; then
        shopt -s dotglob nullglob
        mv "${nested_dir}"/* "$dest_dir"/ || true
        shopt -u dotglob nullglob
      fi
    fi
  fi
}

# Use an existing dir (env) or download a public index into dest_dir.
samovar_use_or_fetch_db() {
  local env_path="$1"
  local dest_dir="$2"
  local url="$3"
  local marker="$4"
  if [[ -n "$env_path" && -e "$env_path" ]]; then
    mkdir -p "$(dirname "$dest_dir")"
    if [[ -d "$env_path" ]]; then
      ln -sfn "$env_path" "$dest_dir"
    else
      mkdir -p "$dest_dir"
      ln -sfn "$env_path" "$dest_dir/$(basename "$env_path")"
      local parent
      parent="$(dirname "$env_path")"
      for extra in nodes.dmp names.dmp; do
        if [[ -f "${parent}/${extra}" ]]; then
          ln -sfn "${parent}/${extra}" "${dest_dir}/${extra}"
        fi
      done
    fi
    echo "Using local database: ${env_path} -> ${dest_dir}"
    return 0
  fi
  samovar_fetch_archive "$url" "$dest_dir" "$marker"
}
