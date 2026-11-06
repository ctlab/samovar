#!/usr/bin/env bash
# Install only the portable bundle in Agents/, never the checkout's .agents/skills.
set -euo pipefail
bundle_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "${SAMOVAR_HARNESS_PYTHON:-}" ]; then
    harness_python="$SAMOVAR_HARNESS_PYTHON"
elif command -v python3 >/dev/null 2>&1; then
    harness_python="$(command -v python3)"
else
    harness_python="$(command -v python)"
fi
exec "$harness_python" "$bundle_root/Agents/install.py" "$@"
