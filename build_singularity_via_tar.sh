#!/usr/bin/env bash

# Backwards-compatible wrapper that exports the micapipe Docker image to a tar
# archive and then converts it to a SIF via build_singularity_auto.sh.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

export CACHE_ROOT="${CACHE_ROOT:-/host/cassio/export03/data/enning}"
exec "${SCRIPT_DIR}/build_singularity_auto.sh" --source docker-archive "$@"
