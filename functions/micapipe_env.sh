#!/bin/bash
# Run a command inside one pinned tool environment, scoped to a subshell.
micapipe_env() (
    local name="$1"; shift
    local prefix="${MICAPIPE_ENV_ROOT:-/opt/environments}/${name}"
    if [[ ! -x "${prefix}/bin/python" ]]; then
        echo "Missing MICAPIPE environment: ${prefix}" >&2
        exit 127
    fi
    unset PYTHONHOME PYTHONPATH
    export PATH="${prefix}/bin:${PATH}"
    export VIRTUAL_ENV="${prefix}"
    local libs="${prefix}/lib"
    # Keep the host NVIDIA driver libraries that singularity --nv provides
    if [[ -d /.singularity.d/libs ]]; then libs="${libs}:/.singularity.d/libs"; fi
    export LD_LIBRARY_PATH="${libs}"
    exec "$@"
)
