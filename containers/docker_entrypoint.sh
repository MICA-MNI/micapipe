#!/bin/bash
# Docker entrypoint: load the environment Singularity reads from .singularity.d/env.
source /opt/micapipe/containers/runtime_env.sh
exec /opt/micapipe/micapipe "$@"
