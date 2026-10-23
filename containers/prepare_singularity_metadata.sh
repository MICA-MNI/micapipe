#!/bin/bash
# Retain Singularity's native startup hooks, including --env injection.
set -euo pipefail
root="${1:?workspace}"
base="${2:-}"
if [[ -z "$base" ]]; then
    stage=$(mktemp -d "$root/tmp/singularity-metadata.XXXXXX")
    base="$stage/rootfs"
    SINGULARITY_TMPDIR="$root/tmp" SINGULARITY_CACHEDIR="$root/cache/singularity" \
        singularity build --sandbox "$base" docker-daemon://micapipe-v1-beta-os:20261022
fi
test -d "$base/.singularity.d/env"
mkdir -p "$root/runtime/.singularity.d"
cp -a "$base/.singularity.d/." "$root/runtime/.singularity.d/"
cp "$root/source/containers/runtime_env.sh" "$root/runtime/.singularity.d/env/90-environment.sh"
printf '#!/bin/bash\nexec /opt/micapipe/micapipe "$@"\n' > "$root/runtime/.singularity.d/runscript"
chmod +x "$root/runtime/.singularity.d/runscript"
