#!/bin/bash
# Pack the prepared filesystem directly; avoid a second full tool-tree copy.
set -euo pipefail
root="${1:?workspace}"
root=$(cd "$root" && pwd)
version=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["candidate"])' "$root/source/containers/versions.json")
target="$root/artifacts/micapipe-${version}.sif"
[[ ! -e "$target" ]] || { echo "Refusing to overwrite $target" >&2; exit 1; }
mkdir -p "$root/artifacts" "$root/tmp"
squash=$(mktemp "$root/tmp/micapipe-rootfs.XXXXXX")
candidate="$root/artifacts/.micapipe-beta-build-$$.sif"
trap 'rm -f "$squash" "$candidate"' EXIT
mksquashfs "$root/runtime" "$squash" -noappend -all-root \
    -processors "${MICAPIPE_BUILD_CPUS:-8}" -mem 4G -comp gzip -Xcompression-level 1 -no-progress
singularity sif new "$candidate"
singularity sif add "$candidate" "$squash" --datatype 4 --parttype 2 --partfs 1 --partarch 2
singularity sif add "$candidate" "$root/runtime/.singularity.d/labels.json" --datatype 3
singularity exec --cleanenv "$candidate" /bin/true
probe=$(singularity exec --cleanenv --env MICAPIPE_FIX_BACKEND=python "$candidate" printenv MICAPIPE_FIX_BACKEND)
[[ "$probe" == python ]] || { echo "Container environment override failed" >&2; exit 1; }
chmod a+r "$candidate"
mv "$candidate" "$target"
sha256sum "$target" > "$root/evidence/container-${version}.sha256"
cp "$root/evidence/container-${version}.sha256" "$root/evidence/container.sha256"
singularity inspect --json "$target" > "$root/evidence/container-inspect.json"
printf '%s\n' "$target"
