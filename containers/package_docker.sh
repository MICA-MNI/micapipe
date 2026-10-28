#!/bin/bash
# Pack the prepared filesystem as a Docker image archive, without a Docker daemon.
# Load it with: docker load -i micapipe-VERSION.docker.tar.gz
set -euo pipefail
root="${1:?workspace}"
version="${2:?container version}"
variant="${3:?nocuda or cuda}"
root=$(cd "$root" && pwd)
target="$root/artifacts/micapipe-${version}.docker.tar.gz"
[[ ! -e "$target" ]] || { echo "Refusing to overwrite $target" >&2; exit 1; }
mkdir -p "$root/artifacts" "$root/tmp"
stage=$(mktemp -d "$root/tmp/docker-archive.XXXXXX")
trap 'rm -rf "$stage" "$target.part"' EXIT
# One layer holding the whole filesystem, owned by root like the SIF
tar --numeric-owner --owner=0 --group=0 --exclude=./.singularity.d --exclude=./.dockerenv \
    -C "$root/runtime" -cf "$stage/layer.tar" .
read -r config layer < <(python3 - "$stage" "$version" "$variant" <<'ARCHIVE'
import hashlib, json, sys
from pathlib import Path
stage, version, variant = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
digest = hashlib.sha256()
with open(stage / "layer.tar", "rb") as layer:
    for block in iter(lambda: layer.read(1 << 24), b""):
        digest.update(block)
layer_id = digest.hexdigest()
(stage / layer_id).mkdir()
(stage / "layer.tar").rename(stage / layer_id / "layer.tar")
env = ["PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
       "LANG=en_US.UTF-8", "LC_ALL=en_US.UTF-8"]
if variant == "cuda":
    env += ["NVIDIA_VISIBLE_DEVICES=all", "NVIDIA_DRIVER_CAPABILITIES=compute,utility"]
config = {
    "architecture": "amd64",
    "os": "linux",
    "config": {"Env": env, "Entrypoint": ["/opt/micapipe/containers/docker_entrypoint.sh"],
               "Labels": {"org.opencontainers.image.version": version}},
    "rootfs": {"type": "layers", "diff_ids": [f"sha256:{layer_id}"]},
    "history": [{"created_by": f"micapipe containers/package_docker.sh {version}"}],
}
raw = json.dumps(config, sort_keys=True).encode()
config_name = hashlib.sha256(raw).hexdigest() + ".json"
(stage / config_name).write_bytes(raw)
manifest = [{"Config": config_name, "RepoTags": [f"micapipe:{version}"], "Layers": [f"{layer_id}/layer.tar"]}]
(stage / "manifest.json").write_text(json.dumps(manifest))
print(config_name, f"{layer_id}/layer.tar")
ARCHIVE
)
(cd "$stage" && tar --numeric-owner --owner=0 --group=0 -cf - manifest.json "$config" "$layer") \
    | pigz -p "${MICAPIPE_BUILD_CPUS:-8}" > "$target.part"
mv "$target.part" "$target"
chmod a+r "$target"
sha256sum "$target" > "$root/evidence/container-${version}.docker.sha256"
printf '%s\n' "$target"
