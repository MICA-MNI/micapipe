#!/bin/bash
# Build on Linux. A normal Docker container installs into a data-disk prefix.
set -euo pipefail
root="${1:?Usage: build_beta_server.sh WORKSPACE}"
root=$(cd "$root" && pwd)
[[ $(uname -s) == Linux ]] || { echo "Build on the Linux server." >&2; exit 1; }
mkdir -p "$root"/{runtime/opt,evidence,artifacts,tmp/home,cache}
(cd "$root" && sha256sum -c source/containers/SOURCES.sha256)
python3 "$root/source/containers/stage_lamareg_models.py"
docker build -t micapipe-v1-beta-os:20261022 "$root/source/containers/os"
snapshot=$(mktemp "$root/tmp/install-runtime.XXXXXX")
cp "$root/source/containers/install_runtime.sh" "$snapshot"
docker run --rm --cpus "${MICAPIPE_BUILD_CPUS:-12}" --memory "${MICAPIPE_BUILD_MEMORY:-48g}" \
  --user "$(id -u):$(id -g)" -e HOME=/build/tmp/home \
  --mount "type=bind,src=$root,dst=/build" \
  --mount "type=bind,src=$root/runtime/opt,dst=/opt" \
  micapipe-v1-beta-os:20261022 bash "/build/tmp/${snapshot##*/}"
cid=$(docker create micapipe-v1-beta-os:20261022)
trap 'docker rm "$cid" >/dev/null' EXIT
docker export "$cid" | tar --no-same-owner -xf - -C "$root/runtime"
bash "$root/source/containers/prepare_singularity_metadata.sh" "$root"
mkdir -p "$root/runtime/.singularity.d/env"
cp "$root/source/containers/runtime_env.sh" "$root/runtime/.singularity.d/env/90-environment.sh"
printf '#!/bin/bash\nexec /opt/micapipe/micapipe "$@"\n' > "$root/runtime/.singularity.d/runscript"
chmod +x "$root/runtime/.singularity.d/runscript"
python3 - "$root" <<'BUILDINFO'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
commit=subprocess.check_output(['git','-C',str(root/'source'),'rev-parse','HEAD'],text=True).strip()
version=json.loads((root/'source/containers/versions.json').read_text())['candidate']
info={'container_version':version,'source_commit':commit,'fix_default':'legacy','scientific_equivalence_established':False}
(root/'runtime/opt/micapipe/BUILDINFO.json').write_text(json.dumps(info,indent=2)+'\n')
(root/'runtime/.singularity.d/labels.json').write_text(json.dumps({'org.opencontainers.image.version':version.lstrip('v'),'org.opencontainers.image.revision':commit}))
BUILDINFO
bash "$root/source/containers/package_beta.sh" "$root"
