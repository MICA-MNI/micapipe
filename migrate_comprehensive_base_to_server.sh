#!/bin/bash
# Sync a curated MICApipe image-build workspace to the remote experiment server.

set -euo pipefail

SOURCE_DIR="${MICAPIPE_SOURCE_DIR:-$PWD}"
LOGIN_HOST="${MICAPIPE_LOGIN_HOST:-login.bic.mni.mcgill.ca}"
COMPUTE_HOST="${MICAPIPE_COMPUTE_HOST:-bb-compxg-01}"
REMOTE_BASE="${MICAPIPE_REMOTE_BASE:-/export03/data/enning}"
LAB_NAME="${MICAPIPE_LAB_NAME:-micapipe-image-lab-$(date +%Y%m%d-%H%M%S)}"

usage() {
    cat <<'EOF'
Usage: ./migrate_comprehensive_base_to_server.sh [options]

Options:
  --source DIR         Local source repository (default: current directory)
  --login-host HOST    Jump host (default: login.bic.mni.mcgill.ca)
  --compute-host HOST  Target compute host (default: bb-compxg-01)
  --remote-base DIR    Base directory on compute host (default: /export03/data/enning)
  --lab-name NAME      Remote experiment folder name
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --source)
            SOURCE_DIR="$2"
            shift 2
            ;;
        --login-host)
            LOGIN_HOST="$2"
            shift 2
            ;;
        --compute-host)
            COMPUTE_HOST="$2"
            shift 2
            ;;
        --remote-base)
            REMOTE_BASE="$2"
            shift 2
            ;;
        --lab-name)
            LAB_NAME="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! -d "$SOURCE_DIR" ]]; then
    echo "Source directory not found: $SOURCE_DIR" >&2
    exit 1
fi

REMOTE_DIR="${REMOTE_BASE}/${LAB_NAME}"
REMOTE_REPO_DIR="${REMOTE_DIR}/repo"

read -r -d '' REMOTE_SETUP <<"EOF" || true
base="$1"
name="$2"
dir="${base}/${name}"
if [ -e "$dir" ]; then
    dir="${base}/${name}-$(date +%H%M%S)"
fi
mkdir -p "$dir/repo" "$dir/logs" "$dir/staging"
printf '%s\n' "$dir"
EOF

echo "🚚 Syncing MICApipe image workspace to server"
echo "📍 Source: ${SOURCE_DIR}"
echo "🌉 Jump host: ${LOGIN_HOST}"
echo "🖥️  Compute host: ${COMPUTE_HOST}"
echo "📁 Remote base: ${REMOTE_BASE}"
echo

ACTUAL_REMOTE_DIR="$(
    ssh "$LOGIN_HOST" \
        "ssh $COMPUTE_HOST 'bash -s' -- '$REMOTE_BASE' '$LAB_NAME'" <<<"$REMOTE_SETUP"
)"
ACTUAL_REMOTE_DIR="$(printf '%s\n' "$ACTUAL_REMOTE_DIR" | tail -n 1)"
REMOTE_REPO_DIR="${ACTUAL_REMOTE_DIR}/repo"

SYNC_ITEMS=(
    Dockerfile.mamba-base
    Dockerfile.main
    build_comprehensive_base_server.sh
    build_main_image_server.sh
    prepare_build_context.sh
    migrate_comprehensive_base_to_server.sh
    check_docker_space.sh
    .dockerignore
    .gitignore
    fix_settings.sh
    micapipe
    micapipe.py
    pyproject.toml
    micapipe_environment.yml
    micapipe_environment_base.yml
    readme.md
    R_config
    parcellations
    surfaces
    MNI152Volumes
    MICs60_T1-atlas
    fsl_conf
    functions
)

(
    cd "$SOURCE_DIR"
    COPYFILE_DISABLE=1 tar -czf - "${SYNC_ITEMS[@]}"
) | ssh "$LOGIN_HOST" "ssh $COMPUTE_HOST 'tar -xzf - -C \"$REMOTE_REPO_DIR\"'"

ssh "$LOGIN_HOST" "ssh $COMPUTE_HOST 'chmod +x \"$REMOTE_REPO_DIR\"/*.sh'"

echo
echo "✅ Remote experiment workspace ready"
echo "📁 Remote repo: ${REMOTE_REPO_DIR}"
echo "🚀 Next steps on the server:"
echo "   cd ${REMOTE_REPO_DIR}"
echo "   ./build_comprehensive_base_server.sh"
echo "   ./build_main_image_server.sh"
