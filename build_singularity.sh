#!/bin/bash
#
# Singularity build from a local Docker image.
#

set -euo pipefail

DOCKER_IMAGE="${MICAPIPE_DOCKER_IMAGE:-ghcr.io/mica-mni/micapipe}"
DOCKER_TAG="${MICAPIPE_DOCKER_TAG:-latest}"
BASE_DIR="${MICAPIPE_SIF_BASE_DIR:-/export03/data/enning}"
OUTPUT_DIR="${MICAPIPE_SIF_OUTPUT_DIR:-${BASE_DIR}/singularity}"
OUTPUT_PATH="${MICAPIPE_SIF_OUTPUT_PATH:-${OUTPUT_DIR}/micapipe_v1_beta.sif}"

usage() {
    cat <<'EOF'
Usage: ./build_singularity.sh [options]

Options:
  --tag TAG           Docker tag to convert (default: latest)
  --image IMAGE       Docker image repository (default: ghcr.io/mica-mni/micapipe)
  --base-dir DIR      Base directory for Singularity cache/tmp/output
  --output PATH       Full output path for the resulting SIF
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)
            DOCKER_TAG="$2"
            shift 2
            ;;
        --image)
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --base-dir)
            BASE_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_PATH="$2"
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

if [[ -z "${MICAPIPE_SIF_OUTPUT_PATH:-}" ]]; then
    OUTPUT_DIR="$(dirname "$OUTPUT_PATH")"
fi

FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

export SINGULARITY_CACHEDIR="${MICAPIPE_SINGULARITY_CACHEDIR:-${BASE_DIR}/.singularity_cache}"
export SINGULARITY_TMPDIR="${MICAPIPE_SINGULARITY_TMPDIR:-${BASE_DIR}/.singularity_tmp}"
export TMPDIR="${MICAPIPE_TMPDIR:-${BASE_DIR}/.tmp}"

log() { echo "[$(date +%H:%M:%S)] $*"; }

log "Building SIF from local Docker: $FULL_DOCKER_IMAGE"
log "Output: $OUTPUT_PATH"
log "Cache dir: $SINGULARITY_CACHEDIR"
log "Temp dir: $SINGULARITY_TMPDIR"

mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$TMPDIR"

if ! command -v singularity >/dev/null 2>&1; then
    log "singularity command not found"
    exit 1
fi

docker image inspect "$FULL_DOCKER_IMAGE" >/dev/null 2>&1 || {
    log "Image not found: $FULL_DOCKER_IMAGE"
    exit 1
}

rm -f "$OUTPUT_PATH"

# Pick a transport. Both docker-daemon:// and docker-archive:// stage the
# image's layers in the Docker daemon's data-root tmp ($DOCKER_ROOT/tmp, e.g.
# /export01/docker/tmp on this server) before writing the SIF/tar, so a tight
# data-root partition fails EITHER way with 'no space left on device'. There is
# no transport that sidesteps this — the only cures are freeing space on that
# partition or moving Docker's data-root to a larger volume. So when the root is
# tight we reclaim what we safely can and, if still short, fail fast with an
# actionable message rather than crashing minutes into a doomed export.
#
# Override with MICAPIPE_SIF_TRANSPORT=docker-archive or =docker-daemon.
TRANSPORT="${MICAPIPE_SIF_TRANSPORT:-auto}"
free_gb() { df -BG --output=avail "$1" 2>/dev/null | tail -1 | tr -dc 0-9 || echo 0; }

if [[ "$TRANSPORT" == "auto" ]]; then
    DOCKER_ROOT="$(docker info --format '{{.DockerRootDir}}' 2>/dev/null || echo "/var/lib/docker")"
    AVAIL_GB="$(free_gb "$DOCKER_ROOT")"
    IMAGE_GB="$(docker image inspect "$FULL_DOCKER_IMAGE" --format '{{.Size}}' | awk '{printf "%d", $1/1024/1024/1024 + 1}')"
    NEEDED_GB=$((IMAGE_GB + 5))  # 5GB buffer
    log "Docker root: $DOCKER_ROOT (${AVAIL_GB}G free); image: ~${IMAGE_GB}G"

    if (( AVAIL_GB < NEEDED_GB )); then
        log "Docker root tight (need ~${NEEDED_GB}G, have ${AVAIL_GB}G) — reclaiming safe space"
        # Safe to drop: build cache, dangling (untagged) images, and stale CI
        # images from previous runs (ci-* tags other than this build's).
        docker builder prune -f >/dev/null 2>&1 || true
        docker image prune -f   >/dev/null 2>&1 || true
        docker images "$DOCKER_IMAGE" --format '{{.Repository}}:{{.Tag}}' 2>/dev/null \
            | grep ':ci-' | grep -v ":${DOCKER_TAG}\$" \
            | xargs -r docker rmi >/dev/null 2>&1 || true
        AVAIL_GB="$(free_gb "$DOCKER_ROOT")"
        log "After reclaim: ${AVAIL_GB}G free"
    fi

    if (( AVAIL_GB >= NEEDED_GB )); then
        TRANSPORT="docker-daemon"
    else
        log "ERROR: $DOCKER_ROOT has ${AVAIL_GB}G free but the SIF build needs ~${NEEDED_GB}G there."
        log "Both transports stage layers in \$DOCKER_ROOT/tmp, so this cannot succeed as-is."
        log "Fix: free that partition (e.g. 'docker image prune -a', or remove old image tags)"
        log "     or move Docker's data-root to a larger volume, then retry."
        exit 1
    fi
fi

case "$TRANSPORT" in
    docker-daemon)
        log "Transport: docker-daemon:// (streaming)"
        singularity build --force "$OUTPUT_PATH" "docker-daemon://${FULL_DOCKER_IMAGE}"
        ;;
    docker-archive)
        TAR="${TMPDIR}/$(basename "$OUTPUT_PATH" .sif).tar"
        log "Transport: docker-archive:// via $TAR"
        rm -f "$TAR"
        trap 'rm -f "$TAR"' EXIT
        docker save -o "$TAR" "$FULL_DOCKER_IMAGE"
        singularity build --force "$OUTPUT_PATH" "docker-archive://${TAR}"
        rm -f "$TAR"
        ;;
    *)
        log "Unknown MICAPIPE_SIF_TRANSPORT: $TRANSPORT (use auto, docker-daemon, or docker-archive)"
        exit 1
        ;;
esac

log "Done: $(du -h "$OUTPUT_PATH" | cut -f1)"
