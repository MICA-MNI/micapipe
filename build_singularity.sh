#!/bin/bash
#
# Singularity build from local Docker image
#

set -e

DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

BASE_DIR="/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"

export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"
export TMPDIR="${BASE_DIR}/.tmp"

log() { echo "[$(date +%H:%M:%S)] $*"; }

log "Building SIF from local Docker: $FULL_DOCKER_IMAGE"
log "Output: $OUTPUT_PATH"

mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "${BASE_DIR}/.tmp"

docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null || { log "Image not found"; exit 1; }

rm -f "$OUTPUT_PATH"

singularity build --force "$OUTPUT_PATH" "docker-daemon://${FULL_DOCKER_IMAGE}"

log "Done: $(du -h "$OUTPUT_PATH" | cut -f1)"
