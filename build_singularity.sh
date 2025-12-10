#!/bin/bash
#
# Singularity build from local Docker image
# Uses docklog "🔧 Building SIF directly from Docker daemon..."r-daemon:// protocol - Singularity reads directly from Docker daemon
# No intermediate tar files, no temp space issues
#

set -e

DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

BASE_DIR="/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"

# Use export03 for ALL temp files to avoid /export01 space issues
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"
export TMPDIR="${BASE_DIR}/.tmp"
export TEMP="${BASE_DIR}/.tmp"
export TMP="${BASE_DIR}/.tmp"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================================
# Pre-flight checks
# ============================================================================
log "🚀 SINGULARITY BUILD (docker-daemon:// method)"
log "📦 Image: $FULL_DOCKER_IMAGE"
log "📍 Output: $OUTPUT_PATH"

# Check local Docker image exists
if ! docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null; then
    log "❌ Local Docker image not found: $FULL_DOCKER_IMAGE"
    log "   Available images:"
    docker images | grep micapipe || echo "   No micapipe images found"
    exit 1
fi

LOCAL_SIZE=$(docker image inspect "$FULL_DOCKER_IMAGE" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')
log "✅ Found local image: $LOCAL_SIZE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "${BASE_DIR}/.tmp"

# Check space
AVAILABLE=$(df -BG "$BASE_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$AVAILABLE" -lt 200 ]; then
    log "❌ Need 200GB+ space, only ${AVAILABLE}GB available"
    exit 1
fi
log "✅ Space check: ${AVAILABLE}GB available"

# Remove existing output
if [ -f "$OUTPUT_PATH" ]; then
    log "⚠️  Removing existing SIF file"
    rm -f "$OUTPUT_PATH"
fi

START_TIME=$(date +%s)

# ============================================================================
# Build using docker-daemon:// protocol
# Singularity reads directly from Docker daemon - no intermediate files needed
# This avoids all temp space issues on /export01
# ============================================================================

log "� Building SIF directly from Docker daemon..."
log "   Using: docker-daemon://${FULL_DOCKER_IMAGE}"

# Monitor SIF creation in background
(
    sleep 10
    while true; do
        if [ -f "$OUTPUT_PATH" ]; then
            SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")
            log "📈 SIF build progress: $SIZE"
        else
            log "⏳ SIF build in progress..."
        fi
        sleep 30
        # Check if parent process is still running
        if ! kill -0 $$ 2>/dev/null; then
            break
        fi
    done
) &
MONITOR_PID=$!

# Build directly from Docker daemon - no tar, no temp files
singularity build --force "$OUTPUT_PATH" "docker-daemon://${FULL_DOCKER_IMAGE}"
SINGULARITY_EXIT=$?

# Stop monitor
kill $MONITOR_PID 2>/dev/null || true
wait $MONITOR_PID 2>/dev/null || true

if [ $SINGULARITY_EXIT -ne 0 ]; then
    log "❌ Singularity build failed with exit code $SINGULARITY_EXIT"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))
DURATION_SEC=$((DURATION % 60))
SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)

log "============================================="
log "✅ BUILD COMPLETE"
log "============================================="
log "📦 File: $OUTPUT_PATH"
log "📊 Size: $SIZE"
log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"
log ""
log "🧪 Test with:"
log "   singularity run $OUTPUT_PATH --help"
log ""
log "🚀 Ready for deployment!"

# Quick verification
if [ -f "$OUTPUT_PATH" ] && [ -s "$OUTPUT_PATH" ]; then
    log "✅ SIF file created successfully"
else
    log "❌ ERROR: SIF file not created or empty"
    exit 1
fi