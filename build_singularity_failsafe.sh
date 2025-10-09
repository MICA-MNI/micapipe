#!/bin/bash
#
# FAILSAFE Singularity build - most reliable method
# Uses separate steps to avoid pipe issues
#

set -e

DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"
TAR_FILE="${BASE_DIR}/micapipe_docker.tar"

# Conservative settings to avoid issues
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"
unset SINGULARITY_MEMORY  # Let singularity manage memory

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================================
# Pre-flight checks
# ============================================================================
log "🛡️  FAILSAFE SINGULARITY BUILD"
log "📦 Image: $FULL_DOCKER_IMAGE"
log "📍 Output: $OUTPUT_PATH"

# Check local Docker image exists
if ! docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null; then
    log "❌ Local Docker image not found: $FULL_DOCKER_IMAGE"
    exit 1
fi

LOCAL_SIZE=$(docker image inspect "$FULL_DOCKER_IMAGE" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')
log "✅ Found local image: $LOCAL_SIZE"

# Create directories
mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Check space (need 2x image size for tar + sif)
AVAILABLE=$(df -BG "$BASE_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$AVAILABLE" -lt 300 ]; then
    log "❌ Need 300GB+ space for failsafe method, only ${AVAILABLE}GB available"
    exit 1
fi
log "✅ Space check: ${AVAILABLE}GB available"

# Clean up any previous attempts
log "🧹 Cleaning up previous attempts..."
rm -f "$OUTPUT_PATH" "$TAR_FILE"

START_TIME=$(date +%s)

# ============================================================================
# Step 1: Export Docker to tar (with progress monitoring)
# ============================================================================
log "📤 Step 1: Exporting Docker image to tar..."
log "📍 Tar file: $TAR_FILE"

# Start size monitoring in background
(
    while true; do
        if [ -f "$TAR_FILE" ]; then
            SIZE=$(du -h "$TAR_FILE" 2>/dev/null | cut -f1 || echo "0")
            log "📈 Tar size: $SIZE"
        else
            log "⏳ Waiting for tar file to appear..."
        fi
        sleep 30
    done
) &
TAR_MONITOR_PID=$!

# Export with timeout
log "🔄 Running: docker save '$FULL_DOCKER_IMAGE' -o '$TAR_FILE'"
if timeout 3600 docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"; then
    kill $TAR_MONITOR_PID 2>/dev/null || true
    wait $TAR_MONITOR_PID 2>/dev/null || true
    
    TAR_SIZE=$(du -h "$TAR_FILE" | cut -f1)
    log "✅ Docker export complete: $TAR_SIZE"
else
    kill $TAR_MONITOR_PID 2>/dev/null || true
    wait $TAR_MONITOR_PID 2>/dev/null || true
    log "❌ Docker export failed or timed out"
    exit 1
fi

# ============================================================================
# Step 2: Build SIF from tar (with progress monitoring)
# ============================================================================
log "🔧 Step 2: Building SIF from tar..."
log "📍 SIF file: $OUTPUT_PATH"

# Start size monitoring in background
(
    while true; do
        if [ -f "$OUTPUT_PATH" ]; then
            SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")
            log "📈 SIF size: $SIZE"
        else
            log "⏳ Waiting for SIF file to appear..."
        fi
        sleep 30
    done
) &
SIF_MONITOR_PID=$!

# Build with timeout
log "🔄 Running: singularity build --force '$OUTPUT_PATH' docker-archive://'$TAR_FILE'"
if timeout 3600 singularity build --force "$OUTPUT_PATH" "docker-archive://$TAR_FILE"; then
    kill $SIF_MONITOR_PID 2>/dev/null || true
    wait $SIF_MONITOR_PID 2>/dev/null || true
    
    SIF_SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)
    log "✅ SIF build complete: $SIF_SIZE"
else
    kill $SIF_MONITOR_PID 2>/dev/null || true
    wait $SIF_MONITOR_PID 2>/dev/null || true
    log "❌ SIF build failed or timed out"
    exit 1
fi

# ============================================================================
# Step 3: Cleanup and verify
# ============================================================================
log "🧹 Cleaning up tar file..."
rm -f "$TAR_FILE"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))
DURATION_SEC=$((DURATION % 60))

log "============================================="
log "✅ FAILSAFE BUILD COMPLETE"
log "============================================="
log "📦 File: $OUTPUT_PATH"
log "📊 Size: $(du -h "$OUTPUT_PATH" | cut -f1)"
log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"
log "🎯 Method: failsafe (tar-based)"
log ""
log "🧪 Test with:"
log "   singularity run $OUTPUT_PATH --help"
log ""
log "🚀 Ready for deployment!"

# Final verification
if [ -f "$OUTPUT_PATH" ] && [ -s "$OUTPUT_PATH" ]; then
    log "✅ SIF file created successfully"
else
    log "❌ ERROR: SIF file not created or empty"
    exit 1
fi