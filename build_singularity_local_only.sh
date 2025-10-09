#!/bin/bash
#
# Local-only fast Singularity build - no registry access needed
# Optimized for local Docker images with 128GB RAM server
#

set -e

DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"

# Performance settings for your 128GB server
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"
export SINGULARITY_MEMORY="64G"  # Use half your RAM
export OMP_NUM_THREADS=$(nproc)

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================================
# Pre-flight checks
# ============================================================================
log "🚀 LOCAL FAST SINGULARITY BUILD"
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
mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

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
# Method 1: Streaming (fastest - no intermediate files)
# ============================================================================
log "⚡ Trying streaming method (fastest)..."

if docker save "$FULL_DOCKER_IMAGE" | singularity build --force "$OUTPUT_PATH" docker-archive:///dev/stdin 2>/dev/null; then
    log "✅ Streaming method succeeded!"
    USED_METHOD="streaming"
else
    log "⚠️  Streaming failed, trying tar method..."
    
    # ============================================================================
    # Method 2: Tar method (more reliable)
    # ============================================================================
    TAR_FILE="${BASE_DIR}/micapipe_docker_$$.tar"
    
    log "📤 Exporting Docker to tar..."
    docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"
    
    TAR_SIZE=$(du -h "$TAR_FILE" | cut -f1)
    log "✅ Export complete: $TAR_SIZE"
    
    log "🔧 Building SIF from tar..."
    singularity build --force \
        "$OUTPUT_PATH" \
        "docker-archive://$TAR_FILE"
    
    log "🧹 Cleaning up tar file..."
    rm -f "$TAR_FILE"
    
    USED_METHOD="tar"
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))
DURATION_SEC=$((DURATION % 60))
SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)

log "============================================="
log "✅ LOCAL BUILD COMPLETE"
log "============================================="
log "📦 File: $OUTPUT_PATH"
log "📊 Size: $SIZE"
log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"
log "🎯 Method: $USED_METHOD"
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