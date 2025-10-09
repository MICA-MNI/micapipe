#!/bin/bash
#
# Smart Singularity build that adapts to your Docker setup
# Automatically detects local vs remote images and chooses best method
#

set -e

DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"

# Performance settings
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================================
# Detect available Docker images and choose best strategy
# ============================================================================
detect_image_source() {
    log "🔍 Detecting best image source..."
    
    # Check local Docker images
    if docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null; then
        LOCAL_SIZE=$(docker image inspect "$FULL_DOCKER_IMAGE" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')
        log "✅ Found local image: $FULL_DOCKER_IMAGE ($LOCAL_SIZE)"
        HAS_LOCAL=true
    else
        log "❌ No local image: $FULL_DOCKER_IMAGE"
        HAS_LOCAL=false
    fi
    
    # Check if registry is accessible (try to pull manifest)
    log "🌐 Checking registry availability..."
    if timeout 10 singularity inspect "docker://$FULL_DOCKER_IMAGE" &>/dev/null; then
        log "✅ Registry accessible: docker://$FULL_DOCKER_IMAGE"
        HAS_REGISTRY=true
    else
        log "❌ Registry not accessible or image not found"
        HAS_REGISTRY=false
    fi
    
    # Choose best strategy
    if [ "$HAS_REGISTRY" = true ]; then
        STRATEGY="registry"
        log "🎯 Strategy: Direct registry pull (fastest)"
    elif [ "$HAS_LOCAL" = true ]; then
        STRATEGY="local"
        log "🎯 Strategy: Local Docker conversion"
    else
        log "❌ ERROR: No image source available!"
        log "   Options:"
        log "   1. Push your local image: docker push $FULL_DOCKER_IMAGE"
        log "   2. Build image locally first"
        exit 1
    fi
}

# ============================================================================
# Strategy 1: Direct registry pull
# ============================================================================
build_from_registry() {
    log "🚀 Building from registry (fastest method)"
    
    time singularity build --force \
        "$OUTPUT_PATH" \
        "docker://$FULL_DOCKER_IMAGE"
}

# ============================================================================
# Strategy 2: Local Docker optimized
# ============================================================================
build_from_local() {
    log "🚀 Building from local Docker (optimized)"
    
    # Method A: Try streaming (fastest)
    log "⚡ Trying streaming method..."
    if timeout 300 bash -c "docker save '$FULL_DOCKER_IMAGE' | singularity build --force '$OUTPUT_PATH' docker-archive:///dev/stdin" 2>/dev/null; then
        log "✅ Streaming method succeeded"
        return 0
    fi
    
    # Method B: Tar method (more reliable)
    log "📦 Falling back to tar method..."
    TAR_FILE="${BASE_DIR}/micapipe_docker_$$.tar"
    
    log "📤 Exporting Docker image..."
    docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"
    
    log "🔧 Building SIF from tar..."
    singularity build --force \
        "$OUTPUT_PATH" \
        "docker-archive://$TAR_FILE"
    
    log "🧹 Cleaning up tar file..."
    rm -f "$TAR_FILE"
}

# ============================================================================
# Main execution
# ============================================================================
log "🚀 SMART SINGULARITY BUILD"
log "📦 Image: $FULL_DOCKER_IMAGE"
log "📍 Output: $OUTPUT_PATH"

# Create directories
mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Check space
AVAILABLE=$(df -BG "$BASE_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "$AVAILABLE" -lt 200 ]; then
    log "❌ Need 200GB+ space, only ${AVAILABLE}GB available"
    exit 1
fi

# Remove existing output
[ -f "$OUTPUT_PATH" ] && rm -f "$OUTPUT_PATH"

# Detect best source and strategy
detect_image_source

START_TIME=$(date +%s)

# Execute chosen strategy
case "$STRATEGY" in
    "registry")
        build_from_registry
        ;;
    "local")
        build_from_local
        ;;
esac

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)

log "============================================="
log "✅ SMART BUILD COMPLETE"
log "============================================="
log "📦 File: $OUTPUT_PATH"
log "📊 Size: $SIZE"
log "⏱️  Time: ${DURATION}s"
log "🎯 Strategy: $STRATEGY"
log ""
log "🧪 Test: singularity run $OUTPUT_PATH --help"

# ============================================================================
# Push recommendation
# ============================================================================
if [ "$STRATEGY" = "local" ] && [ "$HAS_REGISTRY" = false ]; then
    log ""
    log "💡 RECOMMENDATION:"
    log "   Push your image to speed up future builds:"
    log "   docker push $FULL_DOCKER_IMAGE"
    log ""
    log "   This enables the fastest registry-direct method"
fi