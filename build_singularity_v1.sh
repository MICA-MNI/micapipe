#!/bin/bash
#
# Simple Singularity build - Original v1 method
# No fancy optimizations, just the basic reliable approach
#

set -e

DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

# Use standard paths (like v1)
OUTPUT_DIR="/host/cassio/export03/data/enning/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"
TAR_FILE="/host/cassio/export03/data/enning/micapipe_docker.tar"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================================
# Simple build process (v1 style)
# ============================================================================
log "🔧 SIMPLE SINGULARITY BUILD (v1 method)"
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

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Clean up any previous attempts
log "🧹 Cleaning up previous files..."
rm -f "$OUTPUT_PATH" "$TAR_FILE"

START_TIME=$(date +%s)

# ============================================================================
# Step 1: Export Docker to tar
# ============================================================================
log "📤 Step 1: Exporting Docker image to tar..."
docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"

TAR_SIZE=$(du -h "$TAR_FILE" | cut -f1)
log "✅ Export complete: $TAR_SIZE"

# ============================================================================
# Step 2: Build SIF from tar
# ============================================================================
log "🔧 Step 2: Building SIF from tar..."
singularity build "$OUTPUT_PATH" "docker-archive://$TAR_FILE"

SIF_SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)
log "✅ SIF build complete: $SIF_SIZE"

# ============================================================================
# Step 3: Cleanup
# ============================================================================
log "🧹 Cleaning up tar file..."
rm -f "$TAR_FILE"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
DURATION_MIN=$((DURATION / 60))
DURATION_SEC=$((DURATION % 60))

log "============================================="
log "✅ SIMPLE BUILD COMPLETE"
log "============================================="
log "📦 File: $OUTPUT_PATH"
log "📊 Size: $SIF_SIZE"
log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"
log ""
log "🧪 Test with:"
log "   singularity run $OUTPUT_PATH --help"
log ""
log "📋 Next steps for GitHub Actions runner:"
log "   1. Copy SIF to actions-runner directory:"
log "      cp $OUTPUT_PATH /data_/mica1/03_projects/actions-runner/micapipe_v1_beta.sif"
log "   2. Update Dockerfile to copy the SIF"
log "   3. Build new action runner with embedded SIF"
log ""
log "🚀 Ready!"

# Final verification
if [ -f "$OUTPUT_PATH" ] && [ -s "$OUTPUT_PATH" ]; then
    log "✅ SIF file created successfully"
else
    log "❌ ERROR: SIF file not created or empty"
    exit 1
fi