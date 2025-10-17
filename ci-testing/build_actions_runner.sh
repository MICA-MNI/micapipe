#!/bin/bash
#
# Build GitHub Actions runner with embedded micapipe SIF
# This creates a self-contained CI runner with micapipe baked in
#

set -e

# Paths
ACTIONS_RUNNER_DIR="/data_/mica1/03_projects/actions-runner"
SIF_SOURCE="/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif"
SIF_TARGET="$ACTIONS_RUNNER_DIR/micapipe_v1_beta.sif"
DOCKERFILE="$ACTIONS_RUNNER_DIR/Dockerfile"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🏗️  BUILDING GITHUB ACTIONS RUNNER WITH EMBEDDED MICAPIPE"
log "============================================================"

# ============================================================================
# Step 1: Copy SIF to actions-runner directory
# ============================================================================
log "📋 Step 1: Copying SIF file to actions-runner directory..."

if [ ! -f "$SIF_SOURCE" ]; then
    log "❌ SIF file not found: $SIF_SOURCE"
    log "   Run ./build_singularity_v1.sh first to create the SIF file"
    exit 1
fi

SIF_SIZE=$(du -h "$SIF_SOURCE" | cut -f1)
log "📦 Found SIF file: $SIF_SIZE"

log "📥 Copying SIF to actions-runner directory..."
cp "$SIF_SOURCE" "$SIF_TARGET"
log "✅ SIF copied to: $SIF_TARGET"

# ============================================================================
# Step 2: Update Dockerfile to include the SIF
# ============================================================================
log "📋 Step 2: Updating Dockerfile to embed SIF..."

# Check if the COPY line is already uncommented
if grep -q "^COPY micapipe_v1_beta.sif" "$DOCKERFILE"; then
    log "✅ Dockerfile already configured to copy SIF"
else
    # Uncomment the COPY line
    if grep -q "^# COPY micapipe_v0.2.3.sif" "$DOCKERFILE"; then
        log "🔧 Updating COPY line in Dockerfile..."
        sed -i 's|^# COPY micapipe_v0.2.3.sif /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v0.2.3.sif|COPY micapipe_v1_beta.sif /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v1_beta.sif|' "$DOCKERFILE"
        log "✅ Dockerfile updated"
    else
        log "⚠️  Adding COPY line to Dockerfile..."
        # Add the COPY line after the mkdir command
        sed -i '/RUN mkdir -p \/data_\/mica1\/01_programs\/micapipe-v0.2.0/a COPY micapipe_v1_beta.sif /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v1_beta.sif' "$DOCKERFILE"
        log "✅ COPY line added to Dockerfile"
    fi
fi

# ============================================================================
# Step 3: Build the Docker image
# ============================================================================
log "📋 Step 3: Building Docker image with embedded SIF..."

cd "$ACTIONS_RUNNER_DIR"

log "🔄 Building micapipe-runner Docker image..."
docker build -t micapipe-runner .

log "✅ Docker image built successfully"

# ============================================================================
# Step 4: Verification and next steps
# ============================================================================
log "🧪 Verifying the build..."
IMAGE_SIZE=$(docker images micapipe-runner --format "table {{.Size}}" | tail -n1)
log "📊 Docker image size: $IMAGE_SIZE"

log "============================================="
log "✅ GITHUB ACTIONS RUNNER BUILD COMPLETE"
log "============================================="
log "📦 Docker image: micapipe-runner"
log "📊 Image size: $IMAGE_SIZE"
log "🗃️  Embedded SIF: /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v1_beta.sif"
log ""
log "🚀 Next steps:"
log "   1. Test the runner:"
log "      docker run --rm micapipe-runner ls -la /data_/mica1/01_programs/micapipe-v0.2.0/"
log ""
log "   2. Deploy the runner:"
log "      ./run_docker.sh"
log ""
log "   3. The runner will be available for CI jobs with embedded micapipe!"
log ""
log "📋 The runner includes:"
log "   - GitHub Actions runner"
log "   - Singularity 3.10.2-bionic"
log "   - Docker support"
log "   - Embedded micapipe SIF (~${SIF_SIZE})"
log ""
log "🎯 Ready for CI and standalone use!"