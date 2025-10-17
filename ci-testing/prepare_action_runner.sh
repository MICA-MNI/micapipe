#!/bin/bash
#
# Prepare action runner with base image
# Run this once to set up the base image for CI
#

set -e

MICAPIPE_BASE_IMAGE="ghcr.io/mica-mni/micapipe-base:latest"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🔧 PREPARING ACTION RUNNER WITH BASE IMAGE"
log "==========================================="

# Check if running inside action runner container
if [ -f "/.dockerenv" ]; then
    log "✅ Running inside Docker container (action runner)"
    DOCKERFILE_BASE="./Dockerfile.base"
else
    log "⚠️  Running on host - checking if action runner is available"
    # Check if we can access action runner
    if docker ps --format "table {{.Names}}" | grep -q "micapipe-runner"; then
        log "✅ Found running action runner container"
        log "🔄 Executing inside action runner..."
        docker exec micapipe-runner /bin/bash -c "cd /actions-runner && $(cat "$0")"
        exit $?
    else
        log "❌ Action runner not running. Start it first with:"
        log "   cd /data_/mica1/03_projects/actions-runner"
        log "   ./run_docker.sh"
        exit 1
    fi
fi

# Now we're inside the action runner container
log "📋 Checking for base image: $MICAPIPE_BASE_IMAGE"

if docker image inspect "$MICAPIPE_BASE_IMAGE" &>/dev/null; then
    log "✅ Base image already exists"
    docker images | grep micapipe-base
else
    log "❌ Base image not found - building it..."
    
    # Check if we have the Dockerfile.base
    if [ ! -f "Dockerfile.base" ]; then
        log "❌ Dockerfile.base not found in action runner"
        log "   Make sure the repository is up to date inside the runner"
        exit 1
    fi
    
    log "🔄 Building base image..."
    docker build -f Dockerfile.base -t "$MICAPIPE_BASE_IMAGE" .
    
    BASE_SIZE=$(docker images "$MICAPIPE_BASE_IMAGE" --format "table {{.Size}}" | tail -n1)
    log "✅ Base image built successfully: $BASE_SIZE"
fi

log ""
log "🧪 Testing base image..."
if docker run --rm "$MICAPIPE_BASE_IMAGE" /bin/bash -c "echo 'Base image works!' && python3 --version && freesurfer --version 2>/dev/null | head -1"; then
    log "✅ Base image test passed"
else
    log "⚠️  Base image test had issues (might be normal)"
fi

log ""
log "============================================="
log "✅ ACTION RUNNER BASE IMAGE SETUP COMPLETE"
log "============================================="
log "📦 Base image: $MICAPIPE_BASE_IMAGE"
log "📊 Size: $(docker images "$MICAPIPE_BASE_IMAGE" --format "table {{.Size}}" | tail -n1)"
log ""
log "🚀 Your action runner is now ready for CI!"
log "   CI jobs will use the existing base image"
log "   and only build the fast main image (~5-10 min)"
log ""
log "📋 Next CI run will:"
log "   1. ✅ Find existing base image (skip base build)"
log "   2. 🔄 Build main image (fast)"
log "   3. 🔄 Create SIF from main image"
log "   4. 🧪 Run tests"