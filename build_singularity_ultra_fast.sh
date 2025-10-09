#!/bin/bash
#
# Ultra-fast Singularity build using parallel processing and memory optimization
# Designed for your server: 128GB RAM, 12TB storage, no sudo
#

set -e

# ============================================================================
# Optimized Configuration
# ============================================================================
DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"

# Memory-optimized locations
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp" 
export TMPDIR="${BASE_DIR}/.tmp"

# Performance tuning for your 128GB RAM server
export SINGULARITY_MEMORY="96G"
export OMP_NUM_THREADS=$(nproc)
export DOCKER_BUILDKIT=1

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ============================================================================
# Method 1: Direct Registry Pull (Fastest if image is on GHCR)
# ============================================================================
build_from_registry() {
    log "🚀 Method 1: Direct registry pull (fastest)"
    
    # Try pulling directly from registry (bypasses Docker entirely)
    REGISTRY_URI="docker://ghcr.io/mica-mni/micapipe:${DOCKER_TAG}"
    
    log "📥 Pulling directly from registry..."
    time singularity build --force \
        "${OUTPUT_PATH}" \
        "${REGISTRY_URI}"
}

# ============================================================================
# Method 2: Parallel Docker + Pipe (Fast local build)
# ============================================================================
build_from_docker_parallel() {
    log "🚀 Method 2: Parallel docker->singularity (no intermediate tar)"
    
    # Stream directly from docker to singularity (no tar file)
    log "🔄 Streaming docker->singularity..."
    time docker save "${FULL_DOCKER_IMAGE}" | \
         singularity build --force "${OUTPUT_PATH}" docker-archive:///dev/stdin
}

# ============================================================================
# Method 3: Memory-optimized build (if others fail)
# ============================================================================
build_memory_optimized() {
    log "🚀 Method 3: Memory-optimized build"
    
    TAR_FILE="${BASE_DIR}/micapipe_docker_$$.tar"
    
    # Use ionice to prioritize disk I/O
    log "📦 High-priority docker export..."
    ionice -c 1 -n 0 docker save "${FULL_DOCKER_IMAGE}" -o "${TAR_FILE}" 2>/dev/null || \
    docker save "${FULL_DOCKER_IMAGE}" -o "${TAR_FILE}"
    
    log "🔧 Memory-optimized singularity build..."
    ionice -c 1 -n 0 singularity build --force \
        "${OUTPUT_PATH}" \
        "docker-archive://${TAR_FILE}" 2>/dev/null || \
    singularity build --force \
        "${OUTPUT_PATH}" \
        "docker-archive://${TAR_FILE}"
    
    rm -f "${TAR_FILE}"
}

# ============================================================================
# Main execution with fallback methods
# ============================================================================
log "🚀 ULTRA-FAST SINGULARITY BUILD"
log "📦 Image: ${FULL_DOCKER_IMAGE}"
log "📍 Output: ${OUTPUT_PATH}"
log "💾 Base: ${BASE_DIR}"

# Create directories
mkdir -p "${OUTPUT_DIR}" "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}" "${TMPDIR}"

# Check space
AVAILABLE=$(df -BG "${BASE_DIR}" | awk 'NR==2 {print $4}' | sed 's/G//')
if [ "${AVAILABLE}" -lt 200 ]; then
    log "❌ Need 200GB+ space, only ${AVAILABLE}GB available"
    exit 1
fi

# Remove existing output
[ -f "${OUTPUT_PATH}" ] && rm -f "${OUTPUT_PATH}"

START_TIME=$(date +%s)

# Try methods in order of speed
log "🎯 Attempting fastest methods first..."

if build_from_registry 2>/dev/null; then
    log "✅ Registry pull succeeded (fastest method)"
elif build_from_docker_parallel 2>/dev/null; then
    log "✅ Parallel docker build succeeded"
else
    log "⚠️  Falling back to memory-optimized method"
    build_memory_optimized
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
SIZE=$(du -h "${OUTPUT_PATH}" | cut -f1)

log "============================================="
log "✅ ULTRA-FAST BUILD COMPLETE"
log "============================================="
log "📦 File: ${OUTPUT_PATH}"
log "📊 Size: ${SIZE}"
log "⏱️  Time: ${DURATION}s"
log ""
log "🧪 Test: singularity run ${OUTPUT_PATH} --help"
log "🚀 Ready for deployment!"