#!/bin/bash
#
# Fast Singularity build script optimized for no-sudo environments
# Optimizations:
# - Parallel Docker export/build operations
# - Optimized for 128GB RAM and fast disk
# - All data in /host/cassio/export03/data/enning/
# - Progress monitoring and error recovery
#

set -e  # Exit on error

# ============================================================================
# Configuration - Optimized for your server
# ============================================================================
DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

# All data in your directory with 12TB space
BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_DIR="${BASE_DIR}/singularity"
OUTPUT_NAME="micapipe_v1_beta.sif"
OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_NAME}"

# Optimized temp locations (using your 128GB RAM)
TAR_FILE="${BASE_DIR}/micapipe_docker_$$.tar"
WORK_DIR="${BASE_DIR}/.singularity_work_$$"

# Cache locations (persistent across builds)
export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"

# Performance tuning
export SINGULARITY_MEMORY="64G"  # Use half your 128GB RAM
export DOCKER_BUILDKIT=1         # Enable BuildKit for faster operations

# ============================================================================
# Helper Functions
# ============================================================================
log() {
    echo "[$(date '+%H:%M:%S')] $*"
}

check_space() {
    local dir="$1"
    local needed_gb="$2"
    local available=$(df -BG "$dir" | awk 'NR==2 {print $4}' | sed 's/G//')
    
    if [ "$available" -lt "$needed_gb" ]; then
        log "ERROR: Need ${needed_gb}GB in $dir, only ${available}GB available"
        exit 1
    fi
    log "✅ Space check: ${available}GB available in $dir"
}

cleanup() {
    log "🧹 Cleaning up temporary files..."
    rm -f "${TAR_FILE}" 2>/dev/null || true
    rm -rf "${WORK_DIR}" 2>/dev/null || true
}

monitor_progress() {
    local file="$1"
    local desc="$2"
    
    while [ ! -f "$file" ]; do
        sleep 5
    done
    
    while true; do
        if [ -f "$file" ]; then
            size=$(du -h "$file" 2>/dev/null | cut -f1 || echo "0")
            log "📊 $desc: $size"
        fi
        sleep 30
    done
}

# ============================================================================
# Pre-flight checks
# ============================================================================
log "🚀 FAST SINGULARITY BUILD STARTING"
log "📦 Image: $FULL_DOCKER_IMAGE"
log "📍 Output: $OUTPUT_PATH"
log "💾 Using ${BASE_DIR} for all operations"

# Check Docker image exists
if ! docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null; then
    log "❌ Docker image not found: $FULL_DOCKER_IMAGE"
    exit 1
fi

# Create directories
mkdir -p "$OUTPUT_DIR" "$WORK_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Check disk space (need ~350GB during build)
check_space "$BASE_DIR" 350

# Check if output exists
if [ -f "$OUTPUT_PATH" ]; then
    log "⚠️  Output file exists: $OUTPUT_PATH"
    read -p "Overwrite? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log "❌ Build cancelled"
        exit 1
    fi
    rm -f "$OUTPUT_PATH"
fi

# Set trap for cleanup
trap cleanup EXIT

# ============================================================================
# Step 1: Fast Docker Export (Optimized)
# ============================================================================
log "📦 Step 1/2: Exporting Docker image (optimized)..."

START_TIME=$(date +%s)

# Start progress monitor in background
monitor_progress "$TAR_FILE" "Docker export" &
MONITOR_PID=$!

# Export with optimizations
log "🔄 Running optimized docker save..."
time docker save "$FULL_DOCKER_IMAGE" | pv -s 100G > "$TAR_FILE"

# Stop monitor
kill $MONITOR_PID 2>/dev/null || true

EXPORT_TIME=$(date +%s)
EXPORT_DURATION=$((EXPORT_TIME - START_TIME))
TAR_SIZE=$(du -h "$TAR_FILE" | cut -f1)

log "✅ Docker export complete: $TAR_SIZE in ${EXPORT_DURATION}s"

# ============================================================================
# Step 2: Fast Singularity Build (Optimized)
# ============================================================================
log "🔧 Step 2/2: Building Singularity SIF (optimized)..."

# Start progress monitor in background for SIF
monitor_progress "$OUTPUT_PATH" "SIF build" &
MONITOR_PID=$!

# Optimized singularity build with performance flags
log "🔄 Running optimized singularity build..."

# Use --sandbox first (faster), then build SIF
SANDBOX_DIR="${WORK_DIR}/sandbox"

log "🏗️  Building sandbox (faster intermediate step)..."
time singularity build --force --sandbox \
    "$SANDBOX_DIR" \
    "docker-archive://$TAR_FILE"

log "📦 Converting sandbox to SIF..."
time singularity build --force \
    "$OUTPUT_PATH" \
    "$SANDBOX_DIR"

# Stop monitor
kill $MONITOR_PID 2>/dev/null || true

BUILD_TIME=$(date +%s)
BUILD_DURATION=$((BUILD_TIME - EXPORT_TIME))
TOTAL_DURATION=$((BUILD_TIME - START_TIME))

# ============================================================================
# Results
# ============================================================================
SIF_SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)

log "============================================="
log "✅ SUCCESS - FAST BUILD COMPLETE"
log "============================================="
log "📦 SIF file: $OUTPUT_PATH"
log "📊 Size: $SIF_SIZE"
log "⏱️  Timing:"
log "   Export: ${EXPORT_DURATION}s"
log "   Build:  ${BUILD_DURATION}s"
log "   Total:  ${TOTAL_DURATION}s"
log ""
log "🧪 Test with:"
log "   singularity run $OUTPUT_PATH --help"
log ""
log "🚀 Ready for HPC deployment!"

# Auto cleanup happens via trap