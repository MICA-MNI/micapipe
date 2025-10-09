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

# Function to check and kill stuck processes
check_stuck_processes() {
    log "🔍 Checking for stuck processes..."
    
    # Check for singularity processes
    SING_PROCS=$(ps aux | grep singularity | grep -v grep | wc -l)
    if [ "$SING_PROCS" -gt 0 ]; then
        log "🔍 Found $SING_PROCS singularity processes:"
        ps aux | grep singularity | grep -v grep
    fi
    
    # Check for docker save processes  
    DOCKER_PROCS=$(ps aux | grep "docker save" | grep -v grep | wc -l)
    if [ "$DOCKER_PROCS" -gt 0 ]; then
        log "🔍 Found $DOCKER_PROCS docker save processes:"
        ps aux | grep "docker save" | grep -v grep
    fi
    
    # Check temp directory usage
    if [ -d "$SINGULARITY_TMPDIR" ]; then
        TEMP_USAGE=$(du -sh "$SINGULARITY_TMPDIR" 2>/dev/null | cut -f1 || echo "0")
        log "💾 Temp directory: $TEMP_USAGE"
    fi
    
    # Check if output file exists and its size
    if [ -f "$OUTPUT_PATH" ]; then
        OUTPUT_SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)
        log "📁 Output file exists: $OUTPUT_SIZE"
    else
        log "❌ No output file yet"
    fi
}

# Function to kill all related processes
kill_stuck_build() {
    log "🛑 Killing stuck build processes..."
    
    # Kill singularity processes
    pkill -f singularity 2>/dev/null || true
    
    # Kill docker save processes
    pkill -f "docker save" 2>/dev/null || true
    
    # Clean up temp files
    if [ -d "$SINGULARITY_TMPDIR" ]; then
        rm -rf "$SINGULARITY_TMPDIR"/* 2>/dev/null || true
    fi
    
    # Remove partial output
    if [ -f "$OUTPUT_PATH" ]; then
        rm -f "$OUTPUT_PATH"
    fi
    
    log "✅ Cleanup complete"
}

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
log "📊 Progress monitoring: watching output file size..."

# Start progress monitor in background
monitor_progress() {
    while true; do
        if [ -f "$OUTPUT_PATH" ]; then
            SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")
            log "📈 Current size: $SIZE"
        else
            log "⏳ Waiting for output file to appear..."
        fi
        sleep 30
    done
}

# Start monitoring
monitor_progress &
MONITOR_PID=$!

# Run the actual build with verbose output
log "🔄 Starting docker save | singularity build..."
if timeout 7200 bash -c "docker save '$FULL_DOCKER_IMAGE' | singularity build --force '$OUTPUT_PATH' docker-archive:///dev/stdin" 2>&1 | while read line; do
    log "SINGULARITY: $line"
done; then
    # Stop monitor
    kill $MONITOR_PID 2>/dev/null || true
    wait $MONITOR_PID 2>/dev/null || true
    log "✅ Streaming method succeeded!"
    USED_METHOD="streaming"
else
    # Stop monitor  
    kill $MONITOR_PID 2>/dev/null || true
    wait $MONITOR_PID 2>/dev/null || true
    log "⚠️  Streaming failed, trying tar method..."
    
    # ============================================================================
    # Method 2: Tar method (more reliable)
    # ============================================================================
    TAR_FILE="${BASE_DIR}/micapipe_docker_$$.tar"
    
    log "📤 Exporting Docker to tar..."
    
    # Monitor tar creation
    (while [ ! -f "$TAR_FILE" ] || [ $(stat -f%z "$TAR_FILE" 2>/dev/null || echo 0) -eq 0 ]; do
        sleep 5
        log "⏳ Waiting for tar export to start..."
    done
    
    while [ -f "$TAR_FILE" ] && kill -0 $! 2>/dev/null; do
        SIZE=$(du -h "$TAR_FILE" 2>/dev/null | cut -f1 || echo "0")
        log "📈 Tar progress: $SIZE"
        sleep 30
    done) &
    
    TAR_MONITOR_PID=$!
    
    # Create tar file
    docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE" &
    DOCKER_PID=$!
    
    # Wait for docker save to complete
    wait $DOCKER_PID
    
    # Stop tar monitor
    kill $TAR_MONITOR_PID 2>/dev/null || true
    wait $TAR_MONITOR_PID 2>/dev/null || true
    
    TAR_SIZE=$(du -h "$TAR_FILE" | cut -f1)
    log "✅ Export complete: $TAR_SIZE"
    
    log "🔧 Building SIF from tar..."
    
    # Monitor SIF creation
    (while [ ! -f "$OUTPUT_PATH" ] || [ $(stat -f%z "$OUTPUT_PATH" 2>/dev/null || echo 0) -eq 0 ]; do
        sleep 5
        log "⏳ Waiting for SIF build to start..."
    done
    
    while [ -f "$OUTPUT_PATH" ] && kill -0 $! 2>/dev/null; do
        SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")
        log "📈 SIF progress: $SIZE"
        sleep 30
    done) &
    
    SIF_MONITOR_PID=$!
    
    # Build SIF
    singularity build --force \
        "$OUTPUT_PATH" \
        "docker-archive://$TAR_FILE" &
    SINGULARITY_PID=$!
    
    # Wait for singularity build to complete
    wait $SINGULARITY_PID
    
    # Stop SIF monitor
    kill $SIF_MONITOR_PID 2>/dev/null || true
    wait $SIF_MONITOR_PID 2>/dev/null || true
    
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