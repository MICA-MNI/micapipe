#!/bin/bash
#
# Local-only fast Singularity build - no registry access needed
# Optimized for local Docker images with 128GB RAM server
# Uses TAR method to avoid Docker temp space issues
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
# TAR Method: Most reliable - avoids Docker temp directory issues
# The streaming method uses Docker's internal temp (/export01/docker/tmp)
# which is full. TAR method writes directly to our specified path.
# ============================================================================
TAR_FILE="${BASE_DIR}/.tmp/micapipe_docker_$$.tar"

log "📤 Exporting Docker image to tar file..."
log "📍 Tar location: $TAR_FILE"

# Monitor tar creation in background
(
    sleep 10  # Give docker save time to start
    while [ -f "$TAR_FILE" ]; do
        SIZE=$(du -h "$TAR_FILE" 2>/dev/null | cut -f1 || echo "0")
        log "📈 Tar export progress: $SIZE"
        sleep 30
    done
) &
TAR_MONITOR_PID=$!

docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"
# Create tar file - prefer skopeo when DockerRootDir is on /export01 (avoids daemon tmp)
DOCKER_ROOT_DIR=$(docker info --format '{{.DockerRootDir}}' 2>/dev/null || true)
if [[ -n "$DOCKER_ROOT_DIR" && "$DOCKER_ROOT_DIR" == /export01* ]]; then
    log "⚠️  Docker daemon root is on $DOCKER_ROOT_DIR which may be low on space"
    if command -v skopeo >/dev/null 2>&1; then
        log "🔁 Using skopeo to copy image directly to docker-archive (no daemon temporary files)"
        if skopeo copy --override-os linux "docker://$FULL_DOCKER_IMAGE" "docker-archive:$TAR_FILE"; then
            DOCKER_EXIT=0
        else
            DOCKER_EXIT=2
        fi
    else
        log "❌ skopeo not found. Can't safely export image without using Docker daemon temp at $DOCKER_ROOT_DIR"
        log "ℹ️  Options: 1) install skopeo, 2) reconfigure Docker to use /export03, or 3) free space on /export01"
        rm -f "$TAR_FILE"
        exit 1
    fi
else
    # Docker root dir is not on /export01 — use docker save to write tar directly
    docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"
    DOCKER_EXIT=$?
fi

# Stop tar monitor
kill $TAR_MONITOR_PID 2>/dev/null || true
wait $TAR_MONITOR_PID 2>/dev/null || true

if [ $DOCKER_EXIT -ne 0 ]; then
    log "❌ Docker save failed with exit code $DOCKER_EXIT"
    rm -f "$TAR_FILE"
    exit 1
fi

TAR_SIZE=$(du -h "$TAR_FILE" | cut -f1)
log "✅ Docker export complete: $TAR_SIZE"

log "🔧 Building SIF from tar..."

# Monitor SIF creation in background
(
    sleep 10  # Give singularity time to start
    while [ -f "$OUTPUT_PATH" ] || [ ! -f "$OUTPUT_PATH" ]; do
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
SIF_MONITOR_PID=$!

# Build SIF from tar
singularity build --force "$OUTPUT_PATH" "docker-archive://$TAR_FILE"
SINGULARITY_EXIT=$?

# Stop SIF monitor
kill $SIF_MONITOR_PID 2>/dev/null || true
wait $SIF_MONITOR_PID 2>/dev/null || true

# Clean up tar file
log "🧹 Cleaning up tar file..."
rm -f "$TAR_FILE"

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