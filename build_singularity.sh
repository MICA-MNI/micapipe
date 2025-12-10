#!/bin/bash#!/bin/bash#!/bin/bash

#

# Singularity build from local Docker image##

# Uses docker-daemon:// to read directly from local Docker (no network download)

## Singularity build directly from GitHub Container Registry (GHCR)# Singularity build from local Docker image



set -e# Bypasses local Docker daemon entirely - no /export01 temp space needed# Uses docklog "🔧 Building SIF directly from Docker daemon..."r-daemon:// protocol - Singularity reads directly from Docker daemon



DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"## No intermediate tar files, no temp space issues

DOCKER_TAG="${1:-latest}"

FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"#



BASE_DIR="/export03/data/enning"set -e

OUTPUT_DIR="${BASE_DIR}/singularity"

OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"set -e



# Use export03 for Singularity temp/cacheDOCKER_IMAGE="ghcr.io/mica-mni/micapipe"

export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"

export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"DOCKER_TAG="${1:-latest}"DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"

export TMPDIR="${BASE_DIR}/.tmp"

export TEMP="${BASE_DIR}/.tmp"FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"DOCKER_TAG="${1:-latest}"

export TMP="${BASE_DIR}/.tmp"

FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

BASE_DIR="/export03/data/enning"

# ============================================================================

# Pre-flight checksOUTPUT_DIR="${BASE_DIR}/singularity"BASE_DIR="/export03/data/enning"

# ============================================================================

log "🚀 SINGULARITY BUILD (from local Docker image)"OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"OUTPUT_DIR="${BASE_DIR}/singularity"

log "📦 Source: docker-daemon://${FULL_DOCKER_IMAGE}"

log "📍 Output: $OUTPUT_PATH"OUTPUT_PATH="${OUTPUT_DIR}/micapipe_v1_beta.sif"



# Check local Docker image exists# Use export03 for ALL temp/cache files

if ! docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null; then

    log "❌ Local Docker image not found: $FULL_DOCKER_IMAGE"export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"# Use export03 for ALL temp files to avoid /export01 space issues

    docker images | grep micapipe || echo "   No micapipe images found"

    exit 1export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"export SINGULARITY_CACHEDIR="${BASE_DIR}/.singularity_cache"

fi

export TMPDIR="${BASE_DIR}/.tmp"export SINGULARITY_TMPDIR="${BASE_DIR}/.singularity_tmp"

LOCAL_SIZE=$(docker image inspect "$FULL_DOCKER_IMAGE" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')

log "✅ Found local image: $LOCAL_SIZE"export TEMP="${BASE_DIR}/.tmp"export TMPDIR="${BASE_DIR}/.tmp"



# Create directoriesexport TMP="${BASE_DIR}/.tmp"export TEMP="${BASE_DIR}/.tmp"

mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "${BASE_DIR}/.tmp"

export TMP="${BASE_DIR}/.tmp"

# Check space on export03

AVAILABLE=$(df -BG "$BASE_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "✅ Space on /export03: ${AVAILABLE}GB available"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# Check space on /export01 (Docker temp)

DOCKER_SPACE=$(df -BG /export01 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//' || echo "0")# ============================================================================

log "✅ Space on /export01 (Docker temp): ${DOCKER_SPACE}GB available"

# Pre-flight checks# ============================================================================

if [ "$DOCKER_SPACE" -lt 100 ]; then

    log "⚠️  Warning: /export01 has less than 100GB. Docker export may fail."# ============================================================================# Pre-flight checks

    log "   Consider running: docker system prune -f"

filog "🚀 SINGULARITY BUILD (direct from GHCR - no Docker daemon)"# ============================================================================



# Remove existing outputlog "📦 Source: docker://${FULL_DOCKER_IMAGE}"log "🚀 SINGULARITY BUILD (docker-daemon:// method)"

if [ -f "$OUTPUT_PATH" ]; then

    log "⚠️  Removing existing SIF file"log "📍 Output: $OUTPUT_PATH"log "📦 Image: $FULL_DOCKER_IMAGE"

    rm -f "$OUTPUT_PATH"

filog "📍 Output: $OUTPUT_PATH"



START_TIME=$(date +%s)# Create directories



# ============================================================================mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "${BASE_DIR}/.tmp"# Check local Docker image exists

# Build using docker-daemon:// - reads from local Docker, no download needed

# ============================================================================if ! docker image inspect "$FULL_DOCKER_IMAGE" &>/dev/null; then



log "🔧 Building SIF from local Docker daemon..."# Check space on export03    log "❌ Local Docker image not found: $FULL_DOCKER_IMAGE"



# Monitor progress in backgroundAVAILABLE=$(df -BG "$BASE_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')    log "   Available images:"

(

    sleep 30if [ "$AVAILABLE" -lt 200 ]; then    docker images | grep micapipe || echo "   No micapipe images found"

    while true; do

        if [ -f "$OUTPUT_PATH" ]; then    log "❌ Need 200GB+ space, only ${AVAILABLE}GB available"    exit 1

            SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")

            log "📈 SIF build progress: $SIZE"    exit 1fi

        else

            log "⏳ Singularity extracting layers from Docker..."fi

        fi

        sleep 30log "✅ Space check: ${AVAILABLE}GB available on /export03"LOCAL_SIZE=$(docker image inspect "$FULL_DOCKER_IMAGE" --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1024/1024/1024}')

        if ! kill -0 $$ 2>/dev/null; then

            breaklog "✅ Found local image: $LOCAL_SIZE"

        fi

    done# Remove existing output

) &

MONITOR_PID=$!if [ -f "$OUTPUT_PATH" ]; then# Create directories



# Build from local Docker daemon    log "⚠️  Removing existing SIF file"mkdir -p "$OUTPUT_DIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "${BASE_DIR}/.tmp"

singularity build --force "$OUTPUT_PATH" "docker-daemon://${FULL_DOCKER_IMAGE}"

SINGULARITY_EXIT=$?    rm -f "$OUTPUT_PATH"



# Stop monitorfi# Check space

kill $MONITOR_PID 2>/dev/null || true

wait $MONITOR_PID 2>/dev/null || trueAVAILABLE=$(df -BG "$BASE_DIR" | awk 'NR==2 {print $4}' | sed 's/G//')



if [ $SINGULARITY_EXIT -ne 0 ]; thenSTART_TIME=$(date +%s)if [ "$AVAILABLE" -lt 200 ]; then

    log "❌ Singularity build failed with exit code $SINGULARITY_EXIT"

    exit 1    log "❌ Need 200GB+ space, only ${AVAILABLE}GB available"

fi

# ============================================================================    exit 1

END_TIME=$(date +%s)

DURATION=$((END_TIME - START_TIME))# Build directly from GHCR using docker:// protocolfi

DURATION_MIN=$((DURATION / 60))

DURATION_SEC=$((DURATION % 60))# This downloads from the registry directly - Docker daemon is NOT usedlog "✅ Space check: ${AVAILABLE}GB available"

SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)

# ============================================================================

log "============================================="

log "✅ BUILD COMPLETE"# Remove existing output

log "============================================="

log "📦 File: $OUTPUT_PATH"log "🔧 Building SIF directly from GHCR..."if [ -f "$OUTPUT_PATH" ]; then

log "📊 Size: $SIZE"

log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"log "   Singularity will download layers from ghcr.io (not local Docker)"    log "⚠️  Removing existing SIF file"

log ""

log "🧪 Test with:"    rm -f "$OUTPUT_PATH"

log "   singularity run $OUTPUT_PATH --help"

# Monitor SIF creation in backgroundfi

# Verify

if [ -f "$OUTPUT_PATH" ] && [ -s "$OUTPUT_PATH" ]; then(

    log "✅ SIF file created successfully"

else    sleep 30START_TIME=$(date +%s)

    log "❌ ERROR: SIF file not created or empty"

    exit 1    while true; do

fi

        if [ -f "$OUTPUT_PATH" ]; then# ============================================================================

            SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")# Build using docker-daemon:// protocol

            log "📈 SIF build progress: $SIZE"# Singularity reads directly from Docker daemon - no intermediate files needed

        else# This avoids all temp space issues on /export01

            # Check cache size as proxy for download progress# ============================================================================

            if [ -d "$SINGULARITY_CACHEDIR" ]; then

                CACHE_SIZE=$(du -sh "$SINGULARITY_CACHEDIR" 2>/dev/null | cut -f1 || echo "0")log "� Building SIF directly from Docker daemon..."

                log "⏳ Downloading from GHCR... (cache: $CACHE_SIZE)"log "   Using: docker-daemon://${FULL_DOCKER_IMAGE}"

            else

                log "⏳ Starting download from GHCR..."# Monitor SIF creation in background

            fi(

        fi    sleep 10

        sleep 30    while true; do

        # Check if parent process is still running        if [ -f "$OUTPUT_PATH" ]; then

        if ! kill -0 $$ 2>/dev/null; then            SIZE=$(du -h "$OUTPUT_PATH" 2>/dev/null | cut -f1 || echo "0")

            break            log "📈 SIF build progress: $SIZE"

        fi        else

    done            log "⏳ SIF build in progress..."

) &        fi

MONITOR_PID=$!        sleep 30

        # Check if parent process is still running

# Build directly from registry - NO Docker daemon involved        if ! kill -0 $$ 2>/dev/null; then

# docker:// prefix tells Singularity to pull from OCI registry            break

singularity build --force "$OUTPUT_PATH" "docker://${FULL_DOCKER_IMAGE}"        fi

SINGULARITY_EXIT=$?    done

) &

# Stop monitorMONITOR_PID=$!

kill $MONITOR_PID 2>/dev/null || true

wait $MONITOR_PID 2>/dev/null || true# Build directly from Docker daemon - no tar, no temp files

singularity build --force "$OUTPUT_PATH" "docker-daemon://${FULL_DOCKER_IMAGE}"

if [ $SINGULARITY_EXIT -ne 0 ]; thenSINGULARITY_EXIT=$?

    log "❌ Singularity build failed with exit code $SINGULARITY_EXIT"

    exit 1# Stop monitor

fikill $MONITOR_PID 2>/dev/null || true

wait $MONITOR_PID 2>/dev/null || true

END_TIME=$(date +%s)

DURATION=$((END_TIME - START_TIME))if [ $SINGULARITY_EXIT -ne 0 ]; then

DURATION_MIN=$((DURATION / 60))    log "❌ Singularity build failed with exit code $SINGULARITY_EXIT"

DURATION_SEC=$((DURATION % 60))    exit 1

SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)fi



log "============================================="END_TIME=$(date +%s)

log "✅ BUILD COMPLETE"DURATION=$((END_TIME - START_TIME))

log "============================================="DURATION_MIN=$((DURATION / 60))

log "📦 File: $OUTPUT_PATH"DURATION_SEC=$((DURATION % 60))

log "📊 Size: $SIZE"SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)

log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"

log ""log "============================================="

log "🧪 Test with:"log "✅ BUILD COMPLETE"

log "   singularity run $OUTPUT_PATH --help"log "============================================="

log ""log "📦 File: $OUTPUT_PATH"

log "🚀 Ready for deployment!"log "📊 Size: $SIZE"

log "⏱️  Time: ${DURATION_MIN}m ${DURATION_SEC}s"

# Quick verificationlog ""

if [ -f "$OUTPUT_PATH" ] && [ -s "$OUTPUT_PATH" ]; thenlog "🧪 Test with:"

    log "✅ SIF file created successfully"log "   singularity run $OUTPUT_PATH --help"

elselog ""

    log "❌ ERROR: SIF file not created or empty"log "🚀 Ready for deployment!"

    exit 1

fi# Quick verification

if [ -f "$OUTPUT_PATH" ] && [ -s "$OUTPUT_PATH" ]; then
    log "✅ SIF file created successfully"
else
    log "❌ ERROR: SIF file not created or empty"
    exit 1
fi