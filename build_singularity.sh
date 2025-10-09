#!/bin/bash
#
# Build Singularity SIF from micapipe Docker image
#
# This script converts the Docker image to Singularity SIF format
# and stores it in the MICA     if [ $BUILD_EXIT_CODE -eq 0 ]; then
        echo ""
        echo "📦 Moving SIF to final location..."
        mv "${TMP_OUTPUT}" "${OUTPUT_PATH}"
        echo "✅ Build complete"
    firity programs directory.
#
# Usage:
#   ./build_singularity.sh [docker_image_tag]
#
# Example:
#   ./build_singularity.sh                    # Uses :latest
#   ./build_singularity.sh v0.2.3-20251008    # Uses specific tag
#

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================
DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
SINGULARITY_DIR="/data_/mica1/01_programs/singularity"
OUTPUT_NAME="micapipe_v1_beta.sif"

# Alternate user location: /host/cassio/export03/data/enning/singularity
# You can override SINGULARITY_DIR by setting it before running:
#   SINGULARITY_DIR=/host/cassio/export03/data/enning/singularity ./build_singularity.sh

# Get Docker tag from argument or use 'latest'
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"
OUTPUT_PATH="${SINGULARITY_DIR}/${OUTPUT_NAME}"

# ============================================================================
# Header
# ============================================================================
echo "============================================="
echo "🔄 MICAPIPE SINGULARITY BUILD"
echo "============================================="
echo ""
echo "📦 Docker Image: ${FULL_DOCKER_IMAGE}"
echo "📍 Output Path:  ${OUTPUT_PATH}"
echo ""

# ============================================================================
# Pre-flight Checks
# ============================================================================
echo "🔍 Running pre-flight checks..."
echo ""

# Check if singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "❌ ERROR: singularity command not found"
    echo "   Please install Singularity/Apptainer first"
    exit 1
fi
echo "✅ Singularity found: $(singularity --version)"

# Check if Docker image exists
if ! docker image inspect "${FULL_DOCKER_IMAGE}" &> /dev/null; then
    echo "❌ ERROR: Docker image not found: ${FULL_DOCKER_IMAGE}"
    echo ""
    echo "Available micapipe images:"
    docker images "${DOCKER_IMAGE}" --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    exit 1
fi
echo "✅ Docker image found: ${FULL_DOCKER_IMAGE}"

# Get Docker image size
IMAGE_SIZE=$(docker image inspect "${FULL_DOCKER_IMAGE}" --format='{{.Size}}' | awk '{printf "%.1f GB\n", $1/1024/1024/1024}')
echo "   Size: ${IMAGE_SIZE}"

# Check if output directory exists, create if needed
if [ ! -d "${SINGULARITY_DIR}" ]; then
    echo "📁 Creating output directory: ${SINGULARITY_DIR}"
    mkdir -p "${SINGULARITY_DIR}"
fi
echo "✅ Output directory exists: ${SINGULARITY_DIR}"

# Check available disk space
AVAILABLE_SPACE=$(df -BG "${SINGULARITY_DIR}" | awk 'NR==2 {print $4}' | sed 's/G//')
echo "💾 Available disk space: ${AVAILABLE_SPACE} GB"
if [ "${AVAILABLE_SPACE}" -lt 150 ]; then
    echo "⚠️  WARNING: Low disk space (recommended: 150+ GB)"
    read -p "   Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
fi

# Check if output file already exists
if [ -f "${OUTPUT_PATH}" ]; then
    echo ""
    echo "⚠️  WARNING: Output file already exists: ${OUTPUT_PATH}"
    EXISTING_SIZE=$(du -h "${OUTPUT_PATH}" | cut -f1)
    echo "   Existing file size: ${EXISTING_SIZE}"
    read -p "   Overwrite? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
    echo "🗑️  Removing existing file..."
    rm -f "${OUTPUT_PATH}"
fi

echo ""
echo "============================================="
echo "🚀 STARTING SINGULARITY BUILD"
echo "============================================="
echo ""
echo "⏱️  This will take approximately 10-15 minutes"
echo "📊 Progress will be shown below..."
echo ""

# Set Singularity cache to data directory (not home directory)
# Use same base directory as output to avoid nodev mount issues
export SINGULARITY_CACHEDIR="${SINGULARITY_DIR}/.singularity_cache"
export SINGULARITY_TMPDIR="${SINGULARITY_DIR}/.singularity_tmp"

# Create cache directories if needed
mkdir -p "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"

echo "📦 Singularity cache: ${SINGULARITY_CACHEDIR}"
echo "📁 Singularity temp:  ${SINGULARITY_TMPDIR}"
echo ""

# Record start time
START_TIME=$(date +%s)

# ============================================================================
# Build Singularity SIF (using simple approach from CI test)
# ============================================================================

# Check if output filesystem has nodev - if so, build to /tmp first
MOUNT_INFO=$(mount | grep "$(df "${SINGULARITY_DIR}" | awk 'NR==2 {print $1}')" || true)
if echo "$MOUNT_INFO" | grep -q "nodev"; then
    echo "⚠️  WARNING: Output directory is on a 'nodev' mount"
    echo "   Building to /tmp first, then moving to final location..."
    echo ""
    
    # Build to data directory (not /tmp - insufficient space)
    DATA_BASE="/host/cassio/export03/data/enning"
    TMP_OUTPUT="${DATA_BASE}/.tmp_micapipe_build_$$.sif"
    export SINGULARITY_TMPDIR="${DATA_BASE}/.singularity_tmp"
    export SINGULARITY_CACHEDIR="${DATA_BASE}/.singularity_cache"
    mkdir -p "${SINGULARITY_TMPDIR}" "${SINGULARITY_CACHEDIR}"
    
    echo "📦 Temporary output: ${TMP_OUTPUT}"
    echo "📁 Temporary dirs:   ${DATA_BASE}/.singularity_{tmp,cache}"
    echo ""
    echo "🔧 Building SIF with --force flag..."
    
    singularity build --force \
        "${TMP_OUTPUT}" \
        "docker-daemon://${FULL_DOCKER_IMAGE}"
    
    BUILD_EXIT_CODE=$?
    
    if [ $BUILD_EXIT_CODE -eq 0 ]; then
        echo ""
        echo "� Moving SIF to final location..."
        mv "${TMP_OUTPUT}" "${OUTPUT_PATH}"
        echo "🧹 Cleaning up temporary directories..."
        rm -rf "${SINGULARITY_TMPDIR}" "${SINGULARITY_CACHEDIR}"
    fi
else
    echo "�🔧 Building SIF with --force flag..."
    singularity build --force \
        "${OUTPUT_PATH}" \
        "docker-daemon://${FULL_DOCKER_IMAGE}"
    
    BUILD_EXIT_CODE=$?
fi

# Record end time
END_TIME=$(date +%s)
BUILD_TIME=$((END_TIME - START_TIME))
BUILD_TIME_MIN=$((BUILD_TIME / 60))
BUILD_TIME_SEC=$((BUILD_TIME % 60))

echo ""
echo "============================================="

if [ ${BUILD_EXIT_CODE} -eq 0 ]; then
    echo "✅ SINGULARITY BUILD SUCCESSFUL!"
    echo "============================================="
    echo ""
    
    # Get SIF file info
    SIF_SIZE=$(du -h "${OUTPUT_PATH}" | cut -f1)
    
    echo "📍 Output File: ${OUTPUT_PATH}"
    echo "📦 File Size:   ${SIF_SIZE}"
    echo "⏱️  Build Time:  ${BUILD_TIME_MIN}m ${BUILD_TIME_SEC}s"
    echo ""
    
    # Verify SIF file
    echo "🔍 Verifying SIF file..."
    if singularity inspect "${OUTPUT_PATH}" &> /dev/null; then
        echo "✅ SIF file is valid"
        echo ""
        
        # Show SIF metadata
        echo "📋 Image Metadata:"
        singularity inspect --labels "${OUTPUT_PATH}" 2>/dev/null | head -10 || echo "   (metadata not available)"
        echo ""
    else
        echo "⚠️  Warning: Could not verify SIF file"
        echo ""
    fi
    
    echo "============================================="
    echo "🎯 NEXT STEPS"
    echo "============================================="
    echo ""
    echo "Test the image:"
    echo "   singularity run ${OUTPUT_PATH} --help"
    echo ""
    echo "Run micapipe with test data:"
    echo "   singularity run -B /data:/data ${OUTPUT_PATH} \\"
    echo "       -sub <subject> -out /data/derivatives -bids /data/bids"
    echo ""
    echo "Check image details:"
    echo "   singularity inspect ${OUTPUT_PATH}"
    echo ""
    
else
    echo "❌ SINGULARITY BUILD FAILED!"
    echo "============================================="
    echo ""
    echo "Build time: ${BUILD_TIME_MIN}m ${BUILD_TIME_SEC}s"
    echo "Exit code:  ${BUILD_EXIT_CODE}"
    echo ""
    echo "Common issues:"
    echo "  1. Insufficient disk space (need ~150 GB)"
    echo "  2. Docker daemon not running"
    echo "  3. Permissions issue on output directory"
    echo "  4. Corrupted Docker image (try rebuilding)"
    echo ""
    exit 1
fi
