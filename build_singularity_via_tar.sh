#!/bin/bash
#
# Build Singularity SIF from micapipe Docker image via tar export
#
# This method is more reliable for large images than docker-daemon://
# It exports the Docker image to a tar file, then builds from that.
#

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================
DOCKER_IMAGE="ghcr.io/mica-mni/micapipe"
DOCKER_TAG="${1:-latest}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"

SINGULARITY_DIR="/host/cassio/export03/data/enning/singularity"
OUTPUT_NAME="micapipe_v1_beta.sif"
OUTPUT_PATH="${SINGULARITY_DIR}/${OUTPUT_NAME}"

# Use data directory for intermediate tar file (will be deleted after)
TAR_FILE="/host/cassio/export03/data/enning/micapipe_docker_$$.tar"

# ============================================================================
# Header
# ============================================================================
echo "============================================="
echo "🔄 MICAPIPE SINGULARITY BUILD (via tar)"
echo "============================================="
echo ""
echo "📦 Docker Image: ${FULL_DOCKER_IMAGE}"
echo "📍 Output Path:  ${OUTPUT_PATH}"
echo "🗜️  Intermediate: ${TAR_FILE}"
echo ""

# ============================================================================
# Pre-flight Checks
# ============================================================================
echo "🔍 Running pre-flight checks..."
echo ""

# Check if Docker image exists
if ! docker image inspect "${FULL_DOCKER_IMAGE}" &>/dev/null; then
    echo "❌ ERROR: Docker image not found: ${FULL_DOCKER_IMAGE}"
    echo "   Available images:"
    docker images | grep micapipe || echo "   No micapipe images found"
    exit 1
fi
echo "✅ Docker image found: ${FULL_DOCKER_IMAGE}"

# Get image size
IMAGE_SIZE=$(docker image inspect "${FULL_DOCKER_IMAGE}" --format='{{.Size}}' | awk '{printf "%.1f GB\n", $1/1024/1024/1024}')
echo "   Size: ${IMAGE_SIZE}"

# Check /tmp space (need ~110GB for tar file)
DATA_DIR="/host/cassio/export03/data/enning"
TAR_SPACE=$(df -BG "${DATA_DIR}" | awk 'NR==2 {print $4}' | sed 's/G//')
echo ""
echo "💾 Available space in ${DATA_DIR}: ${TAR_SPACE} GB"
if [ "${TAR_SPACE}" -lt 120 ]; then
    echo "❌ ERROR: Not enough space in data directory (need 120+ GB for tar export)"
    echo "   Current: ${TAR_SPACE} GB"
    echo ""
    echo "💡 Options:"
    echo "   1. Clean up old files in ${DATA_DIR}"
    echo "   2. Use a different temp location (edit TAR_FILE variable)"
    exit 1
fi

# Check output directory space
mkdir -p "${SINGULARITY_DIR}"
OUTPUT_SPACE=$(df -BG "${SINGULARITY_DIR}" | awk 'NR==2 {print $4}' | sed 's/G//')
echo "💾 Available space in output dir: ${OUTPUT_SPACE} GB"
if [ "${OUTPUT_SPACE}" -lt 110 ]; then
    echo "⚠️  WARNING: Low disk space in output directory (recommended: 110+ GB)"
    read -p "   Continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
fi

# Check if output file exists
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

# ============================================================================
# Step 1: Export Docker image to tar
# ============================================================================
echo "📦 Step 1/3: Exporting Docker image to tar file..."
echo "   This will take 5-10 minutes for a 107GB image..."
echo ""

START_TIME=$(date +%s)

docker save "${FULL_DOCKER_IMAGE}" -o "${TAR_FILE}"

EXPORT_TIME=$(date +%s)
EXPORT_DURATION=$((EXPORT_TIME - START_TIME))
EXPORT_MIN=$((EXPORT_DURATION / 60))
EXPORT_SEC=$((EXPORT_DURATION % 60))

TAR_SIZE=$(du -h "${TAR_FILE}" | cut -f1)
echo ""
echo "✅ Docker image exported to tar: ${TAR_SIZE}"
echo "   Time: ${EXPORT_MIN}m ${EXPORT_SEC}s"
echo ""

# ============================================================================
# Step 2: Build Singularity SIF from tar
# ============================================================================
echo "🔧 Step 2/3: Building Singularity SIF from tar archive..."
echo "   This will take 10-15 minutes..."
echo ""

# Set Singularity cache/tmp to avoid nodev issues
export SINGULARITY_TMPDIR="/host/cassio/export03/data/enning/.singularity_tmp"
export SINGULARITY_CACHEDIR="/host/cassio/export03/data/enning/.singularity_cache"
mkdir -p "${SINGULARITY_TMPDIR}" "${SINGULARITY_CACHEDIR}"

# Build from tar archive (more reliable than docker-daemon://)
singularity build --force \
    "${OUTPUT_PATH}" \
    "docker-archive://${TAR_FILE}"

BUILD_EXIT_CODE=$?
BUILD_TIME=$(date +%s)
BUILD_DURATION=$((BUILD_TIME - EXPORT_TIME))
BUILD_MIN=$((BUILD_DURATION / 60))
BUILD_SEC=$((BUILD_DURATION % 60))

# ============================================================================
# Step 3: Cleanup
# ============================================================================
echo ""
echo "🧹 Step 3/3: Cleaning up temporary files..."

rm -f "${TAR_FILE}"

echo "✅ Cleanup complete"

# ============================================================================
# Results
# ============================================================================
TOTAL_TIME=$(date +%s)
TOTAL_DURATION=$((TOTAL_TIME - START_TIME))
TOTAL_MIN=$((TOTAL_DURATION / 60))
TOTAL_SEC=$((TOTAL_DURATION % 60))

echo ""
echo "============================================="
if [ $BUILD_EXIT_CODE -eq 0 ]; then
    echo "✅ SUCCESS"
    echo "============================================="
    echo ""
    SIF_SIZE=$(du -h "${OUTPUT_PATH}" | cut -f1)
    echo "📦 Singularity SIF created: ${OUTPUT_PATH}"
    echo "   Size: ${SIF_SIZE}"
    echo ""
    echo "⏱️  Timing:"
    echo "   Docker export: ${EXPORT_MIN}m ${EXPORT_SEC}s"
    echo "   SIF build:     ${BUILD_MIN}m ${BUILD_SEC}s"
    echo "   Total:         ${TOTAL_MIN}m ${TOTAL_SEC}s"
    echo ""
    echo "🧪 Test the SIF with:"
    echo "   singularity run ${OUTPUT_PATH} --help"
    echo ""
    echo "🧪 Or run full pipeline test:"
    echo "   ./test_micapipe.sh singularity ${OUTPUT_PATH}"
else
    echo "❌ BUILD FAILED"
    echo "============================================="
    echo ""
    echo "Exit code: $BUILD_EXIT_CODE"
    echo ""
    echo "💡 Troubleshooting:"
    echo "   - Check disk space: df -h /tmp && df -h ${SINGULARITY_DIR}"
    echo "   - Check permissions: ls -la ${SINGULARITY_DIR}"
    echo "   - Try building to different location"
    exit $BUILD_EXIT_CODE
fi
