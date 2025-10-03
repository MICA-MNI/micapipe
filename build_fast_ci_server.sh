#!/bin/bash
# ============================================================================
# FAST MICAPIPE CI BUILD SCRIPT (SERVER VERSION)
# Uses comprehensive base image for ultra-fast builds in server environment
# Compatible with: /host/cassio/export03/data/enning
# ============================================================================

set -euo pipefail

echo "🚀 Fast MICApipe CI Build (Server)"
echo "=================================="

# Configuration
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
BASE_IMAGE="micapipe-comprehensive-base:latest"
FULL_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE}"
MAIN_IMAGE_NAME="${MICAPIPE_IMAGE_NAME:-micapipe}"
MAIN_IMAGE_TAG="${MICAPIPE_TAG:-latest}"
MAIN_IMAGE="${REGISTRY}/${MAIN_IMAGE_NAME}:${MAIN_IMAGE_TAG}"

echo "Base image: ${FULL_BASE_IMAGE}"
echo "Target image: ${MAIN_IMAGE}"
echo ""

# Check if we're in the right server location
if [[ "$PWD" != *"/host/cassio/export03/data/enning"* ]]; then
    echo "⚠️  Warning: Not in expected server directory"
    echo "   Current: $PWD"
    echo "   Expected: /host/cassio/export03/data/enning/[micapipe_directory]"
    echo ""
fi

# Verify Dockerfile.minimal exists
if [[ ! -f "./Dockerfile.minimal" ]]; then
    echo "❌ Dockerfile.minimal not found in current directory"
    echo "   Please ensure you're in the micapipe source directory"
    echo "   and the comprehensive-base-image branch is checked out"
    exit 1
fi

# Check if base image is available locally
echo "🔍 Checking for base image..."
if docker images | grep -q "micapipe-comprehensive-base"; then
    echo "✅ Base image found locally"
elif docker pull "${FULL_BASE_IMAGE}" 2>/dev/null; then
    echo "✅ Base image pulled from registry"
else
    echo "❌ Base image not available!"
    echo ""
    echo "🔧 Please ensure the base image is available by either:"
    echo "   1. Building it locally: ./build_comprehensive_base_server.sh"
    echo "   2. Pulling from registry: docker pull ${FULL_BASE_IMAGE}"
    echo "   3. Using a different base image with MICAPIPE_REGISTRY env var"
    exit 1
fi

echo ""
echo "🏗️  Building main MICApipe image (should be very fast)..."
echo "⏱️  Expected build time: 3-5 minutes"
echo "💾 Using server temporary directory: /host/cassio/export03/data/enning"

# Build configuration
BUILD_LOG="build_fast_ci_$(date +%Y%m%d_%H%M%S).log"
echo "📝 Build will be logged to: $BUILD_LOG"

# Build main image using minimal Dockerfile with server settings
if docker build \
    --file Dockerfile.minimal \
    --memory=8g \
    --memory-swap=10g \
    --build-arg BASE_IMAGE="${FULL_BASE_IMAGE}" \
    --build-arg CUSTOM_TMPDIR="/host/cassio/export03/data/enning" \
    --cache-from "${MAIN_IMAGE}" \
    --tag "${MAIN_IMAGE}" \
    . 2>&1 | tee "$BUILD_LOG"; then

    echo ""
    echo "✅ Build completed successfully!"
    echo "==============================="
    echo "Final image: ${MAIN_IMAGE}"
    echo "Build log: $BUILD_LOG"
    echo ""
    echo "📊 Final image size:"
    docker images "${MAIN_IMAGE}"
    echo ""
    echo "🎯 Performance Summary:"
    echo "   - Base image (contains): FSL, FreeSurfer, AFNI, ANTs, Python environments"
    echo "   - Main image (contains): MICApipe code, R packages, configuration"  
    echo "   - Build time: ~3-5 minutes (vs 45-90 minutes full build)"
    echo "   - Efficiency gain: 93-95% faster"
    echo ""
    echo "🚀 To push to registry:"
    echo "   docker push ${MAIN_IMAGE}"
    echo ""
    echo "🧪 To test the image:"
    echo "   docker run --rm -it ${MAIN_IMAGE} /bin/bash"
    
else
    echo ""
    echo "❌ Fast CI Build Failed!"
    echo "======================="
    echo "Check the build log for details: $BUILD_LOG"
    echo ""
    echo "🔧 Common issues:"
    echo "   - Base image not available (see error above)"
    echo "   - Insufficient disk space for final image"
    echo "   - Issues with MICApipe source code"
    echo "   - R package installation failures"
    exit 1
fi