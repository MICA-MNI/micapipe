#!/bin/bash
# ============================================================================
# BU#   ./build_main_image_server.sh
#
echo ""

# Verify we're in the right location
if [[ "$PWD" != *"/export03/data/enning/downloads"* ]]; then
    echo "⚠️  WARNING: Not running from expected server location"
    echo "   Current: $PWD"
    echo "   Expected: /export03/data/enning/downloads"MAIN IMAGE (STAGE 2) - SERVER VERSION
# ============================================================================
# This script builds the fast micapipe image using the pre-built base
# 
# PURPOSE: Create the final micapipe image by adding code to base
#   - Uses micapipe-base:latest as foundation
#   - Adds only micapipe source code and configuration
#   - Integrates with pre-installed tools
#
# WHEN TO RUN: Every time micapipe code changes (frequently!)
#!/bin/bash
# Build script for main micapipe image (adds micapipe code to base image)
# This script expects to be run on the server
set -e

echo "=================================="
echo "MICApipe Main Image Build"
echo "Server Build Script"
echo "=================================="
echo ""
echo "This script builds the main micapipe image by adding micapipe code"
echo "to the comprehensive base image."
echo ""
# WHERE TO RUN: On server at /export03/data/enning/downloads
echo "📁 Expected server location: /export03/data/enning/downloads"
echo ""
# USAGE:
#   cd /export03/data/enning/downloads
# EXPECTED TIME: 3-5 minutes (95% faster than full build!)
#
# USAGE:
#   cd /export03/data/enning/downloads
#   ./build_main_image_server.sh [OPTIONS]
#
# OPTIONS:
#   --base-image IMAGE   Base image to use (default: ghcr.io/mica-mni/micapipe-base:latest)
#   --tag TAG            Custom tag for main image (default: date-based)
#   --registry REG       Registry to push to (default: ghcr.io/mica-mni)
# ============================================================================

set -euo pipefail

# Disable Docker Content Trust to avoid SSL certificate issues
export DOCKER_CONTENT_TRUST=0

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_IMAGE="${MICAPIPE_BASE_IMAGE:-ghcr.io/mica-mni/micapipe-base:latest}"
CUSTOM_TAG=""
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
PUSH_TO_REGISTRY="false"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --base-image)
            BASE_IMAGE="$2"
            shift 2
            ;;
        --tag)
            CUSTOM_TAG="$2"
            shift 2
            ;;
        --registry)
            REGISTRY="$2"
            shift 2
            ;;
        --push)
            PUSH_TO_REGISTRY="true"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--base-image IMAGE] [--tag TAG] [--registry REG] [--push]"
            exit 1
            ;;
    esac
done

# ============================================================================
# VALIDATION
# ============================================================================

echo "🚀 MICApipe Main Image Builder (Stage 2 - Fast!)"
echo "================================================"
echo ""
echo "📍 Current directory: $PWD"
echo "🎯 Base image: $BASE_IMAGE"
echo "📦 Registry: $REGISTRY"
echo ""

# Check we're in the right location
if [[ "$PWD" != *"/export03/data/enning/downloads"* ]]; then
    echo "⚠️  Warning: Not in expected server downloads directory"
    echo "   Current: $PWD"
    echo "   Expected: /export03/data/enning/downloads"
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
fi

# Check for Dockerfile.main
if [[ ! -f "./Dockerfile.main" ]]; then
    echo "❌ Dockerfile.main not found in current directory"
    echo "   Please ensure you've copied the latest code from ~/micapipe"
    echo "   Run: ./migrate_comprehensive_base_to_server.sh"
    exit 1
fi

echo "✅ Dockerfile.main found"

# Check for required micapipe files
echo ""
echo "🔍 Checking for micapipe source files..."
REQUIRED_MICAPIPE_FILES=(
    "micapipe"
    "functions"
    "fix_settings.sh"
    "fsl_conf"
    "R_config"
    "surfaces"
    "parcellations"
)

MISSING_FILES=()
for file in "${REQUIRED_MICAPIPE_FILES[@]}"; do
    if [[ -e "$file" ]]; then
        echo "   ✅ $file"
    else
        echo "   ❌ $file: NOT FOUND"
        MISSING_FILES+=("$file")
    fi
done

if [[ ${#MISSING_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "❌ Required micapipe files are missing!"
    echo "   Please run: ./migrate_comprehensive_base_to_server.sh"
    exit 1
fi

echo "✅ All required micapipe files found"

# Check Docker access
if ! docker info >/dev/null 2>&1; then
    echo "❌ Cannot access Docker daemon"
    exit 1
fi

echo "✅ Docker access confirmed"

# Check if base image exists
echo ""
echo "🔍 Checking if base image exists..."
if docker image inspect "$BASE_IMAGE" >/dev/null 2>&1; then
    echo "✅ Base image found locally: $BASE_IMAGE"
    BASE_IMAGE_SIZE=$(docker images "$BASE_IMAGE" --format "{{.Size}}")
    echo "   Size: $BASE_IMAGE_SIZE"
elif docker pull "$BASE_IMAGE" 2>/dev/null; then
    echo "✅ Base image pulled from registry: $BASE_IMAGE"
else
    echo ""
    echo "❌ Base image not found: $BASE_IMAGE"
    echo ""
    echo "You need to build the base image first!"
    echo "Run: ./build_base_image_server.sh"
    echo ""
    echo "Or specify a different base image:"
    echo "   ./build_main_image_server.sh --base-image micapipe-base:latest"
    echo ""
    exit 1
fi

# ============================================================================
# BUILD CONFIGURATION
# ============================================================================

# Generate tag
if [[ -z "$CUSTOM_TAG" ]]; then
    IMAGE_TAG="v0.2.3-$(date +%Y%m%d)"
else
    IMAGE_TAG="$CUSTOM_TAG"
fi

MAIN_IMAGE_NAME="micapipe"
FULL_MAIN_IMAGE="${REGISTRY}/${MAIN_IMAGE_NAME}:${IMAGE_TAG}"
LATEST_MAIN_IMAGE="${REGISTRY}/${MAIN_IMAGE_NAME}:latest"

# Build log
BUILD_LOG="build_main_${IMAGE_TAG}_$(date +%H%M%S).log"

echo ""
echo "🚀 Build Configuration"
echo "====================="
echo "🎯 Using base: ${BASE_IMAGE}"
echo "📦 Main image: ${FULL_MAIN_IMAGE}"
echo "🏷️  Latest tag: ${LATEST_MAIN_IMAGE}"
echo "📝 Build log: ${BUILD_LOG}"
echo "⏱️  Expected time: 3-5 minutes"
echo ""
echo "📋 This build adds:"
echo "   ✅ MICApipe source code"
echo "   ✅ Processing scripts and functions"
echo "   ✅ Configuration files (fix_settings, fsl_conf)"
echo "   ✅ Surface data and parcellations"
echo "   ✅ R configuration"
echo ""

read -p "Proceed with main image build? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ Build cancelled"
    exit 0
fi

# ============================================================================
# DOCKER BUILD
# ============================================================================

echo ""
echo "🏗️  Starting Docker build..."
echo "============================================="

# Build arguments
BUILD_ARGS=(
    --file Dockerfile.main
    --build-arg "BASE_IMAGE=${BASE_IMAGE}"
    --tag "${FULL_MAIN_IMAGE}"
    --tag "${LATEST_MAIN_IMAGE}"
)

# Start build
START_TIME=$(date +%s)

if docker build "${BUILD_ARGS[@]}" . 2>&1 | tee "$BUILD_LOG"; then
    
    END_TIME=$(date +%s)
    BUILD_DURATION=$((END_TIME - START_TIME))
    BUILD_MINUTES=$((BUILD_DURATION / 60))
    BUILD_SECONDS=$((BUILD_DURATION % 60))
    
    echo ""
    echo "🎉 MAIN IMAGE BUILD SUCCESSFUL!"
    echo "============================================="
    echo "✅ Image: ${FULL_MAIN_IMAGE}"
    echo "✅ Latest: ${LATEST_MAIN_IMAGE}"
    echo "📄 Log: ${BUILD_LOG}"
    echo "⏱️  Build time: ${BUILD_MINUTES}m ${BUILD_SECONDS}s"
    echo ""
    echo "📊 Image Details:"
    docker images "${REGISTRY}/${MAIN_IMAGE_NAME}" --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    echo ""
    
    # Show speed improvement
    if [[ $BUILD_MINUTES -lt 10 ]]; then
        echo "🚀 FAST BUILD SUCCESS!"
        echo "   Traditional build: ~60-90 minutes"
        echo "   This build: ${BUILD_MINUTES}m ${BUILD_SECONDS}s"
        echo "   Time saved: ~$((60 - BUILD_MINUTES)) minutes!"
        echo ""
    fi
    
    # Offer to push
    if [[ "$PUSH_TO_REGISTRY" == "true" ]]; then
        echo "📤 Pushing to registry..."
        docker push "${FULL_MAIN_IMAGE}"
        docker push "${LATEST_MAIN_IMAGE}"
        echo "✅ Push complete!"
        echo ""
    else
        echo "💡 To push this image to registry, run:"
        echo "   docker push ${FULL_MAIN_IMAGE}"
        echo "   docker push ${LATEST_MAIN_IMAGE}"
        echo ""
    fi
    
    echo "🎯 IMAGE READY TO USE"
    echo "============================================="
    echo "Test the image locally:"
    echo "   docker run --rm ${LATEST_MAIN_IMAGE} --help"
    echo ""
    echo "Or convert to Singularity:"
    echo "   singularity build micapipe.sif docker-daemon://${LATEST_MAIN_IMAGE}"
    echo ""
    echo "💡 Future code changes can rebuild in ${BUILD_MINUTES}-$((BUILD_MINUTES+2)) minutes!"
    echo ""
    
else
    echo ""
    echo "❌ MAIN IMAGE BUILD FAILED!"
    echo "============================================="
    echo "📄 Check build log: ${BUILD_LOG}"
    echo ""
    echo "Common issues:"
    echo "  1. Base image not accessible"
    echo "  2. Missing micapipe source files"
    echo "  3. Permission issues with files"
    echo ""
    exit 1
fi
