#!/bin/bash
# ============================================================================
# BUILD MICAPIPE BASE IMAGE (STAGE 1) - SERVER VERSION
# ============================================================================
# This script builds the comprehensive base image with ALL neuroimaging tools
# 
# PURPOSE: Create a stable base image containing:
#   - FSL 6.0.2, FreeSurfer 7.4.1, AFNI, ANTs 2.3.4
#   - Miniconda/Mamba with Python environments
#   - MRtrix3 3.0.7, FastSurfer 2.4.2
#   - DESIGNER, LAMAReg, SWM, Synb0/SynBOLD
#   - R environment, c3d tools
#
# WHEN TO RUN: Only when dependencies change (rare - maybe monthly)
# WHERE TO RUN: On server at /host/cassio/export03/data/enning/downloads
# EXPECTED TIME: 45-90 minutes (one-time cost for fast future builds)
#
# USAGE:
#   cd /host/cassio/export03/data/enning/downloads
#   ./build_base_image_server.sh [OPTIONS]
#
# OPTIONS:
#   --enable-cuda    Build with CUDA support (default: CPU-only)
#   --tag TAG        Custom tag for image (default: date-based)
#   --registry REG   Registry to push to (default: ghcr.io/mica-mni)
# ============================================================================

set -euo pipefail

# Disable Docker Content Trust to avoid SSL certificate issues
export DOCKER_CONTENT_TRUST=0

# ============================================================================
# CONFIGURATION
# ============================================================================

ENABLE_CUDA="false"
CUSTOM_TAG=""
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
PUSH_TO_REGISTRY="false"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --enable-cuda)
            ENABLE_CUDA="true"
            shift
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
            echo "Usage: $0 [--enable-cuda] [--tag TAG] [--registry REG] [--push]"
            exit 1
            ;;
    esac
done

# ============================================================================
# VALIDATION
# ============================================================================

echo "🏗️  MICApipe Base Image Builder (Stage 1)"
echo "========================================="
echo ""
echo "📍 Current directory: $PWD"
echo "🎮 CUDA support: $ENABLE_CUDA"
echo "📦 Registry: $REGISTRY"
echo ""

# Check we're in the right location
if [[ "$PWD" != *"/host/cassio/export03/data/enning/downloads"* ]]; then
    echo "⚠️  Warning: Not in expected server downloads directory"
    echo "   Current: $PWD"
    echo "   Expected: /host/cassio/export03/data/enning/downloads"
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
fi

# Check for Dockerfile.base
if [[ ! -f "./Dockerfile.base" ]]; then
    echo "❌ Dockerfile.base not found in current directory"
    echo "   Please ensure you've copied the latest code from ~/micapipe"
    echo "   Run: ./migrate_comprehensive_base_to_server.sh"
    exit 1
fi

echo "✅ Dockerfile.base found"

# Check for pre-downloaded files
echo ""
echo "🔍 Checking for pre-downloaded files..."
REQUIRED_FILES=(
    "fsl-6.0.2-centos6_64.tar.gz"
    "freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
    "afni-linux_openmp_64.tgz"
    "fix-1.068.tar.gz"
    "Miniconda3-py39_22.11.1-1-Linux-x86_64.sh"
)

FILES_FOUND=0
MISSING_FILES=()

for file in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        FILE_SIZE=$(du -h "$file" | cut -f1)
        echo "   ✅ $file: $FILE_SIZE"
        FILES_FOUND=$((FILES_FOUND + 1))
    else
        echo "   ❌ $file: NOT FOUND"
        MISSING_FILES+=("$file")
    fi
done

if [[ ${#MISSING_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "⚠️  Warning: Some pre-downloaded files are missing"
    echo "   Found: $FILES_FOUND/${#REQUIRED_FILES[@]}"
    echo "   Missing files will be downloaded during build (slower, requires internet)"
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
else
    echo ""
    echo "✅ All pre-downloaded files found!"
fi

# Check Docker access
if ! docker info >/dev/null 2>&1; then
    echo "❌ Cannot access Docker daemon"
    echo "   Please ensure Docker is running and you have permissions"
    exit 1
fi

echo "✅ Docker access confirmed"

# ============================================================================
# BUILD CONFIGURATION
# ============================================================================

# Generate tag
if [[ -z "$CUSTOM_TAG" ]]; then
    DATE_TAG=$(date +%Y%m%d)
    if [[ "$ENABLE_CUDA" == "true" ]]; then
        IMAGE_TAG="${DATE_TAG}-cuda"
    else
        IMAGE_TAG="${DATE_TAG}"
    fi
else
    IMAGE_TAG="$CUSTOM_TAG"
fi

BASE_IMAGE_NAME="micapipe-base"
FULL_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:${IMAGE_TAG}"
LATEST_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:latest"

# Build log
BUILD_LOG="build_base_${IMAGE_TAG}_$(date +%H%M%S).log"

echo ""
echo "🚀 Build Configuration"
echo "====================="
echo "📦 Base image: ${FULL_BASE_IMAGE}"
echo "🏷️  Latest tag: ${LATEST_BASE_IMAGE}"
echo "📝 Build log: ${BUILD_LOG}"
echo "⏱️  Expected time: 45-90 minutes"
echo "💾 Using temp directory: /host/cassio/export03/data/enning"
echo ""

if [[ "$ENABLE_CUDA" == "true" ]]; then
    echo "🎮 CUDA-enabled build:"
    echo "   - CUDA 11.8 runtime libraries"
    echo "   - GPU-accelerated PyTorch/TensorFlow"
    echo "   - GPU-enabled FastSurfer"
else
    echo "💻 CPU-only build:"
    echo "   - No CUDA dependencies"
    echo "   - CPU-only PyTorch/TensorFlow"
    echo "   - CPU-only FastSurfer"
fi

echo ""
echo "📋 This base image includes:"
echo "   ✅ FSL 6.0.2 (frozen)"
echo "   ✅ FreeSurfer 7.4.1 (frozen)"
echo "   ✅ FastSurfer 2.4.2 (frozen)"
echo "   ✅ AFNI (latest)"
echo "   ✅ ANTs 2.3.4"
echo "   ✅ MRtrix3 3.0.7 (frozen)"
echo "   ✅ DESIGNER v2 pipeline"
echo "   ✅ LAMAReg + ANTsPy"
echo "   ✅ SWM (Superficial White Matter)"
echo "   ✅ Synb0-DISCO + SynBOLD-DisCo"
echo "   ✅ Miniconda/Mamba + Python environments"
echo "   ✅ R environment with neuroimaging packages"
echo "   ✅ FSL FIX, Connectome Workbench, c3d"
echo ""

read -p "Proceed with base image build? (y/N): " -n 1 -r
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
    --file Dockerfile.base
    --build-arg "ENABLE_CUDA=${ENABLE_CUDA}"
    --build-arg "CUSTOM_TMPDIR=/host/cassio/export03/data/enning"
    --build-arg "DOWNLOADS_DIR=/downloads"
    --tag "${FULL_BASE_IMAGE}"
    --tag "${LATEST_BASE_IMAGE}"
)

# Memory limits for server builds
BUILD_ARGS+=(
    --memory=12g
    --memory-swap=16g
)

# Start build
if docker build "${BUILD_ARGS[@]}" . 2>&1 | tee "$BUILD_LOG"; then
    
    echo ""
    echo "🎉 BASE IMAGE BUILD SUCCESSFUL!"
    echo "============================================="
    echo "✅ Image: ${FULL_BASE_IMAGE}"
    echo "✅ Latest: ${LATEST_BASE_IMAGE}"
    echo "📄 Log: ${BUILD_LOG}"
    echo ""
    echo "📊 Image Details:"
    docker images "${REGISTRY}/${BASE_IMAGE_NAME}" --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
    echo ""
    
    # Offer to push
    if [[ "$PUSH_TO_REGISTRY" == "true" ]]; then
        echo "📤 Pushing to registry..."
        docker push "${FULL_BASE_IMAGE}"
        docker push "${LATEST_BASE_IMAGE}"
        echo "✅ Push complete!"
        echo ""
    else
        echo "💡 To push this image to registry, run:"
        echo "   docker push ${FULL_BASE_IMAGE}"
        echo "   docker push ${LATEST_BASE_IMAGE}"
        echo ""
    fi
    
    echo "🎯 NEXT STEP: Build fast main image"
    echo "============================================="
    echo "Now you can build the main micapipe image quickly using:"
    echo "   ./build_main_image_server.sh"
    echo ""
    echo "Main image builds will be FAST (3-5 min) since they use this base!"
    echo ""
    
else
    echo ""
    echo "❌ BASE IMAGE BUILD FAILED!"
    echo "============================================="
    echo "📄 Check build log: ${BUILD_LOG}"
    echo ""
    echo "Common issues:"
    echo "  1. Network timeout downloading packages"
    echo "  2. Out of memory (try freeing space)"
    echo "  3. Missing pre-downloaded files"
    echo ""
    exit 1
fi
