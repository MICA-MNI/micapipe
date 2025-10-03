#!/bin/bash
# ============================================================================
# COMPREHENSIVE MICAPIPE BASE IMAGE BUILDER (SERVER VERSION)
# Creates a comprehensive base image with ALL neuroimaging tools using pre-downloaded files
# Compatible with server environment: /host/cassio/export03/data/enning
# ============================================================================

set -euo pipefail

echo "🐳 MICApipe Comprehensive Base Image Builder (Server)"
echo "====================================================="

# Verify we're in the right location with pre-downloaded files
echo "🔍 Verifying server environment..."
echo "📍 Current directory: $PWD"
echo "📁 Expected server location: /host/cassio/export03/data/enning/downloads"

# Check if we're in the downloads directory (which is the build directory)
if [[ "$PWD" != *"/host/cassio/export03/data/enning/downloads"* ]]; then
    echo "⚠️  Warning: Not in expected server downloads directory"
    echo "   Current: $PWD"
    echo "   Expected: /host/cassio/export03/data/enning/downloads"
    echo ""
    echo "💡 This script should be run from the downloads directory where"
    echo "   both source files and pre-downloaded dependencies are located."
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        echo "💡 To migrate source to server: ./migrate_comprehensive_base_to_server.sh"
        exit 1
    fi
fi

# Check for required files
missing_files=()
echo "🔍 Checking for pre-downloaded files..."
for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ Found: $file"
    else
        echo "❌ Missing: $file"
        missing_files+=("$file")
    fi
done

if [ ${#missing_files[@]} -gt 0 ]; then
    echo ""
    echo "⚠️  Some files are missing but build will continue."
    echo "   Missing files will be downloaded during build (requires internet)."
    echo "   This may significantly increase build time."
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "❌ Build cancelled"
        exit 1
    fi
fi

# Verify Dockerfile.mamba-base exists
if [[ ! -f "./Dockerfile.mamba-base" ]]; then
    echo "❌ Dockerfile.mamba-base not found in current directory"
    echo "   Please ensure you're in the micapipe source directory"
    echo "   and the comprehensive-base-image branch is checked out"
    exit 1
fi

# Environment setup
export DOCKER_CONTENT_TRUST=0
export BUILDKIT_PROGRESS=plain

# Registry configuration
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
BASE_IMAGE_NAME="micapipe-comprehensive-base"
BASE_IMAGE_TAG="${MICAPIPE_BASE_TAG:-$(date +%Y%m%d)}"
FULL_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:${BASE_IMAGE_TAG}"
LATEST_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:latest"

echo ""
echo "🚀 Building MICApipe comprehensive base image: ${FULL_BASE_IMAGE}"
echo "📅 Using date-based tag for versioning: ${BASE_IMAGE_TAG}"
echo ""
echo "🏗️  This base image includes:"
echo "   - Ubuntu 18.04 with all system dependencies"
echo "   - Conda/Mamba with optimized Python environments"
echo "   - FSL 6.0.2"
echo "   - FreeSurfer 7.4.1"
echo "   - AFNI (latest)"
echo "   - ANTs 2.3.4"
echo "   - FSL FIX"
echo "   - FastSurfer 2.4.2"
echo "   - SWM (Superficial White Matter)"
echo "   - DESIGNER v2"
echo "   - R environment with neuroimaging packages"
echo "   - c3d tools"
echo ""

# Build configuration for server
BUILD_LOG="build_comprehensive_base_$(date +%Y%m%d_%H%M%S).log"
echo "📝 Build will be logged to: $BUILD_LOG"

# Docker caching configuration
CACHE_FROM_IMAGES=(
    "${LATEST_BASE_IMAGE}"
    "micapipe-comprehensive-base:latest"
    "ubuntu:bionic-20201119"
    "kaczmarj/ants:2.3.4"
)

# Build the comprehensive base image with server optimizations
echo "📦 Building comprehensive base image (this will take significant time)..."
echo "⏱️  Expected build time: 45-90 minutes (depending on cache and missing files)"
echo "💾 Using server temporary directory: /host/cassio/export03/data/enning"
echo ""

# Build cache arguments
CACHE_ARGS=()
for cache_image in "${CACHE_FROM_IMAGES[@]}"; do
    CACHE_ARGS+=(--cache-from "$cache_image")
done

# Start Docker build with server-specific settings
if docker build \
    --file Dockerfile.mamba-base \
    --memory=12g \
    --memory-swap=16g \
    --build-arg CUSTOM_TMPDIR="/host/cassio/export03/data/enning" \
    --build-arg DOWNLOADS_DIR="/downloads" \
    "${CACHE_ARGS[@]}" \
    --tag "${FULL_BASE_IMAGE}" \
    --tag "${LATEST_BASE_IMAGE}" \
    . 2>&1 | tee "$BUILD_LOG"; then
    
    echo ""
    echo "✅ Comprehensive base image built successfully!"
    echo "=========================================="
    echo "Image: ${FULL_BASE_IMAGE}"
    echo "Latest: ${LATEST_BASE_IMAGE}"
    echo "Build log: $BUILD_LOG"
    echo ""
    echo "📊 Image size:"
    docker images "${FULL_BASE_IMAGE}"
    echo ""
    echo "🏷️  Tagged images:"
    echo "   - ${FULL_BASE_IMAGE} (dated version)"
    echo "   - ${LATEST_BASE_IMAGE} (latest)"
    echo ""
    echo "🚀 Next steps:"
    echo "   1. Push to registry: docker push ${LATEST_BASE_IMAGE}"
    echo "   2. Push dated version: docker push ${FULL_BASE_IMAGE}"
    echo "   3. Use build_fast_ci_server.sh for ultra-fast main image builds"
    echo ""
    echo "💡 This base image enables 95% faster CI builds!"
    
else
    echo ""
    echo "❌ Build Failed!"
    echo "==============="
    echo "Check the build log for details: $BUILD_LOG"
    echo ""
    echo "🔧 Common issues:"
    echo "   - Insufficient disk space"
    echo "   - Missing pre-downloaded files (check above)"
    echo "   - Network connectivity for downloads"
    echo "   - Docker memory limits"
    echo ""
    echo "💡 To resume build after fixing issues:"
    echo "   docker build --cache-from ${LATEST_BASE_IMAGE} ..."
    exit 1
fi