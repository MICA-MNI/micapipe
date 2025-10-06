#!/bin/bash
# ============================================================================
# SIMPLE COMPREHENSIVE BASE BUILD (No bullshit, just works)
# Based on working docker-container-updates branch approach
# ============================================================================

set -euo pipefail

echo "🐳 Simple MICApipe Comprehensive Base Builder"
echo "============================================="

# Check we're in the right location
if [[ ! -f "./Dockerfile.mamba-base" ]]; then
    echo "❌ Dockerfile.mamba-base not found"
    echo "   Run from: /host/cassio/export03/data/enning/downloads"
    exit 1
fi

echo "📍 Build directory: $(pwd)"

# Check for pre-downloaded files (no copying needed - they're already here!)
echo "📦 Checking pre-downloaded files in current directory..."
FILES_FOUND=0

if [[ -f "./freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]]; then
    echo "   ✅ FreeSurfer: $(du -h ./freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz | cut -f1)"
    FILES_FOUND=1
fi

if [[ -f "./fsl-6.0.2-centos6_64.tar.gz" ]]; then
    echo "   ✅ FSL: $(du -h ./fsl-6.0.2-centos6_64.tar.gz | cut -f1)"
    FILES_FOUND=1
fi

if [[ -f "./afni-linux_openmp_64.tgz" ]]; then
    echo "   ✅ AFNI: $(du -h ./afni-linux_openmp_64.tgz | cut -f1)"
    FILES_FOUND=1
fi

if [[ $FILES_FOUND -eq 0 ]]; then
    echo "   ⚠️  No pre-downloaded files found - will download during build"
fi

# Check Docker access
if ! docker info >/dev/null 2>&1; then
    echo "❌ Cannot access Docker daemon"
    exit 1
fi

echo "✅ Docker access confirmed"

# Build configuration
BUILD_LOG="simple_base_build_$(date +%Y%m%d_%H%M%S).log"
IMAGE_TAG="micapipe-comprehensive-base:$(date +%Y%m%d)"

echo ""
echo "🚀 Starting Docker build..."
echo "📋 Image: $IMAGE_TAG"
echo "📝 Log: $BUILD_LOG"
echo "⏱️  Expected time: 45-90 minutes"
echo ""

# Simple Docker build - no fancy stuff, just build it!
if docker build \
    --file Dockerfile.mamba-base \
    --tag "$IMAGE_TAG" \
    --tag "micapipe-comprehensive-base:latest" \
    --build-arg CUSTOM_TMPDIR="/host/cassio/export03/data/enning" \
    . 2>&1 | tee "$BUILD_LOG"; then
    
    echo ""
    echo "🎉 SUCCESS! Base image built successfully!"
    echo "=================================="
    echo "📋 Image: $IMAGE_TAG"
    echo "📋 Latest: micapipe-comprehensive-base:latest"
    echo "📄 Log: $BUILD_LOG"
    echo ""
    echo "🎯 Next step: Build fast main image"
    echo "   docker build -f Dockerfile.minimal -t micapipe:latest ."
    
else
    echo ""
    echo "❌ BUILD FAILED!"
    echo "Check log: $BUILD_LOG"
    exit 1
fi