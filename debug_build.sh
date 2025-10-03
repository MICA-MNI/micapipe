#!/bin/bash
#
# Debug Docker Build Script
# Helps troubleshoot hanging Docker builds with detailed logging
#

set -euo pipefail

echo "🔍 Docker Build Debugging Tool"
echo "=============================="

# Check if we're in the right directory
if [[ ! -f "Dockerfile" ]]; then
    echo "❌ No Dockerfile found in current directory"
    echo "   Please run this from the build directory"
    exit 1
fi

echo "📋 System Information:"
echo "Docker version: $(docker --version)"
echo "Available memory: $(free -h | grep Mem | awk '{print $7}')"
echo "Available disk: $(df -h . | tail -1 | awk '{print $4}')"
echo ""

echo "🧹 Cleaning up Docker resources..."
# Stop any running containers that might be hanging
docker ps -q | xargs -r docker kill || true

# Remove dangling containers
docker container prune -f || true

# Remove dangling images to free up space
docker image prune -f || true

echo "✅ Cleanup complete"
echo ""

echo "🚀 Starting Docker build with debugging..."
echo "   Using --no-cache to avoid cached layer issues"
echo "   Using --progress=plain for detailed output"
echo ""

# Build with detailed output and no cache
docker build \
    --progress=plain \
    --no-cache \
    --memory=12g \
    --memory-swap=16g \
    --tag micapipe:debug \
    . 2>&1 | tee build_debug.log

BUILD_EXIT_CODE=${PIPESTATUS[0]}

if [[ $BUILD_EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "✅ Debug build completed successfully!"
    echo "🏷️  Tagged as: micapipe:debug"
else
    echo ""
    echo "❌ Debug build failed (exit code: $BUILD_EXIT_CODE)"
    echo "📋 Full log saved to: build_debug.log"
    echo "📋 Last 50 lines of output:"
    tail -50 build_debug.log
fi

exit $BUILD_EXIT_CODE