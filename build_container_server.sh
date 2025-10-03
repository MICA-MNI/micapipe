#!/bin/bash
set -euo pipefail

# MICApipe Container Builder (Server Version)
# ===========================================
# This script runs Docker build from the server location with pre-downloaded dependencies

echo "🐳 MICApipe Container Builder (Server)"
echo "======================================"

# Verify we're in the right location
if [[ ! -f "./Dockerfile" ]]; then
    echo "❌ Dockerfile not found in current directory"
    echo "   Please run this script from the build directory:"
    echo "   cd /host/cassio/export03/data/enning/micapipe_build"
    echo "   ./build_container_server.sh"
    exit 1
fi

# Check if we're in the downloads directory (which is now the build directory)
if [[ "$PWD" != *"/host/cassio/export03/data/enning/downloads"* ]]; then
    echo "⚠️  Warning: Not in expected downloads/build directory"
    echo "   Current: $PWD"
    echo "   Expected: /host/cassio/export03/data/enning/downloads"
    echo ""
fi

# Environment setup
export DOCKER_CONTENT_TRUST=0
export BUILDKIT_PROGRESS=plain

# Docker caching configuration
CACHE_FROM_IMAGES=(
    "micapipe:latest"
    "micapipe:cache"
    "ubuntu:xenial-20161213"
    "kaczmarj/ants:2.3.4"
)

# Build configuration
USE_CACHE=true
FORCE_REBUILD=false
CACHE_TAG="micapipe:cache"

# Parse command line arguments for caching options
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-cache)
            USE_CACHE=false
            echo "🚫 Docker caching disabled"
            shift
            ;;
        --force-rebuild)
            FORCE_REBUILD=true
            echo "🔄 Force rebuild enabled (ignoring all cache)"
            shift
            ;;
        --cache-tag)
            CACHE_TAG="$2"
            echo "🏷️  Using cache tag: $CACHE_TAG"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --no-cache       Disable Docker layer caching"
            echo "  --force-rebuild  Force complete rebuild (ignore all cache)"
            echo "  --cache-tag TAG  Use specific cache tag (default: micapipe:cache)"
            echo "  --help, -h       Show this help"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check Docker access
if ! docker info >/dev/null 2>&1; then
    echo "❌ Cannot access Docker daemon. Solutions:"
    echo "1. Ask admin to add you to docker group: sudo usermod -aG docker $(whoami)"
    echo "2. Use Docker Desktop if available"
    echo "3. Contact system administrator"
    exit 1
fi

echo "✅ Docker access confirmed"

# Check for downloads (should be in current directory since build dir = downloads dir)
echo "📦 Checking downloads in current directory..."
if [[ -f "./fsl-6.0.2-centos6_64.tar.gz" ]]; then
    FSL_SIZE=$(du -h "./fsl-6.0.2-centos6_64.tar.gz" | cut -f1)
    echo "   ✅ FSL: $FSL_SIZE"
else
    echo "   ❌ FSL not found in current directory"
    echo "   Please run ./migrate_to_server.sh first"
    exit 1
fi

if [[ -f "./freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]]; then
    FS_SIZE=$(du -h "./freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" | cut -f1)
    echo "   ✅ FreeSurfer: $FS_SIZE"
else
    echo "   ❌ FreeSurfer not found in current directory"
    echo "   Please run ./migrate_to_server.sh first"
    exit 1
fi

echo "✅ All pre-downloaded files ready for Docker build"

# Create build logs directory
mkdir -p build_logs

# Generate log filename
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BUILD_LOG="build_logs/container_build_${TIMESTAMP}.log"

echo ""
echo "🚀 Starting Docker Build"
echo "========================"
echo "Build context: $PWD"
echo "Log file: $BUILD_LOG"
echo "Cache strategy: $(if [[ "$USE_CACHE" == "true" ]]; then echo "Enabled"; else echo "Disabled"; fi)"
echo "Force rebuild: $(if [[ "$FORCE_REBUILD" == "true" ]]; then echo "Yes"; else echo "No"; fi)"
echo ""

# Prepare cache arguments
CACHE_ARGS=()
if [[ "$FORCE_REBUILD" == "true" ]]; then
    CACHE_ARGS+=("--no-cache")
    echo "🔄 Force rebuild: Ignoring all cache layers"
elif [[ "$USE_CACHE" == "true" ]]; then
    echo "📦 Preparing cache strategy..."
    
    # Check which cache images are available
    AVAILABLE_CACHE=()
    for cache_image in "${CACHE_FROM_IMAGES[@]}"; do
        if docker image inspect "$cache_image" >/dev/null 2>&1; then
            AVAILABLE_CACHE+=("$cache_image")
            echo "   ✅ Cache available: $cache_image"
        else
            echo "   ⚠️  Cache not found: $cache_image"
        fi
    done
    
    # Add cache-from arguments for available images
    for cache_image in "${AVAILABLE_CACHE[@]}"; do
        CACHE_ARGS+=("--cache-from" "$cache_image")
    done
    
    if [[ ${#AVAILABLE_CACHE[@]} -gt 0 ]]; then
        echo "   📊 Using ${#AVAILABLE_CACHE[@]} cache source(s)"
    else
        echo "   ⚠️  No cache images available - building from scratch"
    fi
else
    CACHE_ARGS+=("--no-cache")
    echo "🚫 Cache disabled by user"
fi

# Start Docker build
echo "📝 Building container (output logged to $BUILD_LOG)..."

if docker build \
    --memory=12g \
    --memory-swap=16g \
    --build-arg CUSTOM_TMPDIR="/host/cassio/export03/data/enning" \
    "${CACHE_ARGS[@]}" \
    --tag micapipe:latest \
    --tag "$CACHE_TAG" \
    . 2>&1 | tee "$BUILD_LOG"; then
    
    echo ""
    echo "✅ Build Successful!"
    echo "=================="
    echo "Container: micapipe:latest"
    echo "Cache tag: $CACHE_TAG"
    echo "Build log: $BUILD_LOG"
    echo ""
    
    # Save cache image info
    echo "💾 Caching build for future use..."
    echo "   Main image: micapipe:latest"
    echo "   Cache image: $CACHE_TAG"
    
    # Show cache statistics
    if [[ "$USE_CACHE" == "true" ]] && [[ ${#AVAILABLE_CACHE[@]} -gt 0 ]]; then
        echo "   📊 Cache hits: $(grep -c "Using cache" "$BUILD_LOG" || echo "0")"
        echo "   📊 Cache misses: $(grep -c "Running in" "$BUILD_LOG" || echo "0")"
    fi
    
    echo ""
    echo "To test the container:"
    echo "docker run --rm micapipe:latest micapipe --help"
    
else
    BUILD_EXIT_CODE=${PIPESTATUS[0]}
    echo ""
    echo "❌ Docker build failed (exit code: $BUILD_EXIT_CODE)"
    echo "📋 Check build log: $BUILD_LOG"
    echo ""
    echo "Common issues:"
    echo "1. Pre-downloaded files not accessible"
    echo "2. Network connectivity issues"
    echo "3. Insufficient disk space"
    echo "4. Memory constraints"
    exit $BUILD_EXIT_CODE
fi