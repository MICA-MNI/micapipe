#!/bin/bash

# Quick Build Retry Script
# Applies common fixes and retries the build

echo "🔄 MICApipe Build Retry with Common Fixes"
echo "========================================="

# Check if we're in the right directory
if [ ! -f "Dockerfile" ]; then
    echo "❌ Dockerfile not found. Please run this from the micapipe directory."
    exit 1
fi

echo ""
echo "🧹 Step 1: Cleaning up Docker environment..."

# Clean up Docker system
echo "   Removing old containers and images..."
docker system prune -f > /dev/null 2>&1

# Remove any failed micapipe containers
docker ps -a | grep micapipe | awk '{print $1}' | xargs -r docker rm -f > /dev/null 2>&1

# Remove any partial micapipe images
docker images | grep micapipe | grep -E "(none|intermediate)" | awk '{print $3}' | xargs -r docker rmi -f > /dev/null 2>&1

echo "✅ Docker cleanup completed"

echo ""
echo "⚙️  Step 2: Setting optimized build environment..."

# Set optimal Docker environment
export DOCKER_CONTENT_TRUST=0
export DOCKER_BUILDKIT=1
export DOCKER_BUILDKIT_TIMEOUT=10800  # 3 hours
export BUILDKIT_PROGRESS=plain

echo "✅ Environment configured"

echo ""
echo "🌐 Step 3: Testing network connectivity..."

# Quick network test
if ping -c 2 fsl.fmrib.ox.ac.uk > /dev/null 2>&1; then
    echo "✅ Network connectivity OK"
else
    echo "⚠️  Network connectivity issues detected"
    echo "   Consider running during off-peak hours"
fi

echo ""
echo "💾 Step 4: Checking system resources..."

# Check disk space
DISK_SPACE=$(df . | awk 'NR==2 {print $4}')
DISK_GB=$((DISK_SPACE / 1024 / 1024))

if [ $DISK_GB -lt 20 ]; then
    echo "⚠️  Low disk space: ${DISK_GB}GB available"
    echo "   Recommend at least 20GB for safe build"
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    echo "✅ Sufficient disk space: ${DISK_GB}GB available"
fi

echo ""
echo "🚀 Step 5: Starting optimized build..."

# Determine CUDA support
CUDA_SUPPORT="false"
if command -v nvidia-smi &> /dev/null; then
    echo "   NVIDIA GPU detected - enabling CUDA support"
    CUDA_SUPPORT="true"
else
    echo "   No NVIDIA GPU detected - CPU-only build"
fi

# Choose build method
echo ""
echo "📋 Choose build method:"
echo "1. Quick Singularity build (./build_container.sh)"
echo "2. Comprehensive build with testing (./server_build_test.sh)"
echo "3. Manual Docker build only"
read -p "Choose option (1-3): " -n 1 -r
echo

case $REPLY in
    1)
        echo "🔄 Running quick Singularity build..."
        if [ "$CUDA_SUPPORT" = "true" ]; then
            ./build_container.sh --cuda
        else
            ./build_container.sh
        fi
        ;;
    2)
        echo "🔄 Running comprehensive build with testing..."
        ./server_build_test.sh
        ;;
    3)
        echo "🔄 Running manual Docker build..."
        docker build \
            --progress=plain \
            --network=host \
            --build-arg ENABLE_CUDA=$CUDA_SUPPORT \
            -t micapipe:latest \
            . 2>&1 | tee retry_build.log
        ;;
    *)
        echo "❌ Invalid option selected"
        exit 1
        ;;
esac

echo ""
echo "🎯 Build retry completed!"
echo "📊 Check the output above for results"
echo "📝 Logs saved in: ./build_logs/ (for options 1-2) or retry_build.log (for option 3)"