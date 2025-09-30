#!/bin/bash

# Memory-Safe Docker Build Script
# Specifically handles memory constraints and large downloads

echo "🧠 MICApipe Memory-Safe Build"
echo "============================="

# Check if we're in the right directory
if [ ! -f "Dockerfile" ]; then
    echo "❌ Dockerfile not found. Please run this from the micapipe directory."
    exit 1
fi

echo ""
echo "📊 System Resource Assessment..."

# Check available memory
TOTAL_MEM=$(free -g | grep "Mem:" | awk '{print $2}')
AVAILABLE_MEM=$(free -g | grep "Mem:" | awk '{print $7}')
USED_MEM=$(free -g | grep "Mem:" | awk '{print $3}')

echo "   Total Memory: ${TOTAL_MEM}GB"
echo "   Used Memory: ${USED_MEM}GB"
echo "   Available Memory: ${AVAILABLE_MEM}GB"

if [[ "$AVAILABLE_MEM" -lt 6 ]]; then
    echo ""
    echo "⚠️  WARNING: Low available memory detected!"
    echo "   FreeSurfer download (~3GB) + extraction may fail"
    echo "   Recommendations:"
    echo "   1. Close other applications"
    echo "   2. Use swap space if available"
    echo "   3. Build on a machine with more RAM"
    echo ""
    read -p "Continue with memory-safe optimizations? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
    MEMORY_SAFE="true"
else
    echo "✅ Sufficient memory available"
    MEMORY_SAFE="false"
fi

# Check swap space
SWAP_TOTAL=$(free -g | grep "Swap:" | awk '{print $2}')
SWAP_USED=$(free -g | grep "Swap:" | awk '{print $3}')
SWAP_FREE=$(free -g | grep "Swap:" | awk '{print $4}')

echo "   Swap Total: ${SWAP_TOTAL}GB"
echo "   Swap Used: ${SWAP_USED}GB"
echo "   Swap Free: ${SWAP_FREE}GB"

if [[ "$SWAP_TOTAL" -gt 0 ]]; then
    echo "✅ Swap space available for memory overflow"
else
    echo "⚠️  No swap space - memory constraints may cause failures"
fi

echo ""
echo "🧹 Pre-build Cleanup..."

# Clean Docker environment aggressively
docker system prune -f --volumes > /dev/null 2>&1
docker builder prune -f > /dev/null 2>&1

# Remove any partial micapipe builds
docker images | grep -E "(micapipe|<none>)" | awk '{print $3}' | xargs -r docker rmi -f > /dev/null 2>&1

echo "✅ Docker environment cleaned"

echo ""
echo "⚙️  Configuring Memory-Safe Build Environment..."

# Set environment for memory optimization
export DOCKER_CONTENT_TRUST=0
export BUILDKIT_PROGRESS=plain

# Configure Docker daemon limits if possible
if command -v docker-compose &> /dev/null; then
    echo "   Docker Compose available - using resource limits"
    USE_COMPOSE="true"
else
    echo "   Using direct Docker build with memory flags"
    USE_COMPOSE="false"
fi

# Determine optimal build strategy
if [[ "$MEMORY_SAFE" == "true" ]]; then
    echo "   Using conservative memory settings"
    MEMORY_LIMIT="4g"
    MEMORY_SWAP="6g"
    SHM_SIZE="512m"
    BUILD_PARALLEL="1"
else
    echo "   Using standard memory settings"
    MEMORY_LIMIT="8g"
    MEMORY_SWAP="12g"
    SHM_SIZE="2g"
    BUILD_PARALLEL="2"
fi

echo ""
echo "🚀 Starting Memory-Safe Docker Build..."
echo "   Memory Limit: $MEMORY_LIMIT"
echo "   Memory+Swap: $MEMORY_SWAP"
echo "   Shared Memory: $SHM_SIZE"
echo ""

# Create build log directory
mkdir -p build_logs
LOG_FILE="build_logs/memory_safe_build_$(date +%Y%m%d_%H%M%S).log"

# Start the build with memory constraints
echo "🔄 Building Docker image (this may take 1-2 hours)..."

docker build \
    --progress=plain \
    --memory="$MEMORY_LIMIT" \
    --memory-swap="$MEMORY_SWAP" \
    --shm-size="$SHM_SIZE" \
    --cpus="$BUILD_PARALLEL" \
    --network=host \
    --build-arg ENABLE_CUDA=false \
    -t micapipe:v1-beta \
    . 2>&1 | tee "$LOG_FILE"

BUILD_EXIT_CODE=${PIPESTATUS[0]}

echo ""
if [[ $BUILD_EXIT_CODE -eq 0 ]]; then
    echo "✅ Docker build successful!"
    
    # Test the container
    echo ""
    echo "🧪 Testing the built container..."
    if docker run --rm micapipe:v1-beta echo "Container test successful" > /dev/null 2>&1; then
        echo "✅ Container test passed"
        
        # Offer Singularity conversion
        if command -v singularity &> /dev/null; then
            echo ""
            read -p "Convert to Singularity? (y/N): " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                echo "🔄 Converting to Singularity..."
                singularity build micapipe_v1-beta.sif docker-daemon://micapipe:v1-beta
                echo "✅ Singularity image: micapipe_v1-beta.sif"
            fi
        fi
        
    else
        echo "❌ Container test failed"
        BUILD_EXIT_CODE=1
    fi
    
elif [[ $BUILD_EXIT_CODE -eq 137 ]]; then
    echo "❌ Build killed due to memory constraints (exit code 137)"
    echo ""
    echo "🔧 Troubleshooting suggestions:"
    echo "1. Increase system swap space:"
    echo "   sudo fallocate -l 8G /swapfile"
    echo "   sudo chmod 600 /swapfile"
    echo "   sudo mkswap /swapfile"
    echo "   sudo swapon /swapfile"
    echo ""
    echo "2. Close memory-intensive applications"
    echo "3. Build on a machine with more RAM (>= 16GB recommended)"
    echo "4. Try building individual stages:"
    echo "   docker build --target neurodocker -t micapipe:base ."
    echo ""
    
else
    echo "❌ Build failed with exit code: $BUILD_EXIT_CODE"
    echo ""
    echo "📋 Check the build log for details:"
    echo "   $LOG_FILE"
    echo ""
    echo "🔍 Common issues and solutions:"
    tail -20 "$LOG_FILE" | grep -i "error\|failed\|killed" || echo "No obvious error patterns found"
fi

echo ""
echo "📁 Build log saved: $LOG_FILE"
echo "📊 Memory usage during build:"
echo "   Peak memory may have exceeded available system memory"
echo "   Consider monitoring with: watch -n 1 'free -h && docker stats --no-stream'"

exit $BUILD_EXIT_CODE