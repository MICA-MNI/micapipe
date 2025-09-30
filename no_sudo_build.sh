#!/bin/bash

# No-Sudo Docker Build Script for MICApipe
# Designed for users without administrative privileges

echo "🔐 MICApipe No-Sudo Build"
echo "========================"

# Check if we're in the right directory
if [ ! -f "Dockerfile" ]; then
    echo "❌ Dockerfile not found. Please run this from the micapipe directory."
    exit 1
fi

echo ""
echo "🔍 Checking Docker permissions..."

# Test Docker access without sudo
if ! docker ps >/dev/null 2>&1; then
    echo "❌ Cannot access Docker daemon. Possible solutions:"
    echo ""
    echo "1. Add your user to the docker group (requires admin):"
    echo "   sudo usermod -aG docker $USER"
    echo "   newgrp docker  # or logout/login"
    echo ""
    echo "2. Ask your administrator to add you to the docker group"
    echo ""
    echo "3. Use Docker Desktop (if available)"
    echo ""
    echo "4. Use Singularity instead (if available):"
    echo "   singularity build micapipe.sif docker://path/to/your/image"
    echo ""
    exit 1
fi

echo "✅ Docker access confirmed"

echo ""
echo "📊 System Resource Check..."

# Check available resources (no sudo needed)
TOTAL_MEM=$(free -g 2>/dev/null | grep "Mem:" | awk '{print $2}' || echo "unknown")
AVAILABLE_MEM=$(free -g 2>/dev/null | grep "Mem:" | awk '{print $7}' || echo "unknown")
DISK_SPACE=$(df -h . | tail -1 | awk '{print $4}')

echo "   Total Memory: ${TOTAL_MEM}GB"
echo "   Available Memory: ${AVAILABLE_MEM}GB"
echo "   Available Disk: $DISK_SPACE"

# Warning for low resources
if [[ "$AVAILABLE_MEM" != "unknown" && "$AVAILABLE_MEM" -lt 8 ]]; then
    echo "⚠️  Warning: Less than 8GB available memory"
    echo "   Large downloads may fail - consider closing other applications"
fi

echo ""
echo "🧹 Cleaning Docker environment (user-level)..."

# Clean user-accessible Docker resources
docker container prune -f >/dev/null 2>&1 || echo "   Container cleanup completed"
docker image prune -f >/dev/null 2>&1 || echo "   Image cleanup completed"
docker builder prune -f >/dev/null 2>&1 || echo "   Builder cache cleanup completed"

echo "✅ User-level Docker cleanup completed"

echo ""
echo "⚙️  Setting up build environment..."

# Set environment variables (no sudo needed)
export DOCKER_CONTENT_TRUST=0
export BUILDKIT_PROGRESS=plain

# Check Docker version and supported features
DOCKER_VERSION=$(docker version --format '{{.Server.Version}}' 2>/dev/null || echo "unknown")
echo "   Docker version: $DOCKER_VERSION"

# Test network connectivity
echo ""
echo "🌐 Testing network connectivity..."
if curl -I --max-time 10 https://fsl.fmrib.ox.ac.uk >/dev/null 2>&1; then
    echo "✅ Network connectivity OK"
else
    echo "⚠️  Network issues detected - downloads may be slow"
fi

echo ""
echo "🚀 Starting Docker build..."

# Create build log directory
mkdir -p build_logs
LOG_FILE="build_logs/no_sudo_build_$(date +%Y%m%d_%H%M%S).log"

# Build with minimal flags for maximum compatibility
echo "🔄 Building with compatibility mode..."
echo "   Note: Using basic Docker build flags for older Docker versions"
echo "   Log file: $LOG_FILE"
echo ""

# Start build with basic options
docker build \
    --progress=plain \
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
    if docker run --rm micapipe:v1-beta echo "Container test successful" >/dev/null 2>&1; then
        echo "✅ Container test passed"
        
        # Check if Singularity is available for conversion
        if command -v singularity >/dev/null 2>&1; then
            echo ""
            echo "🔄 Singularity detected. Converting to .sif format..."
            
            # Use user's home directory for .sif file (no sudo needed)
            SIF_PATH="$HOME/micapipe_v1-beta.sif"
            
            if singularity build "$SIF_PATH" docker-daemon://micapipe:v1-beta; then
                echo "✅ Singularity image created: $SIF_PATH"
                echo ""
                echo "📁 Next steps:"
                echo "1. Copy .sif file to your desired location:"
                echo "   cp $SIF_PATH /path/to/your/destination/"
                echo ""
                echo "2. Test the Singularity container:"
                echo "   singularity exec $SIF_PATH micapipe --help"
                echo ""
                echo "3. Run with your data:"
                echo "   singularity exec --bind /path/to/data:/data $SIF_PATH micapipe [options]"
            else
                echo "❌ Singularity conversion failed"
                echo "   You can still use the Docker image: micapipe:v1-beta"
            fi
        else
            echo ""
            echo "📁 Docker image ready: micapipe:v1-beta"
            echo ""
            echo "To save/share this image:"
            echo "1. Save to tar file: docker save micapipe:v1-beta | gzip > micapipe_v1-beta.tar.gz"
            echo "2. Load on another system: docker load < micapipe_v1-beta.tar.gz"
        fi
        
    else
        echo "❌ Container test failed"
        BUILD_EXIT_CODE=1
    fi
    
elif [[ $BUILD_EXIT_CODE -eq 137 ]]; then
    echo "❌ Build killed due to memory constraints (exit code 137)"
    echo ""
    echo "🔧 Troubleshooting for no-sudo users:"
    echo "1. Close memory-intensive applications"
    echo "2. Ask admin to increase swap space"
    echo "3. Build during off-peak hours"
    echo "4. Use a machine with more RAM"
    echo "5. Contact admin about Docker memory limits"
    
elif [[ $BUILD_EXIT_CODE -eq 125 ]]; then
    echo "❌ Docker command error (exit code 125)"
    echo ""
    echo "🔧 Common solutions:"
    echo "1. Check Docker version compatibility"
    echo "2. Ensure user is in docker group"
    echo "3. Try: docker build -t micapipe:v1-beta ."
    echo "4. Contact admin if permissions issue"
    
else
    echo "❌ Build failed with exit code: $BUILD_EXIT_CODE"
    echo ""
    echo "📋 Check build log for details: $LOG_FILE"
    echo ""
    echo "🔍 Common issues for no-sudo users:"
    echo "1. Docker group membership"
    echo "2. Network connectivity"
    echo "3. Disk space in user directory"
    echo "4. Docker daemon permissions"
fi

echo ""
echo "📁 Build log saved: $LOG_FILE"
echo "📊 For additional help:"
echo "   - Review the log file for specific errors"
echo "   - Check available disk space: df -h"
echo "   - Verify Docker access: docker ps"
echo "   - Contact system administrator if needed"

exit $BUILD_EXIT_CODE