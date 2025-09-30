#!/bin/bash

# MICApipe Container Builder
# Builds Docker image and converts to Singularity .sif format
# No sudo privileges required

echo "🐳 MICApipe Container Builder"
echo "============================"

# Check if we're in the right directory
if [ ! -f "Dockerfile" ]; then
    echo "❌ Dockerfile not found. Please run this from the micapipe directory."
    exit 1
fi

# Parse command line arguments
ENABLE_CUDA=false
SINGULARITY_PATH="/data_/mica1/01_programs/singularity"
FALLBACK_PATH="$HOME"
NO_CACHE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --cuda)
            ENABLE_CUDA=true
            echo "🎯 CUDA support enabled"
            shift
            ;;
        --no-cache)
            NO_CACHE=true
            echo "🔄 Building without cache"
            shift
            ;;
        --singularity-path)
            SINGULARITY_PATH="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--cuda] [--no-cache] [--singularity-path PATH]"
            echo "  --cuda                Enable CUDA support"
            echo "  --no-cache            Build without Docker cache"
            echo "  --singularity-path    Custom path for .sif output (default: /data_/mica1/01_programs/singularity)"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo ""
echo "🔍 System Checks..."

# Check Docker access (no sudo required)
if ! docker ps >/dev/null 2>&1; then
    echo "❌ Cannot access Docker daemon. Solutions:"
    echo "1. Ask admin to add you to docker group: sudo usermod -aG docker $USER"
    echo "2. Use Docker Desktop if available"
    echo "3. Contact system administrator"
    exit 1
fi
echo "✅ Docker access confirmed"

# Check Singularity
if ! command -v singularity >/dev/null 2>&1; then
    echo "❌ Singularity not found. Please install Singularity to create .sif files"
    exit 1
fi
echo "✅ Singularity available"

# Check resources
TOTAL_MEM=$(free -g 2>/dev/null | grep "Mem:" | awk '{print $2}' || echo "unknown")
AVAILABLE_MEM=$(free -g 2>/dev/null | grep "Mem:" | awk '{print $7}' || echo "unknown")
DISK_SPACE=$(df -h . | tail -1 | awk '{print $4}')

echo "   Memory: ${AVAILABLE_MEM}GB available / ${TOTAL_MEM}GB total"
echo "   Disk Space: $DISK_SPACE available"

if [[ "$AVAILABLE_MEM" != "unknown" && "$AVAILABLE_MEM" -lt 8 ]]; then
    echo "⚠️  Warning: Less than 8GB available memory - large downloads may fail"
fi

# Determine output path
if [[ -d "$SINGULARITY_PATH" && -w "$SINGULARITY_PATH" ]]; then
    OUTPUT_PATH="$SINGULARITY_PATH"
    echo "✅ Using Singularity path: $OUTPUT_PATH"
elif [[ -w "$FALLBACK_PATH" ]]; then
    OUTPUT_PATH="$FALLBACK_PATH"
    echo "⚠️  Fallback to home directory: $OUTPUT_PATH"
else
    echo "❌ No writable path found for .sif output"
    exit 1
fi

echo ""
echo "🧹 Cleanup..."

# Clean Docker environment (user-level only)
docker container prune -f >/dev/null 2>&1
docker image prune -f >/dev/null 2>&1
docker builder prune -f >/dev/null 2>&1

# Remove any existing micapipe images
docker images | grep micapipe | awk '{print $3}' | xargs -r docker rmi -f >/dev/null 2>&1

echo "✅ Docker cleanup completed"

echo ""
echo "⚙️  Build Configuration..."

# Set environment
export DOCKER_CONTENT_TRUST=0
export BUILDKIT_PROGRESS=plain

# Create log directory
mkdir -p build_logs
LOG_FILE="build_logs/container_build_$(date +%Y%m%d_%H%M%S).log"

echo "   CUDA enabled: $ENABLE_CUDA"
echo "   No cache: $NO_CACHE"
echo "   Output path: $OUTPUT_PATH"
echo "   Log file: $LOG_FILE"

# Test network
if curl -I --max-time 10 https://fsl.fmrib.ox.ac.uk >/dev/null 2>&1; then
    echo "✅ Network connectivity OK"
else
    echo "⚠️  Network issues detected"
fi

echo ""
echo "🚀 Building Docker Image..."

# Build command
DOCKER_CMD="docker build --progress=plain --network=host --build-arg ENABLE_CUDA=$ENABLE_CUDA"

if [[ "$NO_CACHE" == "true" ]]; then
    DOCKER_CMD="$DOCKER_CMD --no-cache"
fi

DOCKER_CMD="$DOCKER_CMD -t micapipe:v1-beta ."

echo "   Command: $DOCKER_CMD"
echo ""

# Execute build
eval "$DOCKER_CMD" 2>&1 | tee "$LOG_FILE"
BUILD_EXIT_CODE=${PIPESTATUS[0]}

if [[ $BUILD_EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "✅ Docker build successful!"
    
    # Test container
    echo ""
    echo "🧪 Testing container..."
    if docker run --rm micapipe:v1-beta echo "Container test successful" >/dev/null 2>&1; then
        echo "✅ Container test passed"
        
        # Convert to Singularity
        echo ""
        echo "🔄 Converting to Singularity..."
        
        SIF_FILE="$OUTPUT_PATH/micapipe_v1-beta.sif"
        
        # Remove existing .sif if present
        [[ -f "$SIF_FILE" ]] && rm -f "$SIF_FILE"
        
        if singularity build "$SIF_FILE" docker-daemon://micapipe:v1-beta; then
            echo "✅ Singularity image created: $SIF_FILE"
            
            # Get file size
            SIF_SIZE=$(du -h "$SIF_FILE" | cut -f1)
            echo "📊 Container size: $SIF_SIZE"
            
            # Cleanup Docker image to save space
            docker rmi micapipe:v1-beta >/dev/null 2>&1
            echo "🧹 Docker image cleaned up"
            
            echo ""
            echo "🎯 Build Complete!"
            echo "📁 Singularity image: $SIF_FILE"
            echo ""
            echo "🚀 Next Steps:"
            echo "1. Test the container:"
            echo "   singularity exec $SIF_FILE micapipe --help"
            echo ""
            echo "2. Run with your data:"
            echo "   singularity exec --bind /path/to/data:/data $SIF_FILE micapipe [options]"
            echo ""
            echo "3. Copy to other systems:"
            echo "   scp $SIF_FILE user@server:/path/to/destination/"
            
        else
            echo "❌ Singularity conversion failed"
            echo "Docker image available: micapipe:v1-beta"
            exit 1
        fi
        
    else
        echo "❌ Container test failed"
        exit 1
    fi
    
else
    echo ""
    echo "❌ Docker build failed (exit code: $BUILD_EXIT_CODE)"
    
    if [[ $BUILD_EXIT_CODE -eq 137 ]]; then
        echo ""
        echo "🔧 Memory constraint detected:"
        echo "1. Close other applications"
        echo "2. Ask admin about increasing swap space"
        echo "3. Build during off-peak hours"
        
    elif [[ $BUILD_EXIT_CODE -eq 125 ]]; then
        echo ""
        echo "🔧 Docker command error:"
        echo "1. Check Docker permissions"
        echo "2. Try: docker build -t micapipe:v1-beta ."
    fi
    
    echo ""
    echo "📋 Check build log: $LOG_FILE"
    exit $BUILD_EXIT_CODE
fi

echo ""
echo "📁 Build log saved: $LOG_FILE"