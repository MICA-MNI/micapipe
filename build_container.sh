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
CACHE_DIR="/data_/mica1/01_programs/micapipe_cache"
DOWNLOADS_DIR=""
CUSTOM_TMPDIR="/host/cassio/export03/data/enning"
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
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --downloads-dir)
            DOWNLOADS_DIR="$2"
            shift 2
            ;;
        --custom-tmpdir)
            CUSTOM_TMPDIR="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--cuda] [--no-cache] [--singularity-path PATH] [--cache-dir PATH] [--downloads-dir PATH] [--custom-tmpdir PATH]"
            echo "  --cuda                Enable CUDA support"
            echo "  --no-cache            Build without Docker cache"
            echo "  --singularity-path    Custom path for .sif output (default: /data_/mica1/01_programs/singularity)"
            echo "  --cache-dir           Path to dependency cache (default: /data_/mica1/01_programs/micapipe_cache)"
            echo "  --downloads-dir       Path to pre-downloaded dependencies (use after running download_dependencies.sh)"
            echo "  --custom-tmpdir       Custom temporary directory for build operations (default: /host/cassio/export03/data/enning)"
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

# Check system memory
echo "💾 Checking system memory..."
if command -v free >/dev/null 2>&1; then
    AVAILABLE_MEM=$(free -m | awk 'NR==2{printf "%.0f", $7/1024}')
    echo "   Available: ${AVAILABLE_MEM}GB"
elif command -v vm_stat >/dev/null 2>&1; then
    # macOS
    FREE_PAGES=$(vm_stat | grep "Pages free" | awk '{print $3}' | sed 's/\.//')
    AVAILABLE_MEM=$((FREE_PAGES * 4096 / 1024 / 1024 / 1024))
    echo "   Available: ${AVAILABLE_MEM}GB (approximate)"
fi

# Set environment
export DOCKER_CONTENT_TRUST=0
export BUILDKIT_PROGRESS=plain

# Ensure cache directory exists if specified and accessible
if [[ "$NO_CACHE" == "false" ]]; then
    if [[ -d "$CACHE_DIR" && -w "$CACHE_DIR" ]]; then
        echo "📦 Using cache directory: $CACHE_DIR"
        # Test if BuildKit --mount is supported
        export DOCKER_BUILDKIT=1
        if docker build --help 2>/dev/null | grep -q "\-\-mount"; then
            DOCKER_CACHE_ARGS="--mount type=bind,source=$CACHE_DIR,target=/cache"
        else
            echo "🔧 BuildKit not available, skipping cache mount"
            DOCKER_CACHE_ARGS=""
        fi
        # Ensure cache script is available
        if [[ -f "./cache_dependencies.sh" ]]; then
            echo "🔧 Running cache setup..."
            chmod +x ./cache_dependencies.sh
            ./cache_dependencies.sh
        fi
    else
        echo "⚠️  Cache directory not accessible: $CACHE_DIR"
        echo "📂 Building without cache mount"
        DOCKER_CACHE_ARGS=""
    fi
else
    DOCKER_CACHE_ARGS=""
fi

# Ensure downloads directory exists if specified and accessible
DOCKER_DOWNLOADS_ARGS=""
COPY_DOWNLOADS=false
if [[ -n "$DOWNLOADS_DIR" ]]; then
    if [[ -d "$DOWNLOADS_DIR" && -r "$DOWNLOADS_DIR" ]]; then
        echo "📦 Using downloads directory: $DOWNLOADS_DIR"
        
        # Test if BuildKit --mount is supported
        export DOCKER_BUILDKIT=1
        if docker build --help 2>/dev/null | grep -q "\-\-mount"; then
            echo "🔧 Using Docker BuildKit for mount support"
            DOCKER_DOWNLOADS_ARGS="--mount type=bind,source=$DOWNLOADS_DIR,target=/downloads"
        else
            # Fall back to copying files into build context
            echo "🔧 BuildKit not available, copying downloads to build context"
            COPY_DOWNLOADS=true
            mkdir -p ./temp_downloads
            if [[ -n "$(ls -A "$DOWNLOADS_DIR" 2>/dev/null)" ]]; then
                cp -r "$DOWNLOADS_DIR"/* ./temp_downloads/ 2>/dev/null || true
                echo "   Copied $(ls -1 ./temp_downloads/ 2>/dev/null | wc -l) files"
            fi
        fi
    else
        echo "⚠️  Downloads directory not accessible: $DOWNLOADS_DIR"
        echo "📂 Building without downloads mount"
    fi
fi

# Create log directory
mkdir -p build_logs
LOG_FILE="build_logs/container_build_$(date +%Y%m%d_%H%M%S).log"

echo "   CUDA enabled: $ENABLE_CUDA"
echo "   No cache: $NO_CACHE"
echo "   Output path: $OUTPUT_PATH"
echo "   Custom tmpdir: $CUSTOM_TMPDIR"
echo "   Log file: $LOG_FILE"

# Test network
if curl -I --max-time 10 https://fsl.fmrib.ox.ac.uk >/dev/null 2>&1; then
    echo "✅ Network connectivity OK"
else
    echo "⚠️  Network issues detected"
fi

echo ""
echo "🚀 Building Docker Image..."

# Build command with memory optimizations
DOCKER_CMD="docker build --progress=plain --network=host"

# Add memory limit if sufficient RAM available
if [[ -n "$AVAILABLE_MEM" ]] && [[ "$AVAILABLE_MEM" -ge 8 ]]; then
    DOCKER_CMD="$DOCKER_CMD --memory=8g"
    echo "   Memory limit: 8GB"
fi

DOCKER_CMD="$DOCKER_CMD --build-arg ENABLE_CUDA=$ENABLE_CUDA"

# Add custom tmpdir build argument
DOCKER_CMD="$DOCKER_CMD --build-arg CUSTOM_TMPDIR=$CUSTOM_TMPDIR"

# Add cache mount if available
if [[ -n "$DOCKER_CACHE_ARGS" ]]; then
    DOCKER_CMD="$DOCKER_CMD $DOCKER_CACHE_ARGS --build-arg MICAPIPE_CACHE_DIR=/cache"
    echo "   Cache: Enabled"
else
    echo "   Cache: Disabled"
fi

# Add downloads mount if available
if [[ -n "$DOCKER_DOWNLOADS_ARGS" ]]; then
    DOCKER_CMD="$DOCKER_CMD $DOCKER_DOWNLOADS_ARGS --build-arg DOWNLOADS_DIR=/downloads"
    echo "   Pre-downloads: Enabled (mounted)"
elif [[ "$COPY_DOWNLOADS" == "true" ]]; then
    DOCKER_CMD="$DOCKER_CMD --build-arg DOWNLOADS_DIR=/temp_downloads"
    echo "   Pre-downloads: Enabled (copied)"
else
    echo "   Pre-downloads: Disabled"
fi

if [[ "$NO_CACHE" == "true" ]]; then
    DOCKER_CMD="$DOCKER_CMD --no-cache"
fi

DOCKER_CMD="$DOCKER_CMD -t micapipe:v1-beta ."

echo "   Command: $DOCKER_CMD"
echo ""

# Execute build
eval "$DOCKER_CMD" 2>&1 | tee "$LOG_FILE"
BUILD_EXIT_CODE=${PIPESTATUS[0]}

# Cleanup temporary downloads directory if it was created
if [[ "$COPY_DOWNLOADS" == "true" && -d "./temp_downloads" ]]; then
    echo "🧹 Cleaning up temporary downloads..."
    rm -rf ./temp_downloads
fi

if [[ $BUILD_EXIT_CODE -eq 137 ]]; then
    echo ""
    echo "❌ Docker build failed (exit code: 137)"
    echo ""
    echo "🔧 Memory constraint detected:"
    echo "1. Close other applications to free RAM"
    echo "2. Increase Docker Desktop memory limit to 8GB+"
    echo "3. Build during off-peak hours"
    echo "4. Try building with explicit memory limit:"
    echo "   docker build --memory=8g -t micapipe:v1-beta ."
    echo ""
    echo "📋 Check build log: $LOG_FILE"
    exit 137
elif [[ $BUILD_EXIT_CODE -ne 0 ]]; then
    echo ""
    echo "❌ Docker build failed (exit code: $BUILD_EXIT_CODE)"
    echo "📋 Check build log: $LOG_FILE"
    exit $BUILD_EXIT_CODE
fi

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