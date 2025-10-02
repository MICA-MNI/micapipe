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

# Check if we should build in custom temp directory with downloads
BUILD_IN_TEMP_DIR=false
if [[ -n "$DOWNLOADS_DIR" && "$DOWNLOADS_DIR" == *"$CUSTOM_TMPDIR"* ]]; then
    BUILD_IN_TEMP_DIR=true
    echo "📦 Downloads are in custom temp directory: $DOWNLOADS_DIR"
    echo "🔧 Will copy build environment to same location for efficiency"
    
    # Create build directory in the same location as downloads
    BUILD_DIR="${CUSTOM_TMPDIR}/micapipe_build_$(date +%Y%m%d_%H%M%S)"
    echo "📁 Build directory: $BUILD_DIR"
    
    # Copy entire source tree to build location
    echo "📋 Copying source code to build directory..."
    mkdir -p "$BUILD_DIR"
    rsync -av --exclude='.git' --exclude='build_logs' --exclude='*.sif' \
          --exclude='temp_downloads' --exclude='downloads_temp' \
          ./ "$BUILD_DIR/"
    
    # Create local downloads directory in build context
    mkdir -p "$BUILD_DIR/downloads"
    ln -sf "$DOWNLOADS_DIR"/* "$BUILD_DIR/downloads/" 2>/dev/null || true
    
    echo "   Source copied to: $BUILD_DIR"
    echo "   Downloads linked: $(ls -1 "$BUILD_DIR/downloads/" 2>/dev/null | wc -l) files"
    
    DOCKER_DOWNLOADS_ARGS="--build-arg DOWNLOADS_DIR=/downloads"
    
elif [[ -n "$DOWNLOADS_DIR" && -d "$DOWNLOADS_DIR" && -r "$DOWNLOADS_DIR" ]]; then
    echo "📦 Using downloads directory: $DOWNLOADS_DIR"
    echo "📦 Files available: $(ls -1 "$DOWNLOADS_DIR" | wc -l)"
    echo "📦 FSL file: $(ls -lh "$DOWNLOADS_DIR"/fsl-*.tar.gz 2>/dev/null || echo "Not found")"
    echo "📦 FreeSurfer file: $(ls -lh "$DOWNLOADS_DIR"/freesurfer-*.tar.gz 2>/dev/null || echo "Not found")"
    
    # Use the custom temp directory for copying files (has space)
    COPY_TARGET="${CUSTOM_TMPDIR}/docker_build_downloads"
    echo "🔧 Copying files to custom temp directory: $COPY_TARGET"
    
    # Create the copy directory in your custom temp space
    mkdir -p "$COPY_TARGET"
    
    # Copy specific files to the custom temp directory
    if [[ -f "$DOWNLOADS_DIR/fsl-6.0.2-centos6_64.tar.gz" ]]; then
        echo "   Copying FSL (4GB) to $COPY_TARGET..."
        cp "$DOWNLOADS_DIR/fsl-6.0.2-centos6_64.tar.gz" "$COPY_TARGET/"
    fi
    
    if [[ -f "$DOWNLOADS_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]]; then
        echo "   Copying FreeSurfer (9GB) to $COPY_TARGET..."
        cp "$DOWNLOADS_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" "$COPY_TARGET/"
    fi
    
    echo "   Copied files: $(ls -1 "$COPY_TARGET" 2>/dev/null | wc -l)"
    echo "   Total size: $(du -sh "$COPY_TARGET" | cut -f1)"
    
    # Create symbolic links in build context pointing to copied files
    mkdir -p ./downloads_temp
    for file in "$COPY_TARGET"/*; do
        if [[ -f "$file" ]]; then
            ln -sf "$file" "./downloads_temp/$(basename "$file")"
        fi
    done
    
    # Remove downloads_temp from .dockerignore temporarily
    if [[ -f ".dockerignore" ]]; then
        cp .dockerignore .dockerignore.backup
        sed -i '/downloads_temp/d' .dockerignore
    fi
    
    DOCKER_DOWNLOADS_ARGS="--build-arg DOWNLOADS_DIR=/downloads_temp"
    echo "🔧 Using symbolic links to copied files in custom temp space"
else
    echo "⚠️  No downloads directory specified or not accessible"
    echo "📂 Building with internet downloads"
    DOCKER_DOWNLOADS_ARGS=""
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

# Add downloads if available
if [[ -n "$DOCKER_DOWNLOADS_ARGS" ]]; then
    DOCKER_CMD="$DOCKER_CMD $DOCKER_DOWNLOADS_ARGS"
    echo "   Pre-downloads: Enabled (copied to custom temp space)"
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
if [[ "$BUILD_IN_TEMP_DIR" == "true" ]]; then
    echo "🚀 Building in custom temp directory: $BUILD_DIR"
    cd "$BUILD_DIR"
    eval "$DOCKER_CMD" 2>&1 | tee "$LOG_FILE"
    BUILD_EXIT_CODE=${PIPESTATUS[0]}
    cd - > /dev/null
else
    eval "$DOCKER_CMD" 2>&1 | tee "$LOG_FILE"
    BUILD_EXIT_CODE=${PIPESTATUS[0]}
fi

# Cleanup temporary files
if [[ -d "./downloads_temp" ]]; then
    echo "🧹 Cleaning up build context symbolic links..."
    rm -rf ./downloads_temp
    
    # Restore .dockerignore if we backed it up
    if [[ -f ".dockerignore.backup" ]]; then
        mv .dockerignore.backup .dockerignore
        echo "   Restored .dockerignore"
    fi
fi

# Cleanup copied files in custom temp directory (after successful build)
if [[ $BUILD_EXIT_CODE -eq 0 && -n "$COPY_TARGET" && -d "$COPY_TARGET" ]]; then
    echo "🧹 Cleaning up copied files in custom temp space..."
    rm -rf "$COPY_TARGET"
    echo "   Removed: $COPY_TARGET"
fi

# Cleanup build directory if we built in temp directory
if [[ "$BUILD_IN_TEMP_DIR" == "true" && -n "$BUILD_DIR" ]]; then
    if [[ $BUILD_EXIT_CODE -eq 0 ]]; then
        echo "🧹 Cleaning up temporary build directory..."
        rm -rf "$BUILD_DIR"
        echo "   Removed: $BUILD_DIR"
    else
        echo "⚠️  Keeping build directory for debugging: $BUILD_DIR"
    fi
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