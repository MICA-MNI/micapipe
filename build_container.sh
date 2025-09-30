#!/bin/bash

# MICApipe Singularity Container Build Script
# Builds Singularity .sif file with Docker as intermediate step

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Configuration
SIF_VERSION="v1-beta"
BUILD_DATE=$(date +%Y%m%d_%H%M%S)
DOCKER_TAG="micapipe:intermediate-${BUILD_DATE}"
DEFAULT_SIF_PATH="/opt/micapipe"
SIF_NAME="micapipe_${SIF_VERSION}.sif"

# Parse command line options
ENABLE_CUDA=false
NO_CACHE=false
SIF_DIR=""

print_info "🚀 MICApipe Singularity Container Build"
print_info "Docker container is intermediate step - final output is .sif file"
echo

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cuda)
            ENABLE_CUDA=true
            SIF_NAME="micapipe_${SIF_VERSION}_cuda.sif"
            shift
            ;;
        --no-cache)
            NO_CACHE=true
            shift
            ;;
        --output|-o)
            SIF_DIR="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --cuda          Enable CUDA support (default: CPU-only)"
            echo "  --no-cache      Build without Docker cache"
            echo "  --output DIR    Custom output directory for .sif file"
            echo "  --help          Show this help"
            echo ""
            echo "Environment Variables:"
            echo "  MICAPIPE_SIF_PATH   Custom path for .sif file (default: /opt/micapipe)"
            echo ""
            echo "Output:"
            echo "  Creates: \${OUTPUT_DIR}/micapipe_${SIF_VERSION}.sif"
            echo "  Docker intermediate container is automatically removed"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Determine output directory
if [ -n "$SIF_DIR" ]; then
    OUTPUT_DIR="$SIF_DIR"
elif [ -n "$MICAPIPE_SIF_PATH" ]; then
    OUTPUT_DIR="$MICAPIPE_SIF_PATH"
else
    # Try default server paths, fallback to user directory
    if [ -d "$DEFAULT_SIF_PATH" ] && [ -w "$DEFAULT_SIF_PATH" ]; then
        OUTPUT_DIR="$DEFAULT_SIF_PATH"
    elif [ -d "/data_/mica1/01_programs/singularity" ] && [ -w "/data_/mica1/01_programs/singularity" ]; then
        OUTPUT_DIR="/data_/mica1/01_programs/singularity"
    else
        OUTPUT_DIR="$HOME/containers"
        print_warning "Using fallback directory: $OUTPUT_DIR"
    fi
fi

SIF_PATH="$OUTPUT_DIR/$SIF_NAME"

print_info "Configuration:"
print_info "  CUDA support: $([ "$ENABLE_CUDA" = true ] && echo "ENABLED" || echo "DISABLED")"
print_info "  Docker cache: $([ "$NO_CACHE" = true ] && echo "DISABLED" || echo "ENABLED")"
print_info "  Output file: $SIF_PATH"
echo

# Step 1: Check prerequisites
print_info "📋 Step 1: Checking prerequisites..."

# Check if we're in the right directory
if [ ! -f "Dockerfile" ]; then
    print_error "Dockerfile not found. Please run this script from the micapipe project directory."
    exit 1
fi

# Check Docker
if ! command -v docker &> /dev/null; then
    print_error "Docker is not installed or not in PATH"
    exit 1
fi

if ! docker info &> /dev/null; then
    print_error "Docker daemon is not running"
    exit 1
fi

# Check Singularity
if ! command -v singularity &> /dev/null; then
    print_error "Singularity is not installed or not in PATH"
    print_error "Install Singularity: https://docs.sylabs.io/guides/latest/user-guide/"
    exit 1
fi

print_success "Docker and Singularity are available"

# Step 2: Setup output directory
print_info "📁 Step 2: Setting up output directory..."

# Create output directory
if ! mkdir -p "$OUTPUT_DIR" 2>/dev/null; then
    print_error "Cannot create output directory: $OUTPUT_DIR"
    print_error "Try: $0 --output ~/my_containers"
    exit 1
fi

# Test write permissions
if ! touch "$OUTPUT_DIR/.test" 2>/dev/null; then
    print_error "No write permission for output directory: $OUTPUT_DIR"
    exit 1
else
    rm -f "$OUTPUT_DIR/.test"
fi

print_success "Output directory ready: $OUTPUT_DIR"

# Step 3: Check disk space
print_info "💾 Step 3: Checking disk space..."

AVAILABLE_SPACE=$(df "$OUTPUT_DIR" | awk 'NR==2 {print $4}')
SPACE_GB=$((AVAILABLE_SPACE / 1024 / 1024))

print_info "Available space: ${SPACE_GB}GB"

if [ $SPACE_GB -lt 15 ]; then
    print_warning "Low disk space (${SPACE_GB}GB). Recommend at least 15GB for container build."
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Try using a different directory: $0 --output /path/with/more/space"
        exit 1
    fi
fi

# Step 4: Set Docker environment
print_info "🔧 Step 4: Setting up Docker environment..."

# Disable Docker Content Trust to avoid certificate issues
export DOCKER_CONTENT_TRUST=0
export DOCKER_BUILDKIT=1

print_success "Docker environment configured"

# Step 5: Build Docker container (intermediate)
print_info "🏗️  Step 5: Building Docker container (intermediate)..."

# Prepare build arguments
BUILD_ARGS=""
if [ "$ENABLE_CUDA" = true ]; then
    BUILD_ARGS="--build-arg ENABLE_CUDA=true"
    print_info "Building with CUDA support enabled"
else
    BUILD_ARGS="--build-arg ENABLE_CUDA=false"
    print_info "Building CPU-only version"
fi

if [ "$NO_CACHE" = true ]; then
    BUILD_ARGS="$BUILD_ARGS --no-cache"
    print_info "Building without cache"
fi

print_info "Docker build starting... (this may take 30-60 minutes)"
print_info "Building intermediate container: $DOCKER_TAG"

# Run Docker build
if docker build $BUILD_ARGS -t "$DOCKER_TAG" . 2>&1 | tee docker_build.log; then
    print_success "Docker container built successfully"
else
    print_error "Docker build failed. Check docker_build.log for details"
    exit 1
fi

# Step 6: Convert to Singularity
print_info "🔄 Step 6: Converting Docker to Singularity..."

print_info "Creating Singularity container: $SIF_PATH"
print_info "This may take 10-20 minutes..."

# Convert Docker to Singularity
if singularity build "$SIF_PATH" docker-daemon://${DOCKER_TAG} 2>&1 | tee singularity_build.log; then
    print_success "Singularity container created successfully"
else
    print_error "Singularity build failed. Check singularity_build.log for details"
    exit 1
fi

# Step 7: Test Singularity container
print_info "🧪 Step 7: Testing Singularity container..."

if singularity exec "$SIF_PATH" echo "Container test successful" > /dev/null 2>&1; then
    print_success "Singularity container test passed"
else
    print_warning "Singularity container test failed (container may still work)"
fi

# Step 8: Cleanup Docker container
print_info "🧹 Step 8: Cleaning up intermediate Docker container..."

if docker rmi "$DOCKER_TAG" > /dev/null 2>&1; then
    print_success "Intermediate Docker container removed"
else
    print_warning "Could not remove intermediate Docker container: $DOCKER_TAG"
    print_info "You can manually remove it with: docker rmi $DOCKER_TAG"
fi

# Step 9: Final summary
echo
print_success "🎉 Build completed successfully!"
echo
print_info "📊 Summary:"
print_info "  Singularity file: $SIF_PATH"
print_info "  File size: $(du -sh "$SIF_PATH" | cut -f1)"
print_info "  CUDA support: $([ "$ENABLE_CUDA" = true ] && echo "Enabled" || echo "Disabled")"
echo
print_info "🚀 Usage:"
print_info "  Interactive: singularity shell $SIF_PATH"
print_info "  Execute: singularity exec $SIF_PATH [command]"
print_info "  With data: singularity exec --bind /data:/data $SIF_PATH [command]"
echo
print_info "📝 Build logs saved:"
print_info "  Docker build: ./docker_build.log"
print_info "  Singularity build: ./singularity_build.log"
echo
print_success "Ready for use! 🎯"