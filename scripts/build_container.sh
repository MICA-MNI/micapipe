#!/bin/bash

# Micapipe Docker Container Build Script
# This script builds the micapipe Docker container for server deployment
# Based on the CI workflow from .github/workflows/ci_test.yml

set -e  # Exit on any error

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
DEFAULT_TAG="micapipe:latest"
DEFAULT_SINGULARITY_DIR="/tmp/singularity_build"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Build micapipe Docker container for server deployment

OPTIONS:
    -t, --tag TAG                Docker tag to use (default: $DEFAULT_TAG)
    -s, --singularity            Also build Singularity image
    -d, --singularity-dir DIR    Directory for Singularity build (default: $DEFAULT_SINGULARITY_DIR)
    -n, --no-cache              Build without using Docker cache
    -c, --cuda                  Enable CUDA support in the container (default: disabled)
    -h, --help                  Show this help message

EXAMPLES:
    # Basic build
    $0

    # Build with custom tag
    $0 --tag micapipe:v0.2.4

    # Build both Docker and Singularity images
    $0 --tag micapipe:v0.2.4 --singularity

    # Build without cache (clean build)
    $0 --no-cache

    # Enable CUDA support
    $0 --cuda --tag micapipe:v0.2.4-cuda

    # Build CUDA-enabled container with Singularity
    $0 --cuda --singularity --tag micapipe:v0.2.4-cuda

EOF
}

# Parse command line arguments
DOCKER_TAG="$DEFAULT_TAG"
BUILD_SINGULARITY=false
SINGULARITY_DIR="$DEFAULT_SINGULARITY_DIR"
NO_CACHE=false
ENABLE_CUDA=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--tag)
            DOCKER_TAG="$2"
            shift 2
            ;;
        -s|--singularity)
            BUILD_SINGULARITY=true
            shift
            ;;
        -d|--singularity-dir)
            SINGULARITY_DIR="$2"
            shift 2
            ;;
        -n|--no-cache)
            NO_CACHE=true
            shift
            ;;
        -c|--cuda)
            ENABLE_CUDA=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Validate prerequisites
print_info "Checking prerequisites..."

# Check if Docker is installed and running
if ! command -v docker &> /dev/null; then
    print_error "Docker is not installed or not in PATH"
    exit 1
fi

if ! docker info &> /dev/null; then
    print_error "Docker daemon is not running"
    exit 1
fi

# Check if Singularity is available (if requested)
if [ "$BUILD_SINGULARITY" = true ]; then
    if ! command -v singularity &> /dev/null; then
        print_error "Singularity is not installed or not in PATH"
        exit 1
    fi
    
    # Create Singularity build directory
    mkdir -p "$SINGULARITY_DIR"
    if [ ! -w "$SINGULARITY_DIR" ]; then
        print_error "Cannot write to Singularity directory: $SINGULARITY_DIR"
        exit 1
    fi
fi

# Check NVIDIA/CUDA availability (if requested)
if [ "$ENABLE_CUDA" = true ]; then
    print_info "Checking CUDA prerequisites..."
    
    # Check if nvidia-smi is available (indicates NVIDIA drivers)
    if command -v nvidia-smi &> /dev/null; then
        print_info "NVIDIA drivers detected: $(nvidia-smi --query-gpu=name --format=csv,noheader,nounits | head -1)"
    else
        print_warning "nvidia-smi not found. CUDA container may not work on this host."
        print_warning "Ensure NVIDIA drivers are installed if you plan to run CUDA containers."
    fi
    
    # Check if nvidia-docker or Docker with GPU support is available
    if docker info 2>/dev/null | grep -i nvidia >/dev/null; then
        print_info "Docker NVIDIA runtime detected"
    else
        print_warning "Docker NVIDIA runtime not detected. You may need to install nvidia-docker2."
        print_warning "CUDA containers may not work without proper GPU runtime."
    fi
fi

# Check if we're in the right directory
if [ ! -f "$PROJECT_ROOT/Dockerfile" ]; then
    print_error "Dockerfile not found. Please run this script from the micapipe project directory."
    exit 1
fi

print_success "Prerequisites check passed"

# Build Docker image
print_info "Building Docker image with tag: $DOCKER_TAG"
print_info "Build directory: $PROJECT_ROOT"

# Prepare Docker build command
DOCKER_BUILD_CMD="docker build"
if [ "$NO_CACHE" = true ]; then
    DOCKER_BUILD_CMD="$DOCKER_BUILD_CMD --no-cache"
    print_info "Building without cache (clean build)"
fi

# Add CUDA build argument
if [ "$ENABLE_CUDA" = true ]; then
    DOCKER_BUILD_CMD="$DOCKER_BUILD_CMD --build-arg ENABLE_CUDA=true"
    print_info "Building with CUDA support enabled"
else
    DOCKER_BUILD_CMD="$DOCKER_BUILD_CMD --build-arg ENABLE_CUDA=false"
    print_info "Building with CPU-only support"
fi

DOCKER_BUILD_CMD="$DOCKER_BUILD_CMD -t $DOCKER_TAG ."

print_info "Running: $DOCKER_BUILD_CMD"

# Change to project root and build
cd "$PROJECT_ROOT"

# Capture build time
BUILD_START=$(date +%s)

if eval "$DOCKER_BUILD_CMD"; then
    BUILD_END=$(date +%s)
    BUILD_TIME=$((BUILD_END - BUILD_START))
    print_success "Docker image built successfully in ${BUILD_TIME} seconds"
    print_success "Image tag: $DOCKER_TAG"
else
    print_error "Docker build failed"
    exit 1
fi

# Build Singularity image if requested
if [ "$BUILD_SINGULARITY" = true ]; then
    # Extract tag name for Singularity file
    SINGULARITY_TAG=$(echo "$DOCKER_TAG" | sed 's/:/_/g')
    SINGULARITY_FILE="$SINGULARITY_DIR/micapipe_${SINGULARITY_TAG}_${TIMESTAMP}.sif"
    
    print_info "Building Singularity image: $SINGULARITY_FILE"
    
    # Set SINGULARITY_TMPDIR if not already set
    export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-$SINGULARITY_DIR/tmp}"
    mkdir -p "$SINGULARITY_TMPDIR"
    
    print_info "Using SINGULARITY_TMPDIR: $SINGULARITY_TMPDIR"
    
    SINGULARITY_BUILD_START=$(date +%s)
    
    if singularity build --force "$SINGULARITY_FILE" docker-daemon://"$DOCKER_TAG"; then
        SINGULARITY_BUILD_END=$(date +%s)
        SINGULARITY_BUILD_TIME=$((SINGULARITY_BUILD_END - SINGULARITY_BUILD_START))
        print_success "Singularity image built successfully in ${SINGULARITY_BUILD_TIME} seconds"
        print_success "Singularity file: $SINGULARITY_FILE"
        
        # Show file size
        SINGULARITY_SIZE=$(du -h "$SINGULARITY_FILE" | cut -f1)
        print_info "Singularity image size: $SINGULARITY_SIZE"
    else
        print_error "Singularity build failed"
        exit 1
    fi
fi

# Show final summary
print_success "Build completed successfully!"
echo ""
print_info "Summary:"
echo "  Docker image: $DOCKER_TAG"
echo "  CUDA support: $([ "$ENABLE_CUDA" = true ] && echo "Enabled" || echo "Disabled")"
if [ "$BUILD_SINGULARITY" = true ]; then
    echo "  Singularity file: $SINGULARITY_FILE"
fi
echo "  Build timestamp: $TIMESTAMP"

# Show image information
print_info "Docker image details:"
docker images "$DOCKER_TAG" --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

# Provide usage examples
echo ""
print_info "Usage examples:"
if [ "$ENABLE_CUDA" = true ]; then
    echo "  # Run Docker container with GPU support:"
    echo "  docker run --gpus all -it --rm $DOCKER_TAG"
else
    echo "  # Run Docker container interactively:"
    echo "  docker run -it --rm $DOCKER_TAG"
fi
echo ""
if [ "$BUILD_SINGULARITY" = true ]; then
    if [ "$ENABLE_CUDA" = true ]; then
        echo "  # Run Singularity container with GPU support:"
        echo "  singularity exec --nv $SINGULARITY_FILE micapipe --help"
    else
        echo "  # Run Singularity container:"
        echo "  singularity exec $SINGULARITY_FILE micapipe --help"
    fi
    echo ""
fi

print_success "All done!"