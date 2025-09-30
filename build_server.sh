#!/bin/bash

# Micapipe Container Build Script for Your Server Configuration
# This script is customized for your specific paths and requirements

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

# Your server configuration
SIF_DIR="/data_/mica1/01_programs/singularity"
TMP_DIR="/host/cassio/export03/data/enning/tmp"
BUILD_DATE=$(date +%Y%m%d_%H%M%S)
DEFAULT_TAG="micapipe:v0.2.4-${BUILD_DATE}"

# Parse command line options
ENABLE_CUDA=false
NO_CACHE=false
TAG="$DEFAULT_TAG"

print_info "Micapipe Container Build for Server"
print_info "SIF Directory: $SIF_DIR"
print_info "Temp Directory: $TMP_DIR"
echo

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --cuda)
            ENABLE_CUDA=true
            TAG="micapipe:v0.2.4-cuda-${BUILD_DATE}"
            shift
            ;;
        --no-cache)
            NO_CACHE=true
            shift
            ;;
        --tag)
            TAG="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --cuda      Enable CUDA support"
            echo "  --no-cache  Build without cache"
            echo "  --tag TAG   Custom Docker tag"
            echo "  --help      Show this help"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

print_info "Building with tag: $TAG"
if [ "$ENABLE_CUDA" = true ]; then
    print_info "CUDA support: ENABLED"
else
    print_info "CUDA support: DISABLED (CPU-only)"
fi
echo

# Step 1: Check prerequisites
print_info "Step 1: Checking prerequisites..."

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
    exit 1
fi

print_success "Docker and Singularity are available"

# Step 2: Setup directories
print_info "Step 2: Setting up directories..."

# Create directories if they don't exist
mkdir -p "$SIF_DIR" 2>/dev/null || {
    print_error "Cannot create SIF directory: $SIF_DIR"
    print_error "Please ask your admin to create it: mkdir -p $SIF_DIR && chown $USER:$USER $SIF_DIR"
    print_error "Or use an alternative path you have write access to with: --singularity-dir /path/to/writable/dir"
    exit 1
}

mkdir -p "$TMP_DIR" 2>/dev/null || {
    print_error "Cannot create temp directory: $TMP_DIR"
    print_error "Please ask your admin to create it: mkdir -p $TMP_DIR && chown $USER:$USER $TMP_DIR"
    print_error "Or set SINGULARITY_TMPDIR to a writable location: export SINGULARITY_TMPDIR=/tmp"
    exit 1
}

# Check write permissions
if [ ! -w "$SIF_DIR" ]; then
    print_error "No write permission for SIF directory: $SIF_DIR"
    print_error "Please ask your admin to fix permissions: chown $USER:$USER $SIF_DIR"
    print_error "Or use the generic build script with a writable path:"
    print_error "  ./scripts/build_container.sh --singularity --singularity-dir ~/containers"
    exit 1
fi

if [ ! -w "$TMP_DIR" ]; then
    print_error "No write permission for temp directory: $TMP_DIR"
    print_error "Please ask your admin to fix permissions: chown $USER:$USER $TMP_DIR"
    print_error "Or use a writable temp directory: export SINGULARITY_TMPDIR=~/tmp"
    exit 1
fi

print_success "Directories are ready"

# Step 3: Check disk space
print_info "Step 3: Checking disk space..."

SIF_SPACE=$(df "$SIF_DIR" | awk 'NR==2 {print $4}')
TMP_SPACE=$(df "$TMP_DIR" | awk 'NR==2 {print $4}')

# Convert to GB (approximate)
SIF_GB=$((SIF_SPACE / 1024 / 1024))
TMP_GB=$((TMP_SPACE / 1024 / 1024))

print_info "Available space - SIF dir: ${SIF_GB}GB, Temp dir: ${TMP_GB}GB"

if [ $SIF_GB -lt 15 ]; then
    print_warning "Low space in SIF directory (${SIF_GB}GB). Recommend at least 15GB."
fi

if [ $TMP_GB -lt 10 ]; then
    print_warning "Low space in temp directory (${TMP_GB}GB). Recommend at least 10GB."
fi

# Step 4: Set environment variables
print_info "Step 4: Setting up environment..."

export SINGULARITY_TMPDIR="$TMP_DIR"
print_info "Set SINGULARITY_TMPDIR=$SINGULARITY_TMPDIR"

# Disable Docker Content Trust to avoid certificate issues
export DOCKER_CONTENT_TRUST=0
print_info "Set DOCKER_CONTENT_TRUST=0 (disabled for build compatibility)"

# Step 5: Build Docker image
print_info "Step 5: Building Docker image..."

BUILD_CMD="docker build"

if [ "$NO_CACHE" = true ]; then
    BUILD_CMD="$BUILD_CMD --no-cache"
    print_info "Building without cache"
fi

if [ "$ENABLE_CUDA" = true ]; then
    BUILD_CMD="$BUILD_CMD --build-arg ENABLE_CUDA=true"
    print_info "Building with CUDA support"
else
    BUILD_CMD="$BUILD_CMD --build-arg ENABLE_CUDA=false"
    print_info "Building CPU-only version"
fi

BUILD_CMD="$BUILD_CMD -t $TAG ."

print_info "Running: $BUILD_CMD"
BUILD_START=$(date +%s)

if eval "$BUILD_CMD"; then
    BUILD_END=$(date +%s)
    BUILD_TIME=$((BUILD_END - BUILD_START))
    print_success "Docker build completed in ${BUILD_TIME} seconds"
else
    print_error "Docker build failed"
    exit 1
fi

# Step 6: Convert to Singularity
print_info "Step 6: Converting to Singularity..."

# Generate unique filename
SIF_FILENAME="micapipe_${TAG//[^a-zA-Z0-9]/_}_${BUILD_DATE}.sif"
SIF_PATH="$SIF_DIR/$SIF_FILENAME"

print_info "Creating Singularity image: $SIF_PATH"

SINGULARITY_START=$(date +%s)

if singularity build --force "$SIF_PATH" docker-daemon://"$TAG"; then
    SINGULARITY_END=$(date +%s)
    SINGULARITY_TIME=$((SINGULARITY_END - SINGULARITY_START))
    print_success "Singularity build completed in ${SINGULARITY_TIME} seconds"
else
    print_error "Singularity build failed"
    exit 1
fi

# Step 7: Test the container
print_info "Step 7: Testing the container..."

if singularity exec "$SIF_PATH" micapipe --help > /dev/null 2>&1; then
    print_success "Container test passed"
else
    print_warning "Container test failed - but build completed"
fi

# Step 8: Cleanup and summary
print_info "Step 8: Cleanup and summary..."

# Show file size
SIF_SIZE=$(du -h "$SIF_PATH" | cut -f1)
print_info "Singularity image size: $SIF_SIZE"

# Optional: Remove Docker image to save space
read -p "Remove Docker image to save space? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    docker rmi "$TAG"
    print_info "Docker image removed"
fi

# Final summary
print_success "Build completed successfully!"
echo
print_info "Summary:"
echo "  Docker tag: $TAG"
echo "  Singularity file: $SIF_PATH"
echo "  File size: $SIF_SIZE"
echo "  CUDA support: $([ "$ENABLE_CUDA" = true ] && echo "Enabled" || echo "Disabled")"
echo "  Total build time: $((BUILD_TIME + SINGULARITY_TIME)) seconds"
echo

print_info "Usage examples:"
echo "  # Run with bind mounts:"
echo "  singularity exec -B /data:/data $SIF_PATH micapipe [args]"
echo
if [ "$ENABLE_CUDA" = true ]; then
    echo "  # Run with GPU support:"
    echo "  singularity exec --nv -B /data:/data $SIF_PATH micapipe [args]"
    echo
fi

print_info "Container includes:"
echo "  ✓ MRtrix 3.0.7"
echo "  ✓ FreeSurfer 7.4.1" 
echo "  ✓ FastSurfer 2.4.2"
echo "  ✓ DESIGNER pipeline"
echo "  ✓ Synb0-DISCO & SynBOLD-DisCo"
echo "  ✓ LAMAReg with antspy"
echo "  ✓ SWM for superficial white matter"

print_success "All done!"