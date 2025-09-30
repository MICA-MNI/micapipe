#!/bin/bash

# Micapipe Container Build Script - No Sudo Version
# This script builds micapipe containers using only user-writable directories

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

# User-writable fallback paths
HOME_CONTAINERS="$HOME/containers"
HOME_TMP="$HOME/tmp/singularity"
BUILD_DATE=$(date +%Y%m%d_%H%M%S)
DEFAULT_TAG="micapipe:v0.2.4-${BUILD_DATE}"

# Parse command line options
ENABLE_CUDA=false
NO_CACHE=false
TAG="$DEFAULT_TAG"
SIF_DIR=""
TMP_DIR=""

print_info "Micapipe Container Build - No Sudo Required"
print_info "Using user-writable directories as fallback"
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
        --sif-dir)
            SIF_DIR="$2"
            shift 2
            ;;
        --tmp-dir)
            TMP_DIR="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --cuda          Enable CUDA support"
            echo "  --no-cache      Build without cache"
            echo "  --tag TAG       Custom Docker tag"
            echo "  --sif-dir DIR   Custom SIF output directory"
            echo "  --tmp-dir DIR   Custom temporary directory"
            echo "  --help          Show this help"
            echo
            echo "This script uses user-writable directories when admin directories are not accessible."
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Determine directories to use
echo "🔍 Determining best directories to use..."

# Try preferred server paths first, fallback to user directories
PREFERRED_SIF="/data_/mica1/01_programs/singularity"
PREFERRED_TMP="/host/cassio/export03/data/enning/tmp"

# Set SIF directory
if [ -n "$SIF_DIR" ]; then
    # User specified custom directory
    print_info "Using custom SIF directory: $SIF_DIR"
elif [ -d "$PREFERRED_SIF" ] && [ -w "$PREFERRED_SIF" ]; then
    SIF_DIR="$PREFERRED_SIF"
    print_success "Using preferred SIF directory: $SIF_DIR"
else
    SIF_DIR="$HOME_CONTAINERS"
    print_warning "Preferred SIF directory not accessible, using: $SIF_DIR"
fi

# Set TMP directory
if [ -n "$TMP_DIR" ]; then
    # User specified custom directory
    print_info "Using custom TMP directory: $TMP_DIR"
elif [ -d "$PREFERRED_TMP" ] && [ -w "$PREFERRED_TMP" ]; then
    TMP_DIR="$PREFERRED_TMP"
    print_success "Using preferred TMP directory: $TMP_DIR"
else
    TMP_DIR="$HOME_TMP"
    print_warning "Preferred TMP directory not accessible, using: $TMP_DIR"
fi

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

# Create directories
if ! mkdir -p "$SIF_DIR" 2>/dev/null; then
    print_error "Cannot create SIF directory: $SIF_DIR"
    print_error "Try: $0 --sif-dir ~/my_containers"
    exit 1
fi

if ! mkdir -p "$TMP_DIR" 2>/dev/null; then
    print_error "Cannot create TMP directory: $TMP_DIR"
    print_error "Try: $0 --tmp-dir ~/my_tmp"
    exit 1
fi

# Test write permissions
if ! touch "$SIF_DIR/.test" 2>/dev/null; then
    print_error "No write permission for SIF directory: $SIF_DIR"
    exit 1
else
    rm -f "$SIF_DIR/.test"
fi

if ! touch "$TMP_DIR/.test" 2>/dev/null; then
    print_error "No write permission for TMP directory: $TMP_DIR"
    exit 1
else
    rm -f "$TMP_DIR/.test"
fi

print_success "Directories are ready"
print_info "SIF output: $SIF_DIR"
print_info "Temp files: $TMP_DIR"

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
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Try using a different directory: $0 --sif-dir /path/with/more/space"
        exit 1
    fi
fi

if [ $TMP_GB -lt 10 ]; then
    print_warning "Low space in temp directory (${TMP_GB}GB). Recommend at least 10GB."
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Try using a different directory: $0 --tmp-dir /path/with/more/space"
        exit 1
    fi
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

# Step 8: Final summary
print_info "Step 8: Final summary..."

# Show file size
SIF_SIZE=$(du -h "$SIF_PATH" | cut -f1)
print_info "Singularity image size: $SIF_SIZE"

# Cleanup option
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

print_info "Directory locations used:"
echo "  SIF files: $SIF_DIR"
echo "  Temp files: $TMP_DIR"
echo

print_success "All done! No sudo was required."