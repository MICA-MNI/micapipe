#!/bin/bash

# Server Docker Build and Test Script
# Comprehensive script for building and testing micapipe container on server

echo "🚀 MICApipe Server Build & Test Script"
echo "======================================"
echo "Date: $(date)"
echo "User: $(whoami)"
echo "Server: $(hostname)"
echo ""

# Configuration
CONTAINER_NAME="micapipe"
BUILD_TAG="latest"
TEST_TAG="test"
SIF_VERSION="v1-beta"
SIF_LOCATION="/opt/micapipe"  # Default location, can be overridden
LOG_DIR="./build_logs"
BACKUP_DIR="./backups"

# Check for custom SIF location from environment or command line
if [ ! -z "$MICAPIPE_SIF_PATH" ]; then
    SIF_LOCATION="$MICAPIPE_SIF_PATH"
fi

# Create necessary directories
mkdir -p "$LOG_DIR" "$BACKUP_DIR"

# Function to log with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/build_session.log"
}

# Function to check system requirements
check_system() {
    log "🔍 Checking system requirements..."
    
    # Check Docker
    if ! command -v docker &> /dev/null; then
        log "❌ Docker not found. Please install Docker first."
        exit 1
    fi
    
    # Check Docker daemon
    if ! docker info &> /dev/null; then
        log "❌ Docker daemon not running. Please start Docker."
        exit 1
    fi
    
    # Check disk space (need at least 10GB)
    AVAILABLE_SPACE=$(df . | awk 'NR==2 {print $4}')
    if [ "$AVAILABLE_SPACE" -lt 10485760 ]; then  # 10GB in KB
        log "⚠️  Warning: Less than 10GB disk space available. Build may fail."
    fi
    
    # Check memory
    AVAILABLE_MEM=$(free -m | awk 'NR==2{print $7}')
    if [ "$AVAILABLE_MEM" -lt 4096 ]; then  # 4GB
        log "⚠️  Warning: Less than 4GB memory available. Build may be slow."
    fi
    
    log "✅ System requirements check completed"
    log "   Docker version: $(docker --version)"
    log "   Available space: $(df -h . | awk 'NR==2 {print $4}')"
    log "   Available memory: $(free -h | awk 'NR==2{print $7}')"
}

# Function to set environment
setup_environment() {
    log "🔧 Setting up build environment..."
    
    # Disable Docker Content Trust to avoid certificate issues
    export DOCKER_CONTENT_TRUST=0
    log "   DOCKER_CONTENT_TRUST=0"
    
    # Set Docker buildkit for better build output
    export DOCKER_BUILDKIT=1
    log "   DOCKER_BUILDKIT=1"
    
    # Set timezone to avoid interactive prompts
    export DEBIAN_FRONTEND=noninteractive
    log "   DEBIAN_FRONTEND=noninteractive"
    
    log "✅ Environment setup completed"
}

# Function to backup existing containers/images
backup_existing() {
    log "💾 Backing up existing containers..."
    
    # Check if container exists
    if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
        log "   Found existing container: $CONTAINER_NAME"
        docker commit "$CONTAINER_NAME" "${CONTAINER_NAME}:backup-$(date +%Y%m%d_%H%M%S)" 2>&1 | tee -a "$LOG_DIR/backup.log"
        log "   Container backed up"
    fi
    
    # Check if image exists
    if docker images --format '{{.Repository}}:{{.Tag}}' | grep -q "^${CONTAINER_NAME}:${BUILD_TAG}$"; then
        log "   Found existing image: ${CONTAINER_NAME}:${BUILD_TAG}"
        docker tag "${CONTAINER_NAME}:${BUILD_TAG}" "${CONTAINER_NAME}:backup-$(date +%Y%m%d_%H%M%S)" 2>&1 | tee -a "$LOG_DIR/backup.log"
        log "   Image backed up"
    fi
    
    log "✅ Backup completed"
}

# Function to run FSL-specific test build
test_fsl_section() {
    log "🧪 Testing FSL section specifically..."
    
    # Create test Dockerfile that stops after FSL
    head -n 150 Dockerfile > Dockerfile.fsl-test
    echo "RUN echo 'FSL test completed successfully'" >> Dockerfile.fsl-test
    echo "CMD [\"/bin/bash\"]" >> Dockerfile.fsl-test
    
    # Build FSL test
    log "   Building FSL test section..."
    docker build \
        --progress=plain \
        --no-cache \
        -f Dockerfile.fsl-test \
        -t "${CONTAINER_NAME}:fsl-test" \
        . 2>&1 | tee "$LOG_DIR/fsl_test_build.log"
    
    FSL_EXIT_CODE=$?
    
    # Check FSL test result
    if [ $FSL_EXIT_CODE -eq 0 ]; then
        log "✅ FSL section builds successfully!"
        
        # Test FSL functionality
        log "   Testing FSL functionality..."
        docker run --rm "${CONTAINER_NAME}:fsl-test" bash -c "source /opt/fsl-6.0.2/etc/fslconf/fsl.sh && echo 'FSLDIR: $FSLDIR' && ls -la $FSLDIR/bin/bet" 2>&1 | tee "$LOG_DIR/fsl_test_function.log"
        
        # Cleanup test image
        docker rmi "${CONTAINER_NAME}:fsl-test" 2>/dev/null || true
        rm -f Dockerfile.fsl-test
        
        return 0
    else
        log "❌ FSL section failed to build"
        log "   Check $LOG_DIR/fsl_test_build.log for details"
        
        # Keep test files for debugging
        mv Dockerfile.fsl-test "$LOG_DIR/"
        
        return 1
    fi
}

# Function to run main build
run_main_build() {
    log "🏗️  Starting main container build..."
    
    # Determine CUDA support
    CUDA_ARG=""
    if command -v nvidia-smi &> /dev/null; then
        log "   NVIDIA GPU detected, enabling CUDA support"
        CUDA_ARG="--build-arg ENABLE_CUDA=true"
    else
        log "   No NVIDIA GPU detected, building CPU-only version"
        CUDA_ARG="--build-arg ENABLE_CUDA=false"
    fi
    
    # Run the build
    log "   Starting Docker build with options: $CUDA_ARG"
    docker build \
        --progress=plain \
        --no-cache \
        $CUDA_ARG \
        -t "${CONTAINER_NAME}:${BUILD_TAG}" \
        . 2>&1 | tee "$LOG_DIR/main_build.log"
    
    BUILD_EXIT_CODE=$?
    
    if [ $BUILD_EXIT_CODE -eq 0 ]; then
        log "✅ Main build completed successfully!"
        return 0
    else
        log "❌ Main build failed"
        log "   Exit code: $BUILD_EXIT_CODE"
        log "   Check $LOG_DIR/main_build.log for details"
        return 1
    fi
}

# Function to test the built container
test_container() {
    log "🧪 Testing built container..."
    
    # Basic container test
    log "   Testing basic container startup..."
    docker run --rm "${CONTAINER_NAME}:${BUILD_TAG}" echo "Container starts successfully" 2>&1 | tee "$LOG_DIR/container_test.log"
    
    if [ $? -ne 0 ]; then
        log "❌ Container startup test failed"
        return 1
    fi
    
    # Test key tools
    log "   Testing key neuroimaging tools..."
    
    # Test FSL
    docker run --rm "${CONTAINER_NAME}:${BUILD_TAG}" bash -c "source /opt/fsl-6.0.2/etc/fslconf/fsl.sh && bet" 2>&1 | grep -q "Usage" && log "   ✅ FSL: OK" || log "   ⚠️  FSL: Issue detected"
    
    # Test FreeSurfer
    docker run --rm "${CONTAINER_NAME}:${BUILD_TAG}" bash -c "source /opt/freesurfer-7.4.1/SetUpFreeSurfer.sh && mri_convert" 2>&1 | grep -q "usage" && log "   ✅ FreeSurfer: OK" || log "   ⚠️  FreeSurfer: Issue detected"
    
    # Test MRtrix
    docker run --rm "${CONTAINER_NAME}:${BUILD_TAG}" bash -c "mrinfo" 2>&1 | grep -q "usage" && log "   ✅ MRtrix: OK" || log "   ⚠️  MRtrix: Issue detected"
    
    # Test Python environment
    docker run --rm "${CONTAINER_NAME}:${BUILD_TAG}" python3 -c "import numpy, nibabel, matplotlib; print('Python packages OK')" 2>&1 | grep -q "OK" && log "   ✅ Python: OK" || log "   ⚠️  Python: Issue detected"
    
    log "✅ Container testing completed"
}

# Function to convert Docker to Singularity
convert_to_singularity() {
    log "🔄 Converting Docker container to Singularity..."
    
    # Check if Singularity is available
    if ! command -v singularity &> /dev/null; then
        log "⚠️  Singularity not found. Skipping .sif conversion."
        log "   Install Singularity to enable .sif creation"
        return 1
    fi
    
    # Create SIF location directory if it doesn't exist
    if [ ! -d "$SIF_LOCATION" ]; then
        log "📁 Creating SIF directory: $SIF_LOCATION"
        mkdir -p "$SIF_LOCATION" 2>/dev/null || {
            log "⚠️  Cannot create $SIF_LOCATION, using current directory"
            SIF_LOCATION="."
        }
    fi
    
    # Check if location is writable
    if [ ! -w "$SIF_LOCATION" ]; then
        log "⚠️  $SIF_LOCATION not writable, using current directory"
        SIF_LOCATION="."
    fi
    
    SIF_FILE="$SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif"
    
    log "   Converting to: $SIF_FILE"
    log "   This may take several minutes..."
    
    # Convert Docker to Singularity
    singularity build "$SIF_FILE" docker-daemon://${CONTAINER_NAME}:${BUILD_TAG} 2>&1 | tee "$LOG_DIR/singularity_build.log"
    
    if [ $? -eq 0 ]; then
        log "✅ Singularity conversion successful!"
        log "   SIF file created: $SIF_FILE"
        log "   Size: $(du -sh "$SIF_FILE" | cut -f1)"
        
        # Test the SIF file
        log "🧪 Testing Singularity container..."
        singularity exec "$SIF_FILE" echo "Singularity container test successful" 2>&1 | tee -a "$LOG_DIR/singularity_test.log"
        
        if [ $? -eq 0 ]; then
            log "✅ Singularity container test passed"
        else
            log "⚠️  Singularity container test failed"
        fi
        
        return 0
    else
        log "❌ Singularity conversion failed"
        log "   Check $LOG_DIR/singularity_build.log for details"
        return 1
    fi
}

# Function to generate summary report
generate_report() {
    log "📊 Generating build report..."
    
    REPORT_FILE="$LOG_DIR/build_report_$(date +%Y%m%d_%H%M%S).txt"
    
    cat > "$REPORT_FILE" << EOF
MICApipe Container Build Report
===============================
Date: $(date)
Server: $(hostname)
User: $(whoami)

Build Configuration:
- Container Name: $CONTAINER_NAME
- Build Tag: $BUILD_TAG
- CUDA Support: $(command -v nvidia-smi &> /dev/null && echo "Enabled" || echo "Disabled")

System Information:
- Docker Version: $(docker --version 2>/dev/null || echo "N/A")
- Available Disk Space: $(df -h . | awk 'NR==2 {print $4}' 2>/dev/null || echo "N/A")
- Available Memory: $(free -h | awk 'NR==2{print $7}' 2>/dev/null || echo "N/A")

Build Results:
$([ -f "$LOG_DIR/fsl_test_build.log" ] && echo "- FSL Test: $(grep -q "successfully" "$LOG_DIR/fsl_test_build.log" && echo "PASSED" || echo "FAILED")" || echo "- FSL Test: SKIPPED")
$([ -f "$LOG_DIR/main_build.log" ] && echo "- Main Build: $(tail -10 "$LOG_DIR/main_build.log" | grep -q "Successfully built" && echo "PASSED" || echo "FAILED")" || echo "- Main Build: NOT ATTEMPTED")
$([ -f "$LOG_DIR/container_test.log" ] && echo "- Container Test: $(grep -q "successfully" "$LOG_DIR/container_test.log" && echo "PASSED" || echo "FAILED")" || echo "- Container Test: NOT ATTEMPTED")

Tool Versions in Container:
$(docker run --rm "${CONTAINER_NAME}:${BUILD_TAG}" bash -c "
echo '- MRtrix: '$(mrinfo --version 2>/dev/null | head -1 || echo 'N/A')
echo '- FreeSurfer: '$(cat /opt/freesurfer-7.4.1/build-stamp.txt 2>/dev/null || echo 'N/A')
echo '- FSL: '$(cat /opt/fsl-6.0.2/etc/fslversion 2>/dev/null || echo 'N/A')
" 2>/dev/null || echo "- Tool versions: Could not retrieve (container may not be running)")

Log Files Generated:
$(ls -la "$LOG_DIR"/*.log 2>/dev/null | awk '{print "- " $9 " (" $5 " bytes)"}' || echo "- No log files found")

Singularity Container:
$(if [ -f "$SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif" ]; then
    echo "✅ Singularity container created: $SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif"
    echo "   Size: $(du -sh "$SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif" | cut -f1)"
    echo "   Usage: singularity exec $SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif [command]"
else
    echo "❌ Singularity container not created"
    if ! command -v singularity &> /dev/null; then
        echo "   Reason: Singularity not installed"
    else
        echo "   Reason: Build or conversion failed"
    fi
fi)

Next Steps:
$(if docker images | grep -q "${CONTAINER_NAME}.*${BUILD_TAG}"; then
    echo "✅ Container built successfully and ready for use"
    echo "   - Run Docker container: docker run -it ${CONTAINER_NAME}:${BUILD_TAG}"
    if [ -f "$SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif" ]; then
        echo "   - Run Singularity container: singularity exec $SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif [command]"
        echo "   - Interactive Singularity: singularity shell $SIF_LOCATION/${CONTAINER_NAME}_${SIF_VERSION}.sif"
    else
        echo "   - Manual Singularity conversion: singularity build ${CONTAINER_NAME}_${SIF_VERSION}.sif docker-daemon://${CONTAINER_NAME}:${BUILD_TAG}"
    fi
else
    echo "❌ Container build failed - check log files for troubleshooting"
    echo "   - Review main build log: $LOG_DIR/main_build.log"
    echo "   - Check FSL test log: $LOG_DIR/fsl_test_build.log"
fi)
EOF

    log "📋 Build report generated: $REPORT_FILE"
    
    # Display summary
    echo ""
    echo "📋 BUILD SUMMARY"
    echo "================"
    cat "$REPORT_FILE"
}

# Main execution function
main() {
    log "🚀 Starting MICApipe container build process..."
    
    # Pre-flight checks
    check_system
    setup_environment
    backup_existing
    
    # Test FSL section first (most problematic)
    if test_fsl_section; then
        log "🎯 FSL test passed, proceeding with full build..."
        
        # Run main build
        if run_main_build; then
            log "🎯 Main build successful, running tests..."
            test_container
            
            # Convert to Singularity
            log "🔄 Converting to Singularity format..."
            convert_to_singularity
        else
            log "❌ Main build failed, skipping tests and Singularity conversion"
        fi
    else
        log "❌ FSL test failed, check FSL installation section"
        log "🔧 Consider reviewing the Dockerfile FSL section for issues"
    fi
    
    # Generate final report
    generate_report
    
    log "🏁 Build process completed. Check logs in $LOG_DIR/"
}

# Command line options
case "${1:-}" in
    "fsl-only")
        check_system
        setup_environment
        test_fsl_section
        ;;
    "no-test")
        check_system
        setup_environment
        backup_existing
        run_main_build
        convert_to_singularity
        generate_report
        ;;
    "test-only")
        check_system
        setup_environment
        test_container
        ;;
    "singularity-only")
        check_system
        setup_environment
        convert_to_singularity
        ;;
    "help"|"-h"|"--help")
        echo "MICApipe Server Build & Test Script"
        echo "=================================="
        echo ""
        echo "Usage: $0 [OPTION]"
        echo ""
        echo "Options:"
        echo "  (no option)     Complete build, test, and Singularity conversion"
        echo "  fsl-only        Test only the FSL installation section"
        echo "  no-test         Build and convert to Singularity without testing"
        echo "  test-only       Test existing container without building"
        echo "  singularity-only Convert existing Docker container to Singularity"
        echo "  help            Show this help message"
        echo ""
        echo "Environment Variables:"
        echo "  MICAPIPE_SIF_PATH   Custom path for .sif file (default: /opt/micapipe)"
        echo ""
        echo "Output:"
        echo "  Docker container:   micapipe:latest"
        echo "  Singularity file:   \${MICAPIPE_SIF_PATH}/micapipe_${SIF_VERSION}.sif"
        echo "  Build logs:         ./build_logs/"
        echo ""
        echo "Examples:"
        echo "  $0                                  # Complete build and conversion"
        echo "  MICAPIPE_SIF_PATH=/data/containers $0  # Custom SIF location"
        echo "  $0 fsl-only                        # Test FSL section only"
        echo ""
        ;;
    *)
        main
        ;;
esac