# MICApipe Server Deployment Guide

## Quick Start on Server

### 1. Upload Files to Server
```bash
# Upload the updated branch to your server
git clone https://github.com/MICA-MNI/micapipe.git
cd micapipe
git checkout docker-container-updates
```

### 2. Run Server Build Script
```bash
# For complete build, test, and Singularity conversion (recommended)
./server_build_test.sh

# For simple build without extensive testing
./build_container.sh

# For help and options
./server_build_test.sh help
```

### 3. Set Custom Singularity Location (Optional)
```bash
# Set custom location for .sif file
export MICAPIPE_SIF_PATH="/your/custom/path"
./server_build_test.sh

# Or inline
MICAPIPE_SIF_PATH="/data/containers" ./server_build_test.sh
```

## Build Script Options

### **server_build_test.sh** (Comprehensive - Recommended)
- Complete build pipeline with testing and validation
- Automatic Singularity conversion
- Comprehensive logging and reporting
- Multiple execution modes

### **build_container.sh** (Simple - Quick Singularity Build)
- Direct Singularity .sif file creation
- Docker container as intermediate step (automatically removed)
- No sudo requirements (automatically detects user directories)
- CUDA auto-detection
- Minimal logging

### Alternative Direct Commands
```bash
# Build Singularity with CUDA support
./build_container.sh --cuda

# Build to custom location
./build_container.sh --output /data/containers

# Build without cache (slower but fresh build)
./build_container.sh --no-cache
```

### 3. Alternative Build Options

#### Option A: Comprehensive Build (Recommended)
```bash
# Complete build, test, and Singularity conversion
./server_build_test.sh
```

#### Option B: Simple Singularity Build (Faster)
```bash
# Direct .sif creation with automatic cleanup
./build_container.sh
```

#### Option C: Direct Docker Build
```bash
# Set environment to avoid certificate issues
export DOCKER_CONTENT_TRUST=0

# Build with CUDA support (if available)
docker build --build-arg ENABLE_CUDA=true -t micapipe:latest .

# Or build CPU-only version
docker build --build-arg ENABLE_CUDA=false -t micapipe:latest .
```

#### Option D: Direct Docker Build (Advanced)
```bash
# Set environment to avoid certificate issues
export DOCKER_CONTENT_TRUST=0

# Build with CUDA support (if available)
docker build --build-arg ENABLE_CUDA=true -t micapipe:latest .

# Or build CPU-only version
docker build --build-arg ENABLE_CUDA=false -t micapipe:latest .

# Then manually convert to Singularity
singularity build /data_/mica1/01_programs/singularity/micapipe_v1-beta.sif docker-daemon://micapipe:latest
```

## Troubleshooting Common Server Issues

### Issue 1: Docker Content Trust Certificate Errors
```bash
export DOCKER_CONTENT_TRUST=0
```

### Issue 2: No Sudo Rights
- All build scripts automatically detect writable directories
- No manual configuration needed for non-sudo environments

### Issue 3: Build Process Concerns
- Use `./build_container.sh` for direct Singularity creation
- Intermediate Docker containers are automatically cleaned up
- Focus on .sif file as primary output format

### Issue 4: Limited Disk Space
```bash
# Clean up Docker before building
docker system prune -a -f

# Check available space
df -h .
```

### Issue 5: Memory Issues During Build
```bash
# Build with reduced parallelism
docker build --memory=4g --memory-swap=8g -t micapipe:latest .
```

## Verification Steps

### Check Tool Versions in Built Container
```bash
# Verify all tools are correctly installed
docker run --rm micapipe:latest bash -c "
echo 'MRtrix:' && mrinfo --version | head -1
echo 'FreeSurfer:' && cat /opt/freesurfer-7.4.1/build-stamp.txt
echo 'FSL:' && cat /opt/fsl-6.0.2/etc/fslversion
echo 'Python packages:' && python3 -c 'import numpy, nibabel, matplotlib; print(\"OK\")'
"
```

### Test Container Functionality
```bash
# Basic test
docker run --rm micapipe:latest echo "Container working"

# Test neuroimaging tools
docker run --rm micapipe:latest bash -c "
source /opt/fsl-6.0.2/etc/fslconf/fsl.sh
source /opt/freesurfer-7.4.1/SetUpFreeSurfer.sh
echo 'FSL FSLDIR:' \$FSLDIR
echo 'FreeSurfer:' \$FREESURFER_HOME
"
```

## Build Logs and Monitoring

All build scripts generate comprehensive logs in `./build_logs/`:
- `build_session.log` - Complete build session log
- `fsl_test_build.log` - FSL-specific test results
- `main_build.log` - Main Docker build output
- `container_test.log` - Container functionality tests
- `build_report_TIMESTAMP.txt` - Final summary report

## Expected Tool Versions

The container includes these updated versions:
- **MRtrix**: 3.0.7 (upgraded from 3.0.1)
- **FreeSurfer**: 7.4.1 (upgraded from 7.3.2)
- **FastSurfer**: v2.4.2 (upgraded)
- **FSL**: 6.0.2 (with improved Python handling)
- **DESIGNER**: Latest from DESIGNER-v2 repository
- **Synb0-DISCO & SynBOLD-DisCo**: For phase encoding
- **LAMAReg**: With antspy dependencies
- **SWM**: Superficial white matter analysis

## Converting to Singularity

The server build script automatically creates a Singularity container:

### Automatic Conversion (Default)
```bash
# The script automatically creates: /data_/mica1/01_programs/singularity/micapipe_v1-beta.sif
./server_build_test.sh
```

### Custom Location
```bash
# Specify custom location for .sif file
MICAPIPE_SIF_PATH="/data/containers" ./server_build_test.sh

# Creates: /data/containers/micapipe_v1-beta.sif
```

### Manual Conversion (if needed)
```bash
# Only if automatic conversion failed
singularity build /data_/mica1/01_programs/singularity/micapipe_v1-beta.sif docker-daemon://micapipe:latest
```

### Using the Singularity Container
```bash
# Interactive session
singularity shell /data_/mica1/01_programs/singularity/micapipe_v1-beta.sif

# Run command
singularity exec /data_/mica1/01_programs/singularity/micapipe_v1-beta.sif [command]

# With data binding
singularity exec --bind /data:/data /data_/mica1/01_programs/singularity/micapipe_v1-beta.sif [command]
```

## Support

If you encounter issues:
1. Check the build logs in `./build_logs/`
2. Review the build report for specific failure points
3. Try the alternative build options listed above
4. Use the FSL fix script if needed: `./fix_fsl_build.sh`