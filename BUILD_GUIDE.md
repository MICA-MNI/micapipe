# Micapipe Container Build Guide

This guide will help you build the micapipe container with your specific server configuration.

## Your Server Configuration

- **SIF Image Location**: `/data_/mica1/01_programs/singularity`
- **Temporary Files**: `/host/cassio/export03/data/enning/tmp`
- **Current Branch**: `docker-container-updates` (with latest updates)

## Prerequisites Check

Before building, ensure you have the required tools and permissions:

```bash
# Check Docker is installed and running
docker --version
docker info

# Check Singularity is available
singularity --version

# Check directory permissions
ls -la /data_/mica1/01_programs/singularity
ls -la /host/cassio/export03/data/enning/tmp

# Create directories if they don't exist
sudo mkdir -p /data_/mica1/01_programs/singularity
sudo mkdir -p /host/cassio/export03/data/enning/tmp

# Ensure you have write permissions
sudo chown $USER:$USER /data_/mica1/01_programs/singularity
sudo chown $USER:$USER /host/cassio/export03/data/enning/tmp
```

## Build Commands

### Option 1: Using Server Paths (requires admin setup)

If your admin has already set up the directories with proper permissions:

```bash
# CPU-only build
./build_server.sh

# CUDA-enabled build  
./build_server.sh --cuda
```

### Option 2: No Sudo Required (user directories)

If you don't have sudo rights or the server directories aren't accessible:

```bash
# CPU-only build using user directories
./build_no_sudo.sh

# CUDA-enabled build using user directories
./build_no_sudo.sh --cuda

# Custom directories
./build_no_sudo.sh --sif-dir ~/my_containers --tmp-dir ~/my_tmp
```

### Option 3: Generic Build Script

Use the original generic script with custom paths:

```bash
# Set your temporary directory
export SINGULARITY_TMPDIR=~/tmp

# Build with custom singularity directory
./scripts/build_container.sh \
  --tag micapipe:v0.2.4 \
  --singularity \
  --singularity-dir ~/containers
```

## Step-by-Step Build Process

### 1. Set Environment Variables

```bash
# Set temporary directory for Singularity builds
export SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/tmp

# Verify the setting
echo "Singularity temp dir: $SINGULARITY_TMPDIR"
```

### 2. Navigate to Project Directory

```bash
cd /path/to/micapipe  # Replace with your actual path
git checkout docker-container-updates
git pull origin docker-container-updates
```

### 3. Choose Your Build Method

**If directories are accessible (admin setup):**
```bash
./build_server.sh
```

**If no sudo rights (user directories):**
```bash
./build_no_sudo.sh
```

**Custom paths:**
```bash
./build_no_sudo.sh --sif-dir ~/containers --tmp-dir ~/tmp
```

### 4. Monitor Build Progress

The script will show:
- Prerequisites validation
- Docker build progress
- Singularity conversion progress
- Final file locations and sizes

## No Sudo Rights Solution

If you don't have sudo access on the server, use the `build_no_sudo.sh` script which:

### Features:
- **Automatically detects** if preferred directories are accessible
- **Falls back** to user-writable directories (`~/containers`, `~/tmp/singularity`)
- **Allows custom paths** via command line options
- **No admin privileges** required

### Usage:
```bash
# Basic build (auto-detects best directories)
./build_no_sudo.sh

# With custom directories
./build_no_sudo.sh --sif-dir ~/my_containers --tmp-dir ~/my_tmp

# CUDA enabled with custom paths
./build_no_sudo.sh --cuda --sif-dir /scratch/$USER/containers
```

### Directory Hierarchy:
1. **Preferred**: `/data_/mica1/01_programs/singularity` (if accessible)
2. **Fallback**: `~/containers` (always works)
3. **Custom**: Any path you specify

## Expected Output Locations

After successful build, you'll find:

### Docker Image
```bash
# List Docker images
docker images | grep micapipe
```

### Singularity Image
```bash
# Check the SIF file
ls -lh /data_/mica1/01_programs/singularity/micapipe_*.sif
```

## Testing the Built Container

### Test Docker Container
```bash
# Test CPU version
docker run -it --rm micapipe:v0.2.4-cpu micapipe --help

# Test CUDA version (if built)
docker run --gpus all -it --rm micapipe:v0.2.4-cuda micapipe --help
```

### Test Singularity Container
```bash
# Test CPU version
singularity exec /data_/mica1/01_programs/singularity/micapipe_*.sif micapipe --help

# Test with GPU (if CUDA enabled)
singularity exec --nv /data_/mica1/01_programs/singularity/micapipe_*.sif micapipe --help
```

## Disk Space Management

### Before Building
```bash
# Check available space
df -h /data_/mica1/01_programs/singularity
df -h /host/cassio/export03/data/enning/tmp

# Clean up old Docker images if needed
docker system prune -f
docker image prune -f
```

### During Build
The build process will use:
- **Docker space**: ~10-15 GB for layers and final image
- **Singularity temp**: ~5-10 GB during conversion
- **Final SIF**: ~8-12 GB depending on CUDA options

### After Build
```bash
# Clean up Docker images to save space (optional)
docker rmi micapipe:v0.2.4-cpu  # Keep only if you need Docker version

# Clean temporary files
rm -rf /host/cassio/export03/data/enning/tmp/singularity_*
```

## Troubleshooting

### Common Issues

#### 1. Permission Denied
```bash
# Fix directory permissions
sudo chown -R $USER:$USER /data_/mica1/01_programs/singularity
sudo chown -R $USER:$USER /host/cassio/export03/data/enning/tmp
```

#### 2. Out of Space
```bash
# Check space usage
du -sh /host/cassio/export03/data/enning/tmp/*
docker system df

# Clean up
docker system prune -a
rm -rf /host/cassio/export03/data/enning/tmp/tmp.*
```

#### 3. Singularity Build Fails
```bash
# Ensure SINGULARITY_TMPDIR is set and accessible
export SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/tmp
mkdir -p $SINGULARITY_TMPDIR
chmod 755 $SINGULARITY_TMPDIR
```

### Build Logs
The build script provides detailed logs. Save them if needed:
```bash
./scripts/build_container.sh [options] 2>&1 | tee build_log_$(date +%Y%m%d_%H%M%S).txt
```

## Complete Example Command

Here's a complete command for your setup:

```bash
#!/bin/bash

# Set environment
export SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/tmp

# Navigate to project
cd /path/to/micapipe
git checkout docker-container-updates

# Build container
./scripts/build_container.sh \
  --tag micapipe:v0.2.4-$(date +%Y%m%d) \
  --singularity \
  --singularity-dir /data_/mica1/01_programs/singularity \
  2>&1 | tee build_log_$(date +%Y%m%d_%H%M%S).txt

# Test the result
echo "Testing built container..."
LATEST_SIF=$(ls -t /data_/mica1/01_programs/singularity/micapipe_*.sif | head -1)
singularity exec "$LATEST_SIF" micapipe --help

echo "Build complete! SIF location: $LATEST_SIF"
```

## New Features in This Build

Your container will include:
- ✅ **MRtrix 3.0.7** (upgraded from 3.0.1)
- ✅ **FreeSurfer 7.4.1** (upgraded from 7.3.2)
- ✅ **FastSurfer 2.4.2** (upgraded and frozen)
- ✅ **DESIGNER pipeline** for diffusion MRI
- ✅ **Synb0-DISCO & SynBOLD-DisCo** for phase encoding
- ✅ **LAMAReg** with antspy for cross-modality registration
- ✅ **SWM** for superficial white matter analysis
- ✅ **Optional CUDA support** (use `--cuda` flag)

The build should take approximately 30-60 minutes depending on your server's performance and network speed.