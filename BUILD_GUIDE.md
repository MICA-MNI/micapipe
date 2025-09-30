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

### Option 1: CPU-Only Container (Default)

Build a CPU-only container for maximum compatibility:

```bash
# Navigate to the micapipe project
cd /path/to/micapipe

# Build CPU-only container with custom paths
./scripts/build_container.sh \
  --tag micapipe:v0.2.4-cpu \
  --singularity \
  --singularity-dir /data_/mica1/01_programs/singularity

# Set temporary directory for Singularity
export SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/tmp
```

### Option 2: CUDA-Enabled Container (If you have GPU)

Build a CUDA-enabled container for GPU acceleration:

```bash
# Build CUDA-enabled container
./scripts/build_container.sh \
  --cuda \
  --tag micapipe:v0.2.4-cuda \
  --singularity \
  --singularity-dir /data_/mica1/01_programs/singularity

# Set temporary directory for Singularity
export SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/tmp
```

### Option 3: Clean Build (No Cache)

If you want to ensure a completely fresh build:

```bash
# Clean build without Docker cache
./scripts/build_container.sh \
  --no-cache \
  --tag micapipe:v0.2.4-fresh \
  --singularity \
  --singularity-dir /data_/mica1/01_programs/singularity

export SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/tmp
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

### 3. Run the Build

```bash
# Example: Build CPU-only container
./scripts/build_container.sh \
  --tag micapipe:v0.2.4-$(date +%Y%m%d) \
  --singularity \
  --singularity-dir /data_/mica1/01_programs/singularity
```

### 4. Monitor Build Progress

The script will show:
- Prerequisites validation
- Docker build progress
- Singularity conversion progress
- Final file locations and sizes

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