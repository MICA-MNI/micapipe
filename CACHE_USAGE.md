# MICApipe Container Caching System

This document explains how to use the local caching system for large dependencies (FSL and FreeSurfer) to speed up Docker builds and avoid repeated downloads.

## Overview

The caching system stores large dependencies locally on the server to avoid downloading them repeatedly during container builds. This is especially useful when:
- Building containers multiple times
- Network bandwidth is limited
- Downloads are slow or unreliable
- Working with large files (FSL ~600MB, FreeSurfer ~3GB)

## Setup

### 1. Cache Directory

Create a cache directory on your server:
```bash
sudo mkdir -p /data_/mica1/01_programs/micapipe_cache
sudo chmod 755 /data_/mica1/01_programs/micapipe_cache
```

### 2. Download Dependencies

Run the cache setup script to download dependencies:
```bash
chmod +x cache_dependencies.sh
./cache_dependencies.sh
```

This will download:
- FSL 6.0.2 (~600MB) to `/data_/mica1/01_programs/micapipe_cache/fsl-6.0.2-centos6_64.tar.gz`
- FreeSurfer 7.4.1 (~3GB) to `/data_/mica1/01_programs/micapipe_cache/freesurfer-linux-centos6_x86_64-7.4.1.tar.gz`

## Usage

### Building with Cache (Recommended)

```bash
# Default build with cache
./build_container.sh

# Build with CUDA and cache
./build_container.sh --cuda

# Custom cache directory
./build_container.sh --cache-dir /path/to/your/cache
```

### Building without Cache

```bash
# Disable caching (downloads files during build)
./build_container.sh --no-cache
```

## Build Script Options

```bash
Usage: ./build_container.sh [--cuda] [--no-cache] [--singularity-path PATH] [--cache-dir PATH]
  --cuda                Enable CUDA support
  --no-cache            Build without Docker cache
  --singularity-path    Custom path for .sif output (default: /data_/mica1/01_programs/singularity)
  --cache-dir           Path to dependency cache (default: /data_/mica1/01_programs/micapipe_cache)
```

## Cache Validation

The system validates cached files by checking:
- File existence
- Minimum file sizes (FSL: 500MB, FreeSurfer: 2GB)
- File integrity during download

If validation fails, the system will re-download the files.

## Cache Management

### Check Cache Status
```bash
ls -lh /data_/mica1/01_programs/micapipe_cache/
```

### Clear Cache
```bash
rm -f /data_/mica1/01_programs/micapipe_cache/*.tar.gz
```

### Update Cache
```bash
./cache_dependencies.sh  # Re-downloads and validates files
```

## Troubleshooting

### Cache Directory Not Accessible
If you see "Cache directory not accessible", check:
1. Directory exists: `ls -ld /data_/mica1/01_programs/micapipe_cache`
2. Write permissions: `touch /data_/mica1/01_programs/micapipe_cache/test && rm /data_/mica1/01_programs/micapipe_cache/test`

### Download Failures
If cache setup fails:
1. Check network connectivity: `curl -I https://fsl.fmrib.ox.ac.uk`
2. Check disk space: `df -h /data_/mica1/01_programs/micapipe_cache`
3. Retry with: `./cache_dependencies.sh`

### Build Issues
If builds still fail with cache:
1. Verify cache files: `ls -lh /data_/mica1/01_programs/micapipe_cache/`
2. Check Docker has access to cache mount
3. Try building without cache: `./build_container.sh --no-cache`

## Benefits

- **Faster Builds**: Eliminates 3.6GB+ of downloads per build
- **Reliability**: Avoids network timeouts and corruption
- **Bandwidth**: Reduces server bandwidth usage
- **Efficiency**: Multiple builds share cached files

## File Locations

- Cache Directory: `/data_/mica1/01_programs/micapipe_cache/`
- FSL Archive: `fsl-6.0.2-centos6_64.tar.gz` (~600MB)
- FreeSurfer Archive: `freesurfer-linux-centos6_x86_64-7.4.1.tar.gz` (~3GB)
- Cache Script: `./cache_dependencies.sh`
- Build Script: `./build_container.sh`