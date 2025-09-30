# Docker Build Memory Issues - Troubleshooting Guide

## Problem: Exit Code 137 (Memory Killed)

When you see exit code 137 during Docker build, it means the process was killed due to insufficient memory. This commonly happens during large file downloads and extractions like FreeSurfer.

## Error Pattern
```
#11 8.337 Downloading FreeSurfer ...
#11 ERROR: executor failed running [...]: exit code: 137
```

## Solutions

### 1. Use Memory-Safe Build Script (Recommended)
```bash
./memory_safe_build.sh
```

This script:
- Checks system memory and swap
- Uses conservative memory limits
- Provides detailed error diagnostics
- Handles cleanup automatically

### 2. Increase Docker Memory Limits
```bash
docker build \
    --memory=8g \
    --memory-swap=12g \
    --shm-size=2g \
    -t micapipe:v1-beta .
```

### 3. Enable System Swap Space

**Linux:**
```bash
# Create 8GB swap file
sudo fallocate -l 8G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile

# Make permanent
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

**macOS:**
Swap is managed automatically, but you can:
- Increase Docker Desktop memory allocation
- Close other memory-intensive applications

### 4. Docker Desktop Settings

1. Open Docker Desktop
2. Go to Settings → Resources → Advanced
3. Increase Memory allocation to 8GB+
4. Increase Swap to 4GB+
5. Apply & Restart

### 5. Build in Stages (Advanced)

For debugging or very constrained systems:
```bash
# Build base neurodocker image first
docker build --target neurodocker -t micapipe:base .

# Then continue with full build
docker build -t micapipe:v1-beta .
```

## Updated Dockerfile Features

The Dockerfile now includes:
- **Memory-efficient FreeSurfer installation**: Downloads to temp file, then extracts
- **Multiple download fallbacks**: FTP → HTTP → Mirror → Graceful failure
- **Robust error handling**: Continues build even if some downloads fail
- **Cleanup after extraction**: Removes large temp files immediately

## Monitoring Build Progress

During build, monitor system resources:
```bash
# Terminal 1: Start build
./memory_safe_build.sh

# Terminal 2: Monitor resources
watch -n 1 'free -h && docker stats --no-stream'
```

## Prevention

- **Minimum Requirements**: 16GB RAM, 8GB swap, 50GB disk space
- **Optimal Setup**: 32GB RAM, 16GB swap, 100GB SSD
- **Network**: Stable connection for large downloads (3-5GB total)

## Getting Help

If memory issues persist:
1. Check build logs in `build_logs/`
2. Run `./diagnose_build_failure.sh` for detailed analysis
3. Try `./retry_build.sh` for automated fixes
4. Consider using a cloud instance with more memory

## Cloud Alternative

For persistent memory issues, consider building on a cloud instance:
```bash
# AWS EC2: c5.2xlarge (8 vCPU, 16GB RAM)
# Google Cloud: n1-standard-4 (4 vCPU, 15GB RAM)  
# Azure: Standard_D4s_v3 (4 vCPU, 16GB RAM)
```