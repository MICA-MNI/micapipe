# MICApipe Two-Stage Build Strategy

## 🎯 Overview

This document describes the **two-stage build strategy** for MICApipe, designed to achieve **95% faster CI builds** by separating heavy dependencies from frequently-changing code.

### Build Time Comparison

| Build Type | Traditional | Two-Stage |
|------------|-------------|-----------|
| **Full build with all tools** | 45-90 minutes | Stage 1: 45-90 min (once) |
| **Code-only rebuild** | 45-90 minutes | Stage 2: 3-5 min ✨ |
| **Time saved per CI run** | - | ~60 minutes! |

---

## 🏗️ Architecture

```
┌─────────────────────────────────────────────────┐
│  STAGE 1: Base Image (Dockerfile.base)         │
│  Build RARELY - only when dependencies change  │
├─────────────────────────────────────────────────┤
│  FROM ubuntu:bionic-20201119                    │
│  ├─ System packages                             │
│  ├─ FSL 6.0.2 (frozen)                          │
│  ├─ FreeSurfer 7.4.1 (frozen)                   │
│  ├─ FastSurfer 2.4.2 (frozen)                   │
│  ├─ AFNI (latest)                               │
│  ├─ ANTs 2.3.4                                  │
│  ├─ MRtrix3 3.0.7 (frozen)                      │
│  ├─ DESIGNER v2 pipeline                        │
│  ├─ LAMAReg + ANTsPy                            │
│  ├─ SWM (Superficial White Matter)             │
│  ├─ Synb0-DISCO + SynBOLD-DisCo                 │
│  ├─ Miniconda/Mamba + Python envs              │
│  ├─ R environment                               │
│  └─ FSL FIX, Workbench, c3d                    │
│                                                 │
│  Result: micapipe-base:latest (stable layer)   │
└─────────────────────────────────────────────────┘
                      ↓
┌─────────────────────────────────────────────────┐
│  STAGE 2: Main Image (Dockerfile.main)         │
│  Build FREQUENTLY - every code change          │
├─────────────────────────────────────────────────┤
│  FROM micapipe-base:latest                     │
│  ├─ MICApipe source code                       │
│  ├─ Processing scripts                         │
│  ├─ Configuration files                        │
│  ├─ Surface data & parcellations               │
│  └─ Integration with base tools                │
│                                                 │
│  Result: micapipe:latest (ready to use!)       │
└─────────────────────────────────────────────────┘
```

---

## 📋 New Requirements Implemented

All new requirements from the issue are implemented in **Dockerfile.base**:

- ✅ **MRtrix3 3.0.7** - Upgraded and frozen
- ✅ **MRtrix3 PATH fix** - Environment path configured correctly
- ✅ **FreeSurfer 7.4.1** - Upgraded and frozen
- ✅ **FastSurfer 2.4.2** - Upgraded, frozen, with model downloads
- ✅ **DESIGNER pipeline** - Added with separate conda environment
- ✅ **Synb0-DISCO + SynBOLD-DisCo** - Added for phase encoding correction
- ✅ **LAMAReg** - Added with ANTsPy dependencies
- ✅ **SWM** - Added for superficial white matter analysis
- ✅ **CUDA build option** - `--build-arg ENABLE_CUDA=true/false` (default: false)

### Note on neurodocker/startup.sh

**Status:** Retained (not removed)

**Reason:** The startup.sh script is essential for:
- Sourcing FSL/FreeSurfer environment variables
- Activating conda environments properly
- Setting up user context

**Future:** Can be modernized to use direct Docker ENV in future versions.

---

## 🚀 Quick Start Guide

### Prerequisites

1. **Local Machine** (`~/micapipe`): Your development environment
2. **Server** (`/host/cassio/export03/data/enning/downloads`): Build location with:
   - Unlimited storage space
   - Pre-downloaded neuroimaging tool files

### Step 1: Migrate Code to Server

On your local machine:

```bash
cd ~/micapipe
git checkout comprehensive-base-image
./migrate_comprehensive_base_to_server.sh
```

This will:
- Copy all source code to server
- Copy Dockerfile.base and Dockerfile.main
- Copy build scripts
- Verify pre-downloaded files exist

### Step 2: Build Base Image (One-Time)

SSH to server and build the base image:

```bash
ssh your-server
cd /host/cassio/export03/data/enning/downloads

# CPU-only build (default)
./build_base_image_server.sh

# OR with CUDA support
./build_base_image_server.sh --enable-cuda

# With custom tag and push to registry
./build_base_image_server.sh --tag 20251006 --registry ghcr.io/mica-mni --push
```

⏱️ **Expected time:** 45-90 minutes  
📦 **Result:** `micapipe-base:latest` with all tools

### Step 3: Build Main Image (Fast!)

Build the main micapipe image:

```bash
cd /host/cassio/export03/data/enning/downloads

# Quick build using local base image
./build_main_image_server.sh

# Using specific base image
./build_main_image_server.sh --base-image micapipe-base:20251006

# With custom tag and push
./build_main_image_server.sh --tag v0.2.3-dev --push
```

⏱️ **Expected time:** 3-5 minutes (95% faster!)  
📦 **Result:** `micapipe:latest` ready to use

---

## 🔄 CI Workflow

### For Code Changes (Most Common)

```bash
# 1. Update code locally
cd ~/micapipe
git pull origin comprehensive-base-image
# ... make your changes ...

# 2. Migrate to server
./migrate_comprehensive_base_to_server.sh

# 3. SSH to server
ssh your-server
cd /host/cassio/export03/data/enning/downloads

# 4. Fast rebuild (3-5 min!)
./build_main_image_server.sh --push
```

### For Dependency Changes (Rare)

When you need to update FSL, FreeSurfer, Python packages, etc.:

```bash
# 1. Update Dockerfile.base locally
cd ~/micapipe
# ... edit Dockerfile.base ...

# 2. Migrate to server
./migrate_comprehensive_base_to_server.sh

# 3. SSH to server
ssh your-server
cd /host/cassio/export03/data/enning/downloads

# 4. Rebuild base (45-90 min)
./build_base_image_server.sh --push

# 5. Then rebuild main (3-5 min)
./build_main_image_server.sh --push
```

---

## 📁 File Structure

### On Local Machine (`~/micapipe/`)

```
micapipe/
├── Dockerfile.base                      # Stage 1: All neuroimaging tools
├── Dockerfile.main                      # Stage 2: MICApipe code only
├── build_base_image_server.sh          # Script to build base image
├── build_main_image_server.sh          # Script to build main image
├── migrate_comprehensive_base_to_server.sh  # Migration script
├── TWO_STAGE_BUILD_GUIDE.md           # This file
├── DOCKERFILE_COMPREHENSIVE_REVIEW.md  # Detailed review
└── ... (micapipe source code)
```

### On Server (`/host/cassio/export03/data/enning/downloads/`)

```
downloads/
├── Dockerfile.base                      # Copied from local
├── Dockerfile.main                      # Copied from local
├── build_base_image_server.sh          # Copied from local
├── build_main_image_server.sh          # Copied from local
├── micapipe/                            # Source code
├── functions/                           # Processing scripts
├── fsl_conf/                            # FSL configuration
├── R_config/                            # R setup
├── surfaces/                            # Surface data
├── parcellations/                       # Parcellation data
└── Pre-downloaded files:
    ├── fsl-6.0.2-centos6_64.tar.gz                      (~1 GB)
    ├── freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz    (~3 GB)
    ├── afni-linux_openmp_64.tgz                         (~500 MB)
    ├── fix-1.068.tar.gz                                 (~50 MB)
    └── Miniconda3-py39_22.11.1-1-Linux-x86_64.sh       (~100 MB)
```

---

## 🔧 Build Script Options

### build_base_image_server.sh

```bash
./build_base_image_server.sh [OPTIONS]

OPTIONS:
  --enable-cuda    Build with CUDA support (default: CPU-only)
  --tag TAG        Custom tag for image (default: date-based YYYYMMDD)
  --registry REG   Registry to push to (default: ghcr.io/mica-mni)
  --push           Automatically push to registry after build

EXAMPLES:
  # CPU-only build
  ./build_base_image_server.sh

  # GPU-enabled build
  ./build_base_image_server.sh --enable-cuda

  # Build and push with custom tag
  ./build_base_image_server.sh --tag 20251006-v1 --push
```

### build_main_image_server.sh

```bash
./build_main_image_server.sh [OPTIONS]

OPTIONS:
  --base-image IMAGE   Base image to use (default: ghcr.io/mica-mni/micapipe-base:latest)
  --tag TAG            Custom tag for main image (default: v0.2.3-YYYYMMDD)
  --registry REG       Registry to push to (default: ghcr.io/mica-mni)
  --push               Automatically push to registry after build

EXAMPLES:
  # Quick build using latest base
  ./build_main_image_server.sh

  # Use specific base image
  ./build_main_image_server.sh --base-image micapipe-base:20251006

  # Build and push with version tag
  ./build_main_image_server.sh --tag v0.2.3 --push
```

---

## 🧪 Testing

### Verify Base Image

```bash
# Check image exists
docker images micapipe-base

# Test basic functionality
docker run --rm micapipe-base:latest bash -c "which fsl && which freesurfer && which mrtrix3"

# Check MRtrix3 version
docker run --rm micapipe-base:latest bash -c "conda activate micapipe && mrinfo -version"
```

### Verify Main Image

```bash
# Check image exists
docker images micapipe

# Test micapipe help
docker run --rm micapipe:latest --help

# Test processing (if you have test data)
docker run --rm -v /path/to/data:/data micapipe:latest -sub 001 -ses 01
```

---

## 🐛 Troubleshooting

### Base Build Issues

**Problem:** "Base image build failed - out of memory"  
**Solution:**
```bash
# Increase Docker memory limit
./build_base_image_server.sh
# The script already uses --memory=12g --memory-swap=16g
# If still failing, close other processes or use a machine with more RAM
```

**Problem:** "FSL/FreeSurfer download timeout"  
**Solution:**
```bash
# Use pre-downloaded files! Check they exist:
ls -lh /host/cassio/export03/data/enning/downloads/*.tar.gz

# If missing, download manually:
cd /host/cassio/export03/data/enning/downloads
wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.2-centos6_64.tar.gz
# Then re-run build
```

### Main Build Issues

**Problem:** "Base image not found"  
**Solution:**
```bash
# Check if base image exists locally
docker images micapipe-base

# If not, either:
# 1. Build it: ./build_base_image_server.sh
# 2. Pull it: docker pull ghcr.io/mica-mni/micapipe-base:latest
# 3. Specify different base: ./build_main_image_server.sh --base-image <image>
```

**Problem:** "Missing micapipe source files"  
**Solution:**
```bash
# Re-run migration from local machine
cd ~/micapipe
./migrate_comprehensive_base_to_server.sh
```

---

## 📊 Performance Benchmarks

### Expected Build Times

| Component | Time | Frequency |
|-----------|------|-----------|
| **Base image (CPU)** | 45-60 min | Monthly |
| **Base image (GPU)** | 60-90 min | Monthly |
| **Main image** | 3-5 min | Daily/Per commit |
| **Code-only rebuild** | 2-3 min | Multiple per day |

### Resource Requirements

**Base Build:**
- CPU: 4+ cores recommended
- RAM: 12+ GB required
- Disk: 20+ GB for image
- Network: Stable connection for downloads

**Main Build:**
- CPU: 2+ cores sufficient
- RAM: 4+ GB sufficient
- Disk: 5+ GB for image
- Network: Not required if base exists locally

---

## 🔒 Security & Best Practices

### Container Registry

Push base images to a private registry:

```bash
# GitHub Container Registry (recommended)
export MICAPIPE_REGISTRY=ghcr.io/mica-mni
./build_base_image_server.sh --push

# Docker Hub
export MICAPIPE_REGISTRY=docker.io/yourusername
./build_base_image_server.sh --push
```

### Version Pinning

All critical dependencies are frozen:
- FSL 6.0.2
- FreeSurfer 7.4.1
- FastSurfer 2.4.2
- MRtrix3 3.0.7
- Python packages with specific versions

### Reproducibility

```bash
# Tag base images with dates for reproducibility
./build_base_image_server.sh --tag $(date +%Y%m%d)

# Reference specific base in CI
./build_main_image_server.sh --base-image micapipe-base:20251006
```

---

## 📚 Additional Resources

- **Comprehensive Review:** `DOCKERFILE_COMPREHENSIVE_REVIEW.md`
- **Original Dockerfile:** `Dockerfile` (v1 branch - single-stage)
- **Migration Script:** `migrate_comprehensive_base_to_server.sh`
- **GitHub Issues:** #133 (CUDA support request)

---

## ✅ Next Steps

1. ✅ Review this guide
2. ✅ Run migration: `./migrate_comprehensive_base_to_server.sh`
3. ✅ Build base: `./build_base_image_server.sh`
4. ✅ Build main: `./build_main_image_server.sh`
5. ✅ Test image: `docker run --rm micapipe:latest --help`
6. ✅ Set up CI to use fast main builds
7. ✅ Celebrate 95% faster builds! 🎉

---

**Questions or Issues?**

Contact: MICA Lab <mica.lab@mcgill.ca>  
Repository: https://github.com/MICA-MNI/micapipe
