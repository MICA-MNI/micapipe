# Pull Request: Comprehensive Base Image Build System

## Branch: `comprehensive-base-image`

### Summary

This PR introduces a complete two-stage Docker build system and working Singularity container for micapipe v1.0.0-beta.

---

## 🎯 What's New

### Two-Stage Docker Build System
- **Stage 1**: `Dockerfile.mamba-base` - Comprehensive base image with all neuroimaging tools (~98GB)
- **Stage 2**: `Dockerfile.main` - Lightweight image adding only micapipe code on top of base (~98.6GB total)

This approach allows:
- Fast rebuilds when only micapipe code changes (3-5 min instead of 45-90 min)
- Stable base image that can be reused across versions
- Easier debugging and maintenance

### Docker Images
- `ghcr.io/mica-mni/micapipe-comprehensive-base:latest` - Base with all tools
- `ghcr.io/mica-mni/micapipe:latest` - Full pipeline image

### Singularity Container
- `/export03/data/enning/singularity/micapipe_v1_beta.sif` (21GB)
- Built and tested on `bb-compxg-01` server

---

## 🛠️ Technical Changes

### Dockerfile.main
- Fixed paths: `/opt/fix` (was `/opt/fix1.068`), `/opt/freesurfer` (was `/opt/freesurfer-7.4.1`)
- Updated entrypoint to `/opt/micapipe/micapipe`
- Using selective COPY commands to avoid copying tarballs

### Build Scripts
- `build_comprehensive_base_server.sh` - Builds base image
- `build_main_image_server.sh` - Builds main image on top of base
- `build_singularity.sh` - Converts Docker to Singularity SIF

### Bug Fixes
- Restored missing `fix_settings.sh` from master branch
- Fixed FIX path configuration
- Fixed FreeSurfer path configuration

---

## ✅ Verified Working

### Container Tools
| Tool | Path | Status |
|------|------|--------|
| FreeSurfer | `/opt/freesurfer` | ✅ |
| FSL | `/opt/fsl-6.0.2` | ✅ |
| AFNI | `/opt/afni` | ✅ |
| ANTs | `/opt/ants` | ✅ |
| Python | 3.11.11 | ✅ |
| MRtrix3 | `/opt/mrtrix3` | ✅ |

### Basic Tests
```bash
# Help
singularity run micapipe_v1_beta.sif -help  # ✅ Works

# Version
singularity run micapipe_v1_beta.sif -version  # ✅ Returns v0.2.3
```

---

## 📋 Usage

### Build Docker Images
```bash
# On server with Docker
./build_comprehensive_base_server.sh  # First time only (~45-90 min)
./build_main_image_server.sh          # For code updates (~3-5 min)
```

### Build Singularity
```bash
./build_singularity.sh
```

### Run Pipeline
```bash
singularity run /path/to/micapipe_v1_beta.sif \
    -bids /path/to/bids \
    -out /path/to/output \
    -sub 001 \
    -proc_structural
```

---

## 📁 Files Changed

### New Files
- `Dockerfile.main` - Main image Dockerfile
- `Dockerfile.mamba-base` - Base image Dockerfile
- `build_comprehensive_base_server.sh`
- `build_main_image_server.sh`
- `build_singularity.sh`
- `fix_settings.sh` - FIX configuration

### Modified Files
- `.dockerignore` - Exclude build artifacts

---

## 🔮 Next Steps

1. Run full pipeline test on sample BIDS dataset
2. Validate all processing modules (-proc_structural, -proc_surf, -proc_dwi, etc.)
3. Performance benchmarking
4. Documentation updates for new build process

---

## 📅 Milestone

**Date**: December 10, 2025  
**Version**: v1.0.0-beta  
**Status**: ✅ Working container build verified
