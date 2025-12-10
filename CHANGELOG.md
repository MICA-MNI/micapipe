# Changelog

All notable changes to micapipe will be documented in this file.

## [v1.0.0-beta] - 2025-12-10

### 🎉 Milestone: Working Container Build

**Branch**: `comprehensive-base-image`

### Added
- Two-stage Docker build system for faster iteration
  - `Dockerfile.mamba-base`: Comprehensive base image with all neuroimaging tools
  - `Dockerfile.main`: Lightweight overlay adding micapipe code
- Build scripts for server deployment
  - `build_comprehensive_base_server.sh`
  - `build_main_image_server.sh`
  - `build_singularity.sh`
- `fix_settings.sh` for FIX classifier configuration

### Fixed
- FreeSurfer path: `/opt/freesurfer` (was `/opt/freesurfer-7.4.1`)
- FIX path: `/opt/fix` (was `/opt/fix1.068`)
- Entrypoint configuration in Dockerfile.main
- Selective COPY commands to avoid including build tarballs

### Container Details
- **Docker Image**: `ghcr.io/mica-mni/micapipe:latest` (98.6 GB)
- **Singularity SIF**: `micapipe_v1_beta.sif` (21 GB)
- **Python**: 3.11.11
- **Build Server**: bb-compxg-01

### Verified Tools
- FreeSurfer 7.4.1
- FSL 6.0.2
- AFNI
- ANTs
- MRtrix3
- Workbench
- FIX

---

## [v0.2.3] - Previous Release

See original repository for previous changelog entries.
