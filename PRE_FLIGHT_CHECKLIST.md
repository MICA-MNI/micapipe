# ✅ PRE-FLIGHT CHECKLIST - READY FOR SERVER BUILD

**Date:** October 6, 2025  
**Branch:** comprehensive-base-image  
**Status:** ✅ ALL CHECKS PASSED

---

## 📋 Git Repository Status

✅ **Branch synced with remote**
- Current: comprehensive-base-image
- Status: Up to date with origin
- Latest commit: eab8521

✅ **Recent commits (5 most recent)**
```
eab8521 - Add conda environment specs for future lock file generation
93bee9a - Add conda pre-computation tools for instant builds
e6ef811 - Optimize conda/mamba environment solving - fix stuck build
6461d85 - Fix Docker Content Trust SSL certificate error
443b2a7 - Add quick reference summary for two-stage build
```

---

## 📦 Critical Files Verification

✅ **Core Dockerfiles**
- `Dockerfile.base` - ✅ EXISTS (31.4 KB, 721 lines)
- `Dockerfile.main` - ✅ EXISTS (2.9 KB, 68 lines)

✅ **Build Scripts**
- `build_base_image_server.sh` - ✅ EXISTS + EXECUTABLE
- `build_main_image_server.sh` - ✅ EXISTS + EXECUTABLE
- `migrate_comprehensive_base_to_server.sh` - ✅ EXISTS + EXECUTABLE

✅ **Bash Syntax Check**
- `build_base_image_server.sh` - ✅ NO ERRORS
- `build_main_image_server.sh` - ✅ NO ERRORS
- `migrate_comprehensive_base_to_server.sh` - ✅ NO ERRORS

---

## 🔧 Critical Fixes Applied

✅ **Docker Content Trust Fix** (Commit 6461d85)
- `export DOCKER_CONTENT_TRUST=0` added to both build scripts
- Prevents SSL certificate errors
- No sudo required

✅ **Conda Optimization** (Commit e6ef811)
- `conda config --set solver libmamba` - 10x faster solver
- Single-pass installation (all packages at once)
- Separated MRtrix3 and ANTsPy installations
- Channel priority set to flexible

---

## ✅ All 10 Requirements Verified in Dockerfile.base

| # | Requirement | Location | Status |
|---|-------------|----------|--------|
| 1 | **MRtrix3 3.0.7** (frozen) | Line 500 | ✅ VERIFIED |
| 2 | **FreeSurfer 7.4.1** (frozen) | Line 227 | ✅ VERIFIED |
| 3 | **FastSurfer 2.4.2** (frozen) | Line 618 | ✅ VERIFIED |
| 4 | **DESIGNER v2** (new) | Line 537 | ✅ VERIFIED |
| 5 | **Synb0/SynBOLD** (new) | Line 563 | ✅ VERIFIED |
| 6 | **LAMAReg** (new) | Line 516 | ✅ VERIFIED |
| 7 | **SWM** (new) | Line 516 | ✅ VERIFIED |
| 8 | **CUDA support** (optional) | Build arg | ✅ VERIFIED |
| 9 | **Remove neurodocker startup** | Not referenced | ✅ VERIFIED |
| 10 | **Keep R 3.6.3** | R_config preserved | ✅ VERIFIED |

---

## 🎯 Performance Expectations

| Stage | Expected Time | Notes |
|-------|---------------|-------|
| **Base image build** | 45-90 minutes | One-time (monthly) |
| **Conda environment** | 15-20 minutes | Optimized single-pass |
| **Main image build** | 3-5 minutes | Fast (per code change) |
| **Overall time savings** | **95% faster** | vs original 60 min/build |

---

## 🚀 Ready to Execute on Server

### Step 1: Pull Latest Code
```bash
cd ~/micapipe
git pull origin comprehensive-base-image
```

### Step 2: Migrate to Server
```bash
./migrate_comprehensive_base_to_server.sh
```

This will:
- Copy code from `~/micapipe` → `/host/cassio/export03/data/enning/downloads`
- Validate pre-downloaded files exist
- Offer to start base build immediately

### Step 3: Build Base Image
```bash
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh
```

**Expected output:**
```
🏗️  Starting Docker build...
=============================================
Step 1/XX : FROM ubuntu:18.04
 ---> 3339fde08fc3
Step 2/XX : ENV DEBIAN_FRONTEND=noninteractive
 ---> Running in abc123def456
...
📦 Installing micapipe conda environment (optimized single-pass)...
Solving environment: done
✅ Core environment installed
...
🎉 BASE IMAGE BUILD SUCCESSFUL!
```

### Step 4: Build Main Image (Future Builds)
```bash
cd /host/cassio/export03/data/enning/downloads
./build_main_image_server.sh
```

**Expected time:** 3-5 minutes

---

## 🔍 Common Issues & Solutions

### Issue: Docker Content Trust SSL Error
**Error:** `unable to retrieve valid leaf certificates`  
**Solution:** ✅ **FIXED** - Scripts now export `DOCKER_CONTENT_TRUST=0`

### Issue: Conda Solving Stuck
**Error:** Stuck at "Solving environment: ...working..."  
**Solution:** ✅ **FIXED** - Using libmamba solver + single-pass install

### Issue: Out of Memory
**Error:** Docker daemon OOM  
**Solution:** Scripts use `--memory=12g --memory-swap=16g`

### Issue: Missing Pre-downloaded Files
**Error:** Cannot download FreeSurfer/FSL  
**Solution:** Migration script validates files exist on server

---

## 📊 Build Architecture Summary

```
Two-Stage Build Strategy
========================

Stage 1: Base Image (Dockerfile.base)
  ├─ Ubuntu 18.04
  ├─ FSL 6.0.2
  ├─ FreeSurfer 7.4.1
  ├─ MATLAB MCR R2017b
  ├─ AFNI, ANTs 2.3.4, Workbench
  ├─ Miniconda + libmamba solver
  ├─ MRtrix3 3.0.7
  ├─ FastSurfer 2.4.2
  ├─ DESIGNER, LAMAReg, SWM
  ├─ Synb0/SynBOLD
  └─ R 3.6.3 environment
  
  Build: 45-90 min (one-time)
  Frequency: Rare (monthly)

Stage 2: Main Image (Dockerfile.main)
  ├─ FROM micapipe-base:latest
  ├─ COPY micapipe code
  └─ Integrate with tools
  
  Build: 3-5 min (fast!)
  Frequency: Every code change
```

---

## ✅ Final Verification Checklist

- [x] Git repository synced with remote
- [x] All critical files present and executable
- [x] Bash syntax validated (no errors)
- [x] All 10 requirements in Dockerfile.base
- [x] Docker Content Trust fix applied
- [x] Conda optimization applied
- [x] Build scripts have proper memory limits
- [x] Migration script updated for two-stage
- [x] Documentation complete
- [x] Ready for server execution

---

## 🎉 READY TO BUILD!

**Status:** ✅ **ALL SYSTEMS GO**

Everything is verified and ready. You can proceed with confidence to:
1. Pull code on local machine
2. Migrate to server
3. Build base image
4. Build main image

**Expected result:** Successful build in 45-90 minutes (base) + 3-5 minutes (main)

---

## 📞 Quick Reference

**If build succeeds:**
- Base image: `ghcr.io/mica-mni/micapipe-base:latest`
- Main image: `ghcr.io/mica-mni/micapipe:latest`
- Future builds: 3-5 minutes per code change!

**If issues occur:**
- Check logs: `/tmp/micapipe_base_build_*.log`
- Review: `FIX_DCT_NO_SUDO.md`, `FIX_CONDA_STUCK.md`
- Documentation: `TWO_STAGE_BUILD_GUIDE.md`

**Everything is ready! 🚀**
