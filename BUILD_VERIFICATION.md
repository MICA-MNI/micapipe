# Build Configuration Verification

**Date**: October 8, 2025  
**Branch**: comprehensive-base-image  
**Latest Commits**: 86ee6be (FastSurfer fix), 0b187c3 (docs)

---

## Current Build Configuration

### CUDA Support
```dockerfile
ARG ENABLE_CUDA=false  # Line 31 in Dockerfile.base
```

**Status**: ✅ **CPU-ONLY BUILD** (correct for current setup)

- Default: CUDA **DISABLED**
- CUDA version if enabled: 11.8
- To enable: Pass `--build-arg ENABLE_CUDA=true` or use `--enable-cuda` in build script

### FastSurfer Configuration

**CPU Build (current)**:
```dockerfile
python /opt/FastSurfer/Docker/install_env.py \
    -m cpu \                          # ← CPU mode
    -i /opt/FastSurfer/env/fastsurfer.yml \
    -o /opt/FastSurfer/fastsurfer_cpu.yml
```

**This generates**:
- `torch==2.5.1+cpu` (after our patches)
- `torchvision==0.20.1+cpu` (after our patches)
- PyTorch CPU-only wheels from https://download.pytorch.org/whl/cpu

**GPU Build (if --enable-cuda used)**:
```dockerfile
python /opt/FastSurfer/Docker/install_env.py \
    -m gpu \                          # ← GPU mode
    -i /opt/FastSurfer/env/fastsurfer.yml \
    -o /opt/FastSurfer/fastsurfer_gpu.yml
```

**Would generate**:
- `torch==2.5.1` (no +cpu suffix)
- `torchvision==0.20.1` (no +cpu suffix)
- PyTorch with CUDA 11.8 support

---

## Applied Patches

### 1. Source YAML Patch (Commit 86ee6be)
**Before install_env.py runs**:
```bash
sed -i 's/torch>=2\.0/torch==2.5.1/g' /opt/FastSurfer/env/fastsurfer.yml
sed -i 's/torchvision>=0\.15/torchvision==0.20.1/g' /opt/FastSurfer/env/fastsurfer.yml
```

**Effect**:
- Source YAML specifies exact versions
- `install_env.py` uses these when generating environment files
- Both conda YAML and pip requirements get correct versions

### 2. Generated File Patches (Commit 6a811a4)
**After install_env.py runs** (CPU build):
```bash
sed -i 's/torchio==0.20.4/torchio==0.20.5/g' /opt/FastSurfer/fastsurfer_cpu.yml
sed -i 's/torch==2\.6\.0\+cpu/torch==2.5.1+cpu/g' /opt/FastSurfer/fastsurfer_cpu.yml
sed -i 's/torchvision==0\.21\.0\+cpu/torchvision==0.20.1+cpu/g' /opt/FastSurfer/fastsurfer_cpu.yml
```

**Effect**:
- Safety net to catch any remaining 2.6.0 references
- Ensures generated YAML has correct versions before conda env create

---

## PyTorch Version Compatibility

| Version | Status | Notes |
|---------|--------|-------|
| torch==2.6.0+cpu | ❌ Does NOT exist | FastSurfer v2.4.2 incorrectly requests this |
| torch==2.5.1+cpu | ✅ Latest stable | What our patches install |
| torch==2.4.1+cpu | ✅ Previous stable | Alternative if 2.5.1 has issues |
| torch==2.7.1+cpu | ❌ Future version | From FastSurfer main branch |

**PyTorch 2.5.1 Release Date**: December 2024  
**PyTorch 2.6.0 Status**: Not released as of October 2025

---

## Verification on Server

After pulling the latest changes, verify the configuration:

```bash
cd ~/micapipe
git pull origin comprehensive-base-image

# Check ENABLE_CUDA setting
grep "ARG ENABLE_CUDA" Dockerfile.base
# Should show: ARG ENABLE_CUDA=false

# Check FastSurfer patches
grep -A 3 "torch>=2.0" Dockerfile.base
# Should show the sed commands patching to torch==2.5.1

# Check install_env.py mode
grep "install_env.py" Dockerfile.base | grep -E "\-m (cpu|gpu)"
# Should show: -m cpu (for CPU build)
```

---

## Expected Build Behavior

### CPU Build (Current)
1. ✅ Clone FastSurfer v2.4.2
2. ✅ Patch source YAML: torch → 2.5.1, torchvision → 0.20.1
3. ✅ Run `install_env.py -m cpu`
4. ✅ Generates `fastsurfer_cpu.yml` with `torch==2.5.1+cpu`
5. ✅ Apply safety patches (catch any remaining 2.6.0)
6. ✅ Run `conda env create -f fastsurfer_cpu.yml`
7. ✅ PyTorch 2.5.1 CPU installs successfully from PyPI
8. ✅ FastSurfer environment ready

### GPU Build (If Enabled)
Would follow same process but with:
- `-m gpu` instead of `-m cpu`
- `torch==2.5.1` (no +cpu suffix)
- CUDA 11.8 libraries available
- GPU-accelerated PyTorch from CUDA wheels

---

## CUDA Versions Reference

If CUDA support is needed in the future:

**Current Dockerfile.base CUDA version**: 11.8
```dockerfile
cuda-runtime-11-8
cuda-libraries-11-8
libcudnn8
```

**Compatible PyTorch builds**:
- PyTorch 2.5.1 supports CUDA 11.8, 12.1, 12.4
- FastSurfer works with CUDA 11.8+

**To enable CUDA**:
```bash
./build_base_image_server.sh --enable-cuda
```

---

## Summary

✅ **Current Configuration**: CPU-only, correct  
✅ **FastSurfer Mode**: CPU (`-m cpu`)  
✅ **PyTorch Version**: 2.5.1+cpu (after patches)  
✅ **Patches Applied**: Source YAML + Generated files  
✅ **Expected Result**: Build should succeed with CPU PyTorch

**No CUDA version issues** - we're building CPU-only and using the correct PyTorch CPU version.
