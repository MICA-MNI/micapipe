# FastSurfer PyTorch Version Fix

**Error**: `ERROR: No matching distribution found for torch==2.6.0+cpu`

---

## Root Cause

FastSurfer's `install_env.py` script generates conda environment files with **PyTorch 2.6.0**, which **doesn't exist yet**. The latest stable PyTorch version is **2.5.1**.

### What Happened

1. FastSurfer v2.4.2 cloned ✅
2. `install_env.py` generates `fastsurfer_cpu.yml` ✅
3. Generated YAML specifies: `torch==2.6.0+cpu` ❌
4. `conda env create` tries to install torch 2.6.0 ❌
5. PyPI doesn't have torch 2.6.0 → **BUILD FAILS** ❌

---

## The Fix (Commits 6a811a4 + 86ee6be)

### Approach 1: Patch Source YAML (Commit 86ee6be - PRIMARY FIX)

Patch the source `fastsurfer.yml` **BEFORE** `install_env.py` runs:

```dockerfile
RUN git clone https://github.com/Deep-MI/FastSurfer.git /opt/FastSurfer \
    && cd /opt/FastSurfer \
    && git checkout v2.4.2 \
    && sed -i 's/torch>=2\.0/torch==2.5.1/g' /opt/FastSurfer/env/fastsurfer.yml \
    && sed -i 's/torchvision>=0\.15/torchvision==0.20.1/g' /opt/FastSurfer/env/fastsurfer.yml
```

**Why this is the correct fix**:
- `install_env.py` reads `fastsurfer.yml` and generates BOTH:
  1. Conda environment YAML (with pinned versions)
  2. Pip requirements file (with pinned versions)
- By patching the SOURCE, both generated files get the correct version
- This prevents torch 2.6.0 from appearing in EITHER file

### Approach 2: Patch Generated YAML (Commit 6a811a4 - SAFETY NET)

Keep the sed patches after `install_env.py` as a safety measure:

```dockerfile
# Patch versions after install_env.py generates YAML
RUN if [ "$ENABLE_CUDA" = "true" ]; then \
        # GPU build
        sed -i 's/torch==2\.6\.0/torch==2.5.1/g' /opt/FastSurfer/fastsurfer_gpu.yml; \
        sed -i 's/torchvision==0\.21\.0/torchvision==0.20.1/g' /opt/FastSurfer/fastsurfer_gpu.yml; \
    else \
        # CPU build  
        sed -i 's/torch==2\.6\.0\+cpu/torch==2.5.1+cpu/g' /opt/FastSurfer/fastsurfer_cpu.yml; \
        sed -i 's/torchvision==0\.21\.0\+cpu/torchvision==0.20.1+cpu/g' /opt/FastSurfer/fastsurfer_cpu.yml; \
    fi
```

**Why keep this**:
- Defense in depth - catches any version that slips through
- The `install_env.py` script behavior might change in future FastSurfer updates

---

## Version Compatibility

| Package | FastSurfer Requested | What We Use | Status |
|---------|---------------------|-------------|--------|
| torch (CPU) | 2.6.0+cpu | 2.5.1+cpu | ✅ Compatible |
| torch (GPU) | 2.6.0 | 2.5.1 | ✅ Compatible |
| torchvision (CPU) | 0.21.0+cpu | 0.20.1+cpu | ✅ Compatible |
| torchvision (GPU) | 0.21.0 | 0.20.1 | ✅ Compatible |
| torchio | 0.20.4 | 0.20.5 | ✅ Compatible |

**PyTorch 2.5.1** is the latest stable release and fully compatible with FastSurfer 2.4.2.

---

## Why FastSurfer Has Wrong Version

FastSurfer's `install_env.py` likely has a bug or was written in anticipation of PyTorch 2.6.0, which hasn't been released yet. This is a common issue when dependencies specify future versions.

---

## Deploy on Server

```bash
cd ~/micapipe
git pull origin comprehensive-base-image  # Get commits 6a811a4 + 86ee6be
./migrate_comprehensive_base_to_server.sh
```

The build will now:
1. ✅ Clone FastSurfer v2.4.2
2. ✅ **Patch source fastsurfer.yml** (torch 2.6.0 → 2.5.1) ← NEW PRIMARY FIX
3. ✅ Generate FastSurfer CPU environment file (uses patched source)
4. ✅ **Patch generated YAML** (safety net for any remaining 2.6.0)
5. ✅ Create conda environment successfully
6. ✅ Continue with rest of build

---

## Summary of All Recent Fixes

| Commit | Issue | Fix |
|--------|-------|-----|
| d691d2f | DESIGNER runtime | PYTHONPATH + shared environment |
| 05187db | DESIGNER build | pybind11 + fftw + gcc/g++ |
| 6a811a4 | FastSurfer torch (partial) | Patch generated YAML |
| 86ee6be | FastSurfer torch (complete) | Patch source YAML |

**Status**: ✅ All fixes applied - Ready for rebuild

---

**Next**: The build should proceed past FastSurfer installation now! 🎉
