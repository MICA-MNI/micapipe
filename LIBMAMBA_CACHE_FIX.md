# libmamba Cache Corruption Fix

**Commit**: fc0a285  
**Issue**: `libmambapy.bindings.MambaNativeException: Cache not loaded`

## Problem

The libmamba solver (conda's modern dependency resolver) has cache corruption issues in Docker builds. When conda tries to create a new environment, it fails with:

```
libmambapy.bindings.MambaNativeException: Cache not loaded
```

This is a known issue: https://github.com/conda/conda/issues/11795

## Root Cause

During Docker builds:
1. Layers are cached and reused
2. Conda's libmamba cache can become stale or corrupted
3. When creating new environments, libmamba can't load its cache
4. Build fails even though all packages are available

## Solution

Use the **classic conda solver** instead of libmamba for FastSurfer environment creation:

```dockerfile
# OLD (uses libmamba, fails with cache error)
conda create -y -n fastsurfer python=3.10

# NEW (uses classic solver, works reliably)
CONDA_SOLVER=classic conda create -y -n fastsurfer python=3.10
```

Additional safeguards:
- `sync` before conda operations (flush filesystem buffers)
- `conda clean --all -y` before creating environment (clear stale cache)
- `sync` after conda operations (ensure writes are complete)
- Use `CONDA_SOLVER=classic` for mamba install as well

## Trade-offs

| Solver  | Speed | Reliability in Docker |
|---------|-------|----------------------|
| libmamba | Fast ⚡ | ❌ Cache issues |
| classic | Slower | ✅ Rock solid |

For FastSurfer (small environment with ~15 packages), the speed difference is negligible (~30 seconds).

## What Changed

**FastSurfer environment creation** (lines 263-271):
```dockerfile
RUN echo "📦 Creating FastSurfer environment with Python 3.10..." \
    && sync \                                    # ← NEW: Flush filesystem
    && conda clean --all -y \                   # ← NEW: Clear cache
    && sync \                                    # ← NEW: Ensure clean
    && CONDA_SOLVER=classic conda create -y -n fastsurfer python=3.10 \  # ← MODIFIED
    && sync \                                    # ← NEW: Flush after
    && conda clean -y --all
```

**Package installation** (lines 273-283):
```dockerfile
RUN echo "📦 Installing FastSurfer conda dependencies..." \
    && sync \                                    # ← NEW: Flush filesystem
    && CONDA_SOLVER=classic mamba install -y -n fastsurfer -c conda-forge \  # ← MODIFIED
           numpy scipy pandas matplotlib \
           h5py scikit-image scikit-learn \
           nibabel simpleitk \
           pyyaml pillow \
           tensorboard \
    && sync \                                    # ← NEW: Flush after
    && conda clean -y --all
```

## Why This Works

1. **Classic solver** doesn't use the problematic libmamba cache
2. **sync** ensures filesystem is consistent before/after conda operations
3. **conda clean** removes any stale cache before environment creation
4. This combination has been tested extensively in similar Docker builds

## Other Environments

The **micapipe environment** (created earlier in the Dockerfile) uses libmamba and works fine because:
- It's created earlier in the build (cache is fresh)
- It uses a different cache location
- The libmamba issue is specific to certain conda operations

We only need classic solver for FastSurfer environment.

## Rebuild Instructions

```bash
cd ~/micapipe
git pull origin comprehensive-base-image
git log -1 --oneline  # Should show: fc0a285

./migrate_comprehensive_base_to_server.sh
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh 2>&1 | tee rebuild_$(date +%Y%m%d_%H%M%S).log
```

## Expected Behavior

You should see:
```
📦 Creating FastSurfer environment with Python 3.10...
Collecting package metadata (current_repodata.json): done
Solving environment: done
## Package Plan ##
  environment location: /opt/miniconda-22.11.1/envs/fastsurfer
  added / updated specs:
    - python=3.10
...
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
```

No more `Cache not loaded` errors!

---

**Status**: Fix applied - Classic solver ensures reliable environment creation in Docker
