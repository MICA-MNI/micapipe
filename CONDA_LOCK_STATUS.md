# NOTE: Conda-lock Pre-computation Status

## Current Status: ⏳ Not Yet Implemented

**Why:** Generating the conda-lock file takes 10-30 minutes even on a local machine due to:
- Complex dependency graph (30+ packages)
- MRtrix3 from custom channel
- ANTsPy has complex C++ dependencies
- Conda solver is inherently slow for this combination

## Current Solution: Optimized Single-Pass Install ✅

The Dockerfile.base currently uses an **optimized single-pass installation** which is much better than the original multi-pass approach:

```dockerfile
# Set libmamba solver for speed
RUN conda config --set solver libmamba

# Single-pass install (much faster than multiple mamba install commands)
RUN mamba install -y -n micapipe -c conda-forge \
       numpy scipy pandas matplotlib ... (all packages)
```

**Performance:**
- Original (multiple installs): 60+ min or stuck
- Current (single-pass): 15-20 min ✅
- Theoretical (with lock file): 2-3 min

## Future Implementation Plan

When time allows, generate lock file on server:

```bash
# On server with more resources
cd /host/cassio/export03/data/enning/downloads

# Install conda-lock
pip install conda-lock

# Generate lock file (may take 15-30 min)
conda-lock lock \
  --file micapipe_environment.yml \
  --platform linux-64 \
  --kind explicit \
  --filename-template "micapipe-{platform}.lock"

# Commit the generated lock file
git add micapipe-linux-64.lock
git commit -m "Add pre-computed conda lock file"
git push
```

Then update Dockerfile to use it.

## Alternative: Accept Current Performance

The current **15-20 minute** conda install time is acceptable for:
- Base image built rarely (monthly)
- Main image builds are fast (3-5 min)
- 95% overall time savings vs original approach

**Recommendation:** Keep current optimized approach, revisit lock file later if needed.

## Files Included

- `micapipe_environment.yml` - Full environment specification (for reference)
- `micapipe_environment_base.yml` - Base packages only (for lock file generation)
- `generate_conda_lock.sh` - Script to generate lock file
- This README documenting status

These files are ready for future use when lock file generation is feasible.
