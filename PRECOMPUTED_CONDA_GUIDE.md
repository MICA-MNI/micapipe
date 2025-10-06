# HOW TO USE PRE-COMPUTED CONDA ENVIRONMENT IN DOCKERFILE

## Step 1: Generate the pre-computed environment (ONE TIME)

Run this on a machine with conda/mamba installed:

```bash
cd ~/micapipe
./precompute_conda_environment.sh
```

This will create:
- `micapipe_environment_explicit.txt` - Pre-solved package list with URLs

## Step 2: Replace the conda installation in Dockerfile.base

Find this section in Dockerfile.base:
```dockerfile
# ============================================================================
# MICAPIPE CONDA ENVIRONMENT - OPTIMIZED INSTALLATION
# ============================================================================

# Install ALL packages in ONE command to avoid dependency conflicts
RUN echo "📦 Installing micapipe conda environment (optimized single-pass)..." \
    && mamba install -y -n micapipe -c conda-forge \
           "numpy==1.21.5" \
           ... lots of packages ...
    && mamba clean -y --all
```

Replace it with:
```dockerfile
# ============================================================================
# MICAPIPE CONDA ENVIRONMENT - PRE-COMPUTED (INSTANT!)
# ============================================================================

# Copy pre-computed environment specification
COPY micapipe_environment_explicit.txt /tmp/

# Install from explicit spec (NO SOLVING - instant!)
RUN echo "📦 Installing micapipe conda environment (pre-computed)..." \
    && conda create -n micapipe --file /tmp/micapipe_environment_explicit.txt \
    && mamba clean -y --all \
    && rm /tmp/micapipe_environment_explicit.txt \
    && echo "✅ Pre-computed environment installed instantly!"

# Install pip packages that need to be fresh
RUN echo "📦 Installing fresh pip packages..." \
    && mamba run -n micapipe pip install --no-cache-dir \
           git+https://github.com/MICA-MNI/LAMAReg.git \
           git+https://github.com/MICA-MNI/ENIGMA.git \
    && rm -rf ~/.cache/pip/* \
    && echo "✅ Pip packages installed"
```

## Step 3: Update build process

The pre-computed file needs to be available when building:

```bash
# On server
cd /host/cassio/export03/data/enning/downloads

# Make sure micapipe_environment_explicit.txt is present
ls -lh micapipe_environment_explicit.txt

# Build (will be MUCH faster!)
./build_base_image_server.sh
```

## Performance Comparison

| Method | Solving Time | Install Time | Total |
|--------|--------------|--------------|-------|
| **Current (multiple installs)** | 30+ min (stuck) | 10 min | 40+ min |
| **Optimized (single install)** | 5-8 min | 10 min | 15-18 min |
| **Pre-computed (explicit spec)** | 0 sec ⚡ | 3-5 min | 3-5 min |

## Benefits

1. ⚡ **No solving** - Dependencies already computed
2. ⚡ **Deterministic** - Exact same packages every time
3. ⚡ **Faster** - Direct download and install
4. ⚡ **Reproducible** - Lock file ensures consistency
5. ⚡ **No stuck builds** - No dependency solving means no hanging

## When to Regenerate

Regenerate `micapipe_environment_explicit.txt` when:
- You add new packages to `micapipe_environment.yml`
- You update package versions
- You want to update to latest compatible versions

## Alternative: Use conda-lock

For even better reproducibility across platforms:

```bash
# Install conda-lock
pip install conda-lock

# Generate lock file for linux-64
conda-lock -f micapipe_environment.yml -p linux-64

# This creates: conda-linux-64.lock
# Use in Dockerfile:
# RUN mamba create -n micapipe --file conda-linux-64.lock
```

## Quick Start

```bash
# 1. Generate pre-computed environment
./precompute_conda_environment.sh

# 2. Commit the explicit spec file
git add micapipe_environment_explicit.txt
git commit -m "Add pre-computed conda environment for instant builds"

# 3. Update Dockerfile.base (see Step 2 above)

# 4. Build will now skip dependency solving!
```
