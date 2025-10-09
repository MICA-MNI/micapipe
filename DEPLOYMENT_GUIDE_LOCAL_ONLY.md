# Server Deployment Instructions - LOCAL IMAGES ONLY (No Registry)

This guide explains how to deploy the updated micapipe CI system using only local Docker images (no ghcr.io or external registries required).

## Understanding the Image Setup

- ❌ **ghcr.io**: GitHub Container Registry (we're NOT using this)
- ✅ **Local images**: Built and stored locally on your server
- ✅ **Image names**: `micapipe-base:latest`, `micapipe:latest` (local only)

## Prerequisites

- Access to server: `/data_/mica1/03_projects/actions-runner/`
- Docker service running (managed by admin)
- Singularity installed system-wide
- Your updated micapipe code on the `comprehensive-base-image` branch

## Option 1: Quick Setup (Recommended for Local Images)

### Step 1: Update Action Runner for Local-Only
```bash
# Navigate to your action runner directory
cd /data_/mica1/03_projects/actions-runner/

# Copy the local-only Dockerfile (no registry dependencies)
cp /path/to/micapipe/github-actions-runner-reference/Dockerfile.local-only ./Dockerfile
```

### Step 2: Build Local-Only Action Runner
```bash
# Build the action runner (much smaller, no registry pulls)
export DOCKER_CONTENT_TRUST=0
docker build -t micapipe-runner-local .

# This runner will build everything locally as needed
```

### Step 3: Replace Current Runner
```bash
# Stop existing runner
docker stop micapipe-runner || echo "No existing runner"
docker rm micapipe-runner || echo "No existing runner to remove"

# Start local-only runner
docker run -d --name micapipe-runner --restart always --privileged \
    -v /data/mica1/03_projects/enning/BIDS_CI/sing_cache:/data/mica1/03_projects/enning/BIDS_CI/sing_cache \
    -v /export02/local/singularity_tmp/:/export02/local/singularity_tmp/ \
    -v /data/mica3/BIDS_CI/rawdata:/data/mica3/BIDS_CI/rawdata \
    -v /data_/mica3/BIDS_CI:/data_/mica3/BIDS_CI \
    -v /data/mica1/03_projects/enning/BIDS_CI:/data/mica1/03_projects/enning/BIDS_CI \
    -v /tmp:/tmp \
    -v /var/run/docker.sock:/var/run/docker.sock \
    micapipe-runner-local
```

### Step 4: Update GitHub Workflow for Local Images
```bash
# In your micapipe repository, use the local images workflow
cp .github/workflows/ci_test_local_images.yml .github/workflows/ci_test.yml

# Commit and push
git add .github/workflows/ci_test.yml
git commit -m "Update CI for local-only Docker images"
git push origin comprehensive-base-image
```

## Option 2: Pre-build Base Image for Maximum Speed

### Step 1: Build Base Image Once (45-60 minutes)
```bash
# Navigate to your micapipe source
cd /path/to/micapipe/

# Build the base image locally (one-time setup)
docker build -f Dockerfile.base -t micapipe-base:latest .

# Verify the base image
docker images | grep micapipe-base
```

### Step 2: Test Build Speed
```bash
# Test the main image build (should be ~3-5 minutes now)
docker build -f Dockerfile.main --build-arg "BASE_IMAGE=micapipe-base:latest" -t micapipe:test .

# Time this - should be very fast compared to full build
```

### Step 3: Follow Option 1 Steps 1-4
The CI will now find your pre-built base image and use it for fast builds.

## How the Local CI Works

### First Run (No Base Image)
```bash
# CI detects no base image exists
echo "❌ No base image found - will build new base"

# Builds base image (45-60 minutes)
docker build -f Dockerfile.base -t micapipe-base:latest .

# Builds main image (3-5 minutes)
docker build -f Dockerfile.main --build-arg "BASE_IMAGE=micapipe-base:latest" -t micapipe:latest .
```

### Subsequent Runs (Base Image Cached)
```bash
# CI detects existing base image
echo "✅ Found local base image: micapipe-base:latest"

# Skips base build, only builds main image (3-5 minutes)
docker build -f Dockerfile.main --build-arg "BASE_IMAGE=micapipe-base:latest" -t micapipe:latest .
```

### Legacy Fallback (Single Dockerfile)
```bash
# If no Dockerfile.base exists, falls back to v1 behavior
docker build -t micapipe:test .
```

## Image Lifecycle Management

### Local Images Created
```bash
# After first run, you'll have:
micapipe-base:latest    # Heavy tools (FSL, FreeSurfer, etc.) - KEEP THIS
micapipe:latest         # Micapipe code on top of base - Rebuilt each time

# Check images
docker images | grep micapipe
```

### Cleanup Strategy
```bash
# The CI automatically:
# ✅ Keeps base image (for speed)
# ✅ Removes main image after use
# ✅ Cleans up dangling images

# Manual cleanup if needed:
docker image prune -f
```

## Troubleshooting Local Setup

### Issue: "Image not found" in CI
```bash
# Check what images exist locally
docker images | grep micapipe

# If none, the CI will build them from scratch (expected on first run)
```

### Issue: Split-stage build fails
```bash
# CI will automatically fall back to single-stage
# Check if you have both Dockerfile.base and Dockerfile.main
ls -la Dockerfile*
```

### Issue: Base image build takes too long
```bash
# This is normal for first build (45-60 minutes)
# Subsequent builds reuse this image and are much faster

# Monitor build progress:
docker logs -f $(docker ps -q --filter ancestor=micapipe-runner-local)
```

## Expected Performance (Local Images)

| Build Stage | First Run | Subsequent Runs | 
|-------------|-----------|-----------------|
| Base Image Build | 45-60 min | 0 min (cached) |
| Main Image Build | 3-5 min | 3-5 min |
| Singularity Convert | 5-10 min | 5-10 min |
| **Total Time** | **55-75 min** | **10-15 min** |

## Verification Commands

### Check Action Runner
```bash
docker logs micapipe-runner

# Should show connection to GitHub
```

### Check Local Images
```bash
docker images | grep micapipe

# Should show:
# micapipe-base    latest    ...    (large size, ~8-12GB)
# micapipe         latest    ...    (smaller, recent)
```

### Monitor CI Build
```bash
# Create test commit
echo "# Test local build" >> README.md
git add README.md
git commit -m "Test local CI build"
git push origin comprehensive-base-image

# Watch in GitHub Actions - first run ~60min, subsequent ~15min
```

## Key Differences from Registry Approach

| Aspect | Registry (ghcr.io) | Local Only |
|--------|-------------------|------------|
| Base Image Storage | GitHub Container Registry | Local Docker |
| Internet Dependency | Required for pulls | None after setup |
| Action Runner Size | Larger (pre-pulls) | Smaller |
| Build Speed | Instant (if published) | Fast (if cached locally) |
| Setup Complexity | Higher (needs publishing) | Lower (local only) |
| **Recommended For** | Production/Teams | Development/Single Server |

## Next Steps

1. ✅ Use the local-only action runner (`Dockerfile.local-only`)
2. ✅ Use the local images CI workflow (`ci_test_local_images.yml`)  
3. ✅ Test with a commit (first run will be slow, builds base image)
4. ✅ Verify subsequent runs are much faster (~15 minutes vs 60+ minutes)

This approach keeps everything local and doesn't depend on any external registries like ghcr.io.