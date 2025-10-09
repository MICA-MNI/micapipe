# Server Deployment Instructions (No Sudo Required)

This guide explains how to deploy the updated micapipe CI system on the server without sudo rights, compatible with the new comprehensive base image architecture.

## Prerequisites

- Access to server: `/data_/mica1/03_projects/actions-runner/`
- Docker service running (managed by admin)
- Singularity installed system-wide
- Your updated micapipe code on the `comprehensive-base-image` branch

## Option 1: Recommended - Use Pre-built Base Image (Fastest)

### Step 1: Prepare the Action Runner
```bash
# Navigate to your action runner directory
cd /data_/mica1/03_projects/actions-runner/

# Copy the optimized Dockerfile from the reference files
cp /path/to/micapipe/github-actions-runner-reference/Dockerfile.action-runner-optimized ./Dockerfile

# Copy the original run script
cp /path/to/micapipe/github-actions-runner-reference/run_docker.sh ./
```

### Step 2: Build Updated Action Runner (Pre-pulls Base Image)
```bash
# Build the updated action runner that pre-pulls the base image
export DOCKER_CONTENT_TRUST=0
docker build -t micapipe-runner-optimized .

# This will:
# 1. Install all dependencies
# 2. Pre-pull ghcr.io/mica-mni/micapipe-base:latest
# 3. Setup the runner for instant CI builds
```

### Step 3: Stop Current Runner and Start New One
```bash
# Stop existing runner
docker stop micapipe-runner || echo "No existing runner"
docker rm micapipe-runner || echo "No existing runner to remove"

# Start optimized runner
docker run -d --name micapipe-runner --restart always --privileged \
    -v /data/mica1/03_projects/enning/BIDS_CI/sing_cache:/data/mica1/03_projects/enning/BIDS_CI/sing_cache \
    -v /export02/local/singularity_tmp/:/export02/local/singularity_tmp/ \
    -v /data/mica3/BIDS_CI/rawdata:/data/mica3/BIDS_CI/rawdata \
    -v /data_/mica3/BIDS_CI:/data_/mica3/BIDS_CI \
    -v /data/mica1/03_projects/enning/BIDS_CI:/data/mica1/03_projects/enning/BIDS_CI \
    -v /tmp:/tmp \
    -v /var/run/docker.sock:/var/run/docker.sock \
    micapipe-runner-optimized
```

### Step 4: Update GitHub Workflow
```bash
# In your micapipe repository, replace .github/workflows/ci_test.yml with the new version
cp .github/workflows/ci_test_action_runner_compatible.yml .github/workflows/ci_test.yml

# Commit and push
git add .github/workflows/ci_test.yml
git commit -m "Update CI for comprehensive base image compatibility"
git push origin comprehensive-base-image
```

## Option 2: Manual Base Image Build (If Registry Not Available)

### Step 1: Build Base Image Manually
```bash
# Navigate to your micapipe source
cd /path/to/micapipe/

# Build the base image (45-60 minutes, one-time)
docker build -f Dockerfile.base -t ghcr.io/mica-mni/micapipe-base:latest .

# Verify the base image
docker images | grep micapipe-base
```

### Step 2: Follow Option 1 Steps 2-4
The action runner will now find the locally built base image and use it for fast CI builds.

## Option 3: Fallback - Single-Stage Build (If Split Build Fails)

### Step 1: Use Original Action Runner Setup
```bash
# Use the original Dockerfile from reference
cp /path/to/micapipe/github-actions-runner-reference/Dockerfile ./

# The CI will automatically fall back to single-stage builds
```

## Verification Steps

### 1. Check Action Runner Status
```bash
docker logs micapipe-runner

# Should show:
# "✅ Connected to GitHub"
# "🔍 Listening for Jobs..."
```

### 2. Test CI Manually
```bash
# Create a test commit on comprehensive-base-image branch
cd /path/to/micapipe/
echo "# Test" >> README.md
git add README.md
git commit -m "Test CI build"
git push origin comprehensive-base-image

# Monitor in GitHub Actions tab - should complete in ~10-15 minutes (vs 60+ minutes before)
```

### 3. Check Build Speed
- **First run**: 45-60 minutes (builds base image)
- **Subsequent runs**: 10-15 minutes (uses cached base)
- **Main image build**: 3-5 minutes
- **Singularity conversion**: 5-10 minutes

## Troubleshooting

### Issue: Base Image Not Found
```bash
# Check if base image exists
docker images | grep micapipe-base

# If not, manually build it:
docker build -f Dockerfile.base -t ghcr.io/mica-mni/micapipe-base:latest .
```

### Issue: Docker Permission Denied
```bash
# Check if your user is in docker group (ask admin)
groups $USER

# If not in docker group, ask admin to add you:
# sudo usermod -aG docker $USER
```

### Issue: Singularity Build Fails
```bash
# Check singularity version
singularity --version

# Should be 3.10.2 or compatible
# Check available space in /export02/local/singularity_tmp/
df -h /export02/local/singularity_tmp/
```

### Issue: Old CI Still Running
```bash
# The new CI workflow supports both architectures
# It will automatically detect which build method to use
# Old CI will continue to work until you update the workflow file
```

## Expected Performance Improvements

| Build Stage | v1 Time | New Time | Improvement |
|-------------|---------|----------|-------------|
| Base Image (first time) | N/A | 45-60 min | One-time setup |
| Base Image (cached) | 45-60 min | 0 min | ✅ 100% faster |
| Main Image | N/A | 3-5 min | New stage |
| Full Build | 45-60 min | 10-15 min | ✅ 75% faster |

## Next Steps

1. **Test the new CI** on a feature branch first
2. **Monitor the first build** (will take 45-60 minutes to build base)
3. **Verify subsequent builds** are much faster
4. **Update any other CI workflows** to use the same pattern

## Files Created/Modified

- ✅ `.github/workflows/ci_test_action_runner_compatible.yml` - New CI workflow
- ✅ `github-actions-runner-reference/Dockerfile.action-runner-optimized` - Updated action runner
- ✅ This deployment guide

All files are backward compatible - old CI will continue working until you switch workflows.