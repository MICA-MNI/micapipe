# 🔍 ERROR PREVENTION & SERVER INSTRUCTIONS

## 🚨 CRITICAL ERRORS TO PREVENT

Based on analysis of docker-container-updates branch, here are the key errors to prevent:

### 1. Missing Pre-Downloaded Files Error
**Problem**: Docker build fails because pre-downloaded files not found
**Prevention**: Verify files exist before building
```bash
# Check on server:
ls -la /host/cassio/export03/data/enning/downloads/*.tar.gz
ls -la /host/cassio/export03/data/enning/downloads/*.tgz
```

### 2. Wrong Directory Error  
**Problem**: Building from wrong location
**Prevention**: Always build from downloads directory where files exist
```bash
# Must be in this exact directory:
cd /host/cassio/export03/data/enning/downloads
pwd  # Should show: /host/cassio/export03/data/enning/downloads
```

### 3. Dockerfile Not Found Error
**Problem**: Dockerfile.mamba-base doesn't exist in build directory  
**Prevention**: Ensure git pull brought all files
```bash
# Check files exist:
ls -la Dockerfile.mamba-base build_comprehensive_base_server.sh
```

### 4. Docker Daemon Access Error
**Problem**: Cannot connect to Docker daemon
**Prevention**: Check Docker access before building
```bash
# Test Docker access:
docker info
docker images
```

### 5. Memory/Space Errors
**Problem**: Insufficient memory or disk space during build
**Prevention**: Check resources before building
```bash
# Check memory:
free -h
# Check disk space:
df -h /host/cassio/export03/data/enning
```

## ✅ CORRECT SERVER INSTRUCTIONS

### Step 1: Navigate to Server Build Directory
```bash
cd /host/cassio/export03/data/enning/downloads
```

### Step 2: Update Code (Git Pull)
```bash
# Make sure you're on the right branch:
git branch --show-current
# Should show: comprehensive-base-image

# Pull latest changes:
git pull origin comprehensive-base-image
```

### Step 3: Verify Pre-Downloaded Files
```bash
# Check required files exist:
echo "🔍 Checking pre-downloaded files..."
ls -la freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz 2>/dev/null && echo "✅ FreeSurfer found" || echo "❌ FreeSurfer missing"
ls -la fsl-6.0.2-centos6_64.tar.gz 2>/dev/null && echo "✅ FSL found" || echo "❌ FSL missing"  
ls -la afni-linux_openmp_64.tgz 2>/dev/null && echo "✅ AFNI found" || echo "❌ AFNI missing"
```

### Step 4: Verify Docker Access
```bash
# Test Docker daemon:
docker info >/dev/null 2>&1 && echo "✅ Docker accessible" || echo "❌ Docker not accessible"
```

### Step 5: Build Comprehensive Base Image
```bash
# Build base image (one-time, ~45-77 min):
docker build -f Dockerfile.mamba-base -t micapipe-comprehensive-base .

# OR use the provided script:
./build_comprehensive_base_server.sh
```

### Step 6: Build Main Image (Fast Daily Builds)
```bash
# Build main image (~3-5 min):
docker build -f Dockerfile.minimal -t micapipe:latest .

# OR use fast CI script:
./build_fast_ci_server.sh
```

## 🎯 EXPECTED RESULTS

### Successful Base Build Output:
```
✅ Copying FreeSurfer...
✅ Copying FSL... 
✅ Copying AFNI...
✅ Using pre-downloaded Miniconda installer
```

### Failed Build Indicators:
```
❌ FSL not found, will download during build
❌ FreeSurfer not found, will download during build
❌ Cannot connect to Docker daemon
```

## 🚀 PERFORMANCE EXPECTATIONS

- **Base Image Build**: 45-77 minutes (one-time)
- **Main Image Build**: 3-5 minutes (daily)
- **Total Size**: ~8-12GB (comprehensive neuroimaging environment)
- **Cache Hit Rate**: 80-95% for subsequent builds

## 🔧 TROUBLESHOOTING COMMANDS

### If Pre-Downloaded Files Missing:
```bash
# Re-download dependencies:
./download_dependencies.sh
```

### If Docker Permission Denied:
```bash
# Add user to docker group (requires admin):
sudo usermod -aG docker $USER
# Then logout and login again
```

### If Build Runs Out of Space:
```bash
# Clean Docker cache:
docker system prune -f
docker builder prune -f
```

### If Build Fails with Memory Error:
```bash
# Build with memory limits:
docker build --memory=8g --memory-swap=12g -f Dockerfile.mamba-base -t micapipe-base .
```

---

**🎯 SUCCESS CRITERIA**: Base build completes without downloading FSL/FreeSurfer (uses pre-downloaded files), main image builds in under 5 minutes.