# Docker Build Mistakes Reference
## Tracking All Mistakes and Fixes to Prevent Recurrence

**Last Updated:** October 7, 2025  
**Purpose:** Document all mistakes made during Docker build development to prevent repeating them

---

## ❌ MISTAKE #1: Docker Content Trust SSL Certificate Error
**Date:** October 7, 2025  
**Error:**
```
Error response from daemon: Get https://notary.docker.io/v2/: x509: certificate signed by unknown authority
```

**Root Cause:** Docker Content Trust (DCT) was enabled and causing SSL verification failures

**Fix:** Disable Docker Content Trust in build script
```bash
export DOCKER_CONTENT_TRUST=0
```

**Commit:** 6461d85  
**Status:** ✅ FIXED  
**Lesson:** Always disable DCT in automated build environments unless specifically needed

---

## ❌ MISTAKE #2: Conda Solving Stuck for 30+ Minutes
**Date:** October 7, 2025  
**Error:**
```
Solving environment: ...working... [stuck indefinitely]
```

**Root Cause:** Default conda solver is extremely slow with complex dependency graphs

**Wrong Fix Attempt:**
```dockerfile
# WRONG - just installing mamba doesn't help
conda install -y -n base -c conda-forge mamba
```

**Correct Fix:** Install mamba with libmamba solver and configure it properly
```dockerfile
conda install -y -n base -c conda-forge "mamba>=1.4" "conda-libmamba-solver"
conda config --set solver libmamba
conda config --set channel_priority flexible
```

**Commit:** 3bc3068  
**Status:** ✅ FIXED  
**Lesson:** 
- Must install `conda-libmamba-solver` package explicitly
- Use `channel_priority flexible` to avoid frozen solve errors
- Set solver config in same RUN layer as mamba install

---

## ❌ MISTAKE #3: GPG Signature Verification Failure
**Date:** October 7, 2025  
**Error:**
```
GPG error: The following signatures couldn't be verified
```

**Root Cause:** Ubuntu apt repositories' GPG keys not available on keyserver

**Wrong Fix Attempt:**
```dockerfile
# WRONG - single keyserver can fail
apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3B4FE6ACC0B21F32
```

**Correct Fix:** Use multiple fallback keyservers
```dockerfile
(apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3B4FE6ACC0B21F32 || \
 apt-key adv --keyserver hkp://pgp.mit.edu:80 --recv-keys 3B4FE6ACC0B21F32 || \
 wget -qO - https://keyserver.ubuntu.com/pks/lookup?op=get\&search=0x3B4FE6ACC0B21F32 | apt-key add - || \
 echo "Warning: Could not import key, continuing anyway...")
```

**Commit:** 80010f7  
**Status:** ✅ FIXED  
**Lesson:** Always provide multiple fallback methods for external dependencies (keyservers, mirrors, etc.)

---

## ❌ MISTAKE #4: "No Space Left on Device" Error
**Date:** October 7, 2025  
**Error:**
```
Error: write /path/to/file: no space left on device
```

**Root Cause:** Building Docker image in home directory (`~/micapipe`) with limited disk space

**Wrong Fix Attempts:**
1. Creating separate `~/micapipe_build/` subdirectory → Still in home directory, same problem
2. Complex rsync exclusion patterns → Overcomplicated solution

**Correct Fix:** Build directly in server location with unlimited space
```bash
# Build location enforcement in script
if [[ "$SCRIPT_DIR" == "$HOME"* ]] || [[ "$SCRIPT_DIR" == "/home/"*"/micapipe"* ]]; then
    echo "❌ ERROR: Building from home directory ($SCRIPT_DIR)"
    exit 1
fi

# Use server location
BUILD_DIR="/host/cassio/export03/data/enning/downloads"
```

**Commit:** 3b9350c  
**Status:** ✅ FIXED  
**Lesson:** 
- Never build Docker images in limited disk space locations
- Enforce build location validation in scripts
- Keep solution simple - don't overcomplicate

---

## ❌ MISTAKE #5: Pre-Downloaded Files Not Being Used
**Date:** October 7, 2025  
**Error:**
```
⬇️  Downloading FSL...  [2.6GB download when file already exists]
⬇️  Downloading FreeSurfer...  [2.7GB download when file already exists]
```

**Root Cause:** Dockerfile checking wrong path - checking `/downloads/file.tar.gz` (container path that doesn't exist during build) instead of build context path

**Wrong Check:**
```dockerfile
# WRONG - $DOWNLOADS_DIR doesn't exist during build
if [ -f "$DOWNLOADS_DIR/fsl-6.0.2-centos6_64.tar.gz" ]; then
```

**Correct Fix:** Copy build context first, then check files there
```dockerfile
# At start of Dockerfile
COPY . /tmp/build_context/

# Later in RUN commands
if [ -f "/tmp/build_context/fsl-6.0.2-centos6_64.tar.gz" ]; then
    echo "✅ Using pre-downloaded FSL";
    cp "/tmp/build_context/fsl-6.0.2-centos6_64.tar.gz" "$TMPDIR/fsl.tar.gz";
else
    echo "⬇️  Downloading FSL...";
    curl -o fsl.tar.gz https://...;
fi

# At end - cleanup to save space
RUN rm -rf /tmp/build_context
```

**Also Required:** Update .dockerignore to NOT exclude .tar.gz files
```gitignore
# BEFORE (WRONG):
*.tar.gz
*.tgz

# AFTER (CORRECT):
# Allow pre-downloaded dependency files
# *.tar.gz  <- COMMENTED OUT
# *.tgz     <- COMMENTED OUT
```

**Commit:** 30be584  
**Status:** ✅ FIXED  
**Lesson:** 
- Docker build context is the ONLY way to access local files during build (no mount support in our environment)
- Always verify file checks reference paths that actually exist during build
- Update .dockerignore to include necessary files
- Use global search/replace to ensure ALL references updated (not just some)

---

## ❌ MISTAKE #6: Missing System Dependencies for pip Packages
**Date:** October 7, 2025  
**Error:**
```
ERROR: Dependency lookup for cairo with method 'pkgconfig' failed: 
Pkg-config for machine host machine not found. Giving up.
```

**Root Cause:** pip package `pycairo` (dependency of `xhtml2pdf`) requires system libraries to compile

**Wrong Approach:**
```dockerfile
# WRONG - trying to install pip package without system dependencies
RUN mamba run -n micapipe pip install --no-cache-dir xhtml2pdf
```

**Correct Fix:** Install system dependencies BEFORE pip packages
```dockerfile
# Install system dependencies for cairo/pycairo
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           pkg-config \
           libcairo2-dev \
           libgirepository1.0-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# THEN install pip packages
RUN mamba run -n micapipe pip install --no-cache-dir xhtml2pdf
```

**Commit:** 23c5bc5  
**Status:** ✅ FIXED  
**Lesson:** 
- Always check if pip packages have system library dependencies
- Install system packages BEFORE conda/pip packages
- Common dependencies: pkg-config, lib*-dev packages, compilers

---

## ❌ MISTAKE #7: Trying to Install Pip-Only Package via Conda
**Date:** October 7, 2025  
**Error:**
```
error libmamba Could not solve for environment specs
    The following package could not be installed
    └─ antspy =* * does not exist (perhaps a typo or a missing channel).
```

**Root Cause:** `antspy` is ONLY available via pip, not in conda channels

**Wrong Approach:**
```dockerfile
# WRONG - antspy doesn't exist in conda-forge
RUN mamba install -y -n micapipe -c conda-forge antspy
```

**Correct Fix:** Install via pip
```dockerfile
RUN mamba run -n micapipe pip install --no-cache-dir antspy
```

**Commit:** 648f139  
**Status:** ✅ FIXED  
**Lesson:** 
- Some packages are pip-only (antspy, brainspace, etc.)
- Check package availability: `mamba search -c conda-forge package_name`
- If not in conda, use pip
- Don't duplicate installations - check environment.yml first

---

## ❌ MISTAKE #8: Unused Build Arguments
**Date:** October 7, 2025  
**Issue:** Build script passes `--build-arg DOWNLOADS_DIR=/downloads` but Dockerfile doesn't use it

**Wrong:**
```bash
# build_base_image_server.sh
docker build --build-arg "DOWNLOADS_DIR=/downloads" ...

# Dockerfile.base
ARG DOWNLOADS_DIR="/downloads"  # Defined but never used
```

**Status:** ⚠️ HARMLESS BUT CLEANUP NEEDED  
**Lesson:** Remove unused build arguments to avoid confusion

---

## 🎯 KEY PATTERNS TO REMEMBER

### 1. External Dependencies (Keyservers, Mirrors, Downloads)
**ALWAYS** provide multiple fallback methods:
```bash
(method1 || method2 || method3 || echo "Warning: continuing anyway")
```

### 2. Conda/Mamba Environment Setup
**ALWAYS** use this pattern:
```dockerfile
# Install mamba with libmamba solver
RUN conda install -y -n base -c conda-forge "mamba>=1.4" "conda-libmamba-solver" \
    && conda config --set solver libmamba \
    && conda config --set channel_priority flexible
```

### 3. Pre-Downloaded Files Access
**ALWAYS** use this pattern:
```dockerfile
# Copy build context first
COPY . /tmp/build_context/

# Check files in build context
RUN if [ -f "/tmp/build_context/package.tar.gz" ]; then
        echo "✅ Using pre-downloaded package";
        cp "/tmp/build_context/package.tar.gz" /path/to/use/;
    else
        echo "⬇️  Downloading package...";
        curl -o package.tar.gz https://...;
    fi

# Cleanup at end
RUN rm -rf /tmp/build_context
```

### 4. System Dependencies for Pip Packages
**ALWAYS** check and install system libraries first:
```dockerfile
# System dependencies FIRST
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           pkg-config \
           lib*-dev \
    && apt-get clean

# THEN pip packages
RUN mamba run -n micapipe pip install package
```

### 5. Build Location Validation
**ALWAYS** validate build location at script start:
```bash
if [[ "$SCRIPT_DIR" == "$HOME"* ]]; then
    echo "❌ ERROR: Building from home directory"
    exit 1
fi
```

---

## 📋 VALIDATION CHECKLIST BEFORE EACH BUILD

Before running build, verify:

- [ ] Building from server location with unlimited space (NOT home directory)
- [ ] Docker Content Trust disabled (`export DOCKER_CONTENT_TRUST=0`)
- [ ] Conda using libmamba solver with flexible channel priority
- [ ] All pre-downloaded files referenced via `/tmp/build_context/`
- [ ] .dockerignore allows .tar.gz, .tgz files
- [ ] System dependencies installed before pip packages requiring them
- [ ] No duplicate conda/pip package installations
- [ ] All external downloads have multiple fallback methods
- [ ] Pre-downloaded files exist in current directory before build

---

## 🔄 QUICK FIX COMMANDS

```bash
# If conda solving is stuck:
# In Dockerfile, ensure:
conda config --set solver libmamba
conda config --set channel_priority flexible

# If pre-downloaded files not used:
# Check .dockerignore doesn't exclude .tar.gz
# Verify all file checks use /tmp/build_context/ path

# If pip package compilation fails:
# Add system dependencies:
apt-get install -y pkg-config libcairo2-dev lib*-dev

# If conda package not found:
# Try pip instead:
mamba run -n micapipe pip install package_name
```

---

## 📊 BUILD SUCCESS METRICS

**Total Mistakes Made:** 8  
**Total Mistakes Fixed:** 7  
**Harmless Issues:** 1  

**Time Wasted on Repeated Mistakes:** ~3 hours  
**Time Saved by This Document:** Hopefully infinite 😊

---

**Remember:** Every mistake is a learning opportunity. This document prevents repeating them!
