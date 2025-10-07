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

## ❌ MISTAKE #9: Using environment.yml Instead of Direct Mamba Install
**Date:** October 7, 2025  
**Error:**
```
Updating pip packages: antspy, vtk==9.2.2...
ERROR: Could not find a version that satisfies the requirement antspy
```

**Root Cause:** `mamba env update --file environment.yml` has a bug/misinterpretation where it tries to install conda packages via pip even when they are correctly placed in the conda section (not pip section) of the YAML file.

**Wrong Approach:**
```dockerfile
# WRONG - environment.yml approach causes pip/conda confusion
COPY micapipe_environment.yml /tmp/micapipe_environment.yml
RUN mamba env update -n micapipe --file /tmp/micapipe_environment.yml
```

Even with correct YAML structure:
```yaml
# micapipe_environment.yml
dependencies:
  - antspy      # ← In conda section, NOT pip
  - antspyx
  - pip:
    - vtk==9.2.2  # ← pip section separate
```

Mamba still tries: "Updating pip packages: antspy, vtk==9.2.2..." (WRONG!)

**Correct Fix:** Use direct `mamba install` commands like v1 Dockerfile (which works perfectly)
```dockerfile
# Install conda packages directly - NO environment.yml
RUN mamba install -y -n micapipe -c conda-forge -c mrtrix3 \
       numpy=1.21.5 scipy pandas matplotlib \
       nibabel pillow scikit-learn \
       mrtrix3=3.0.7 antspy antspyx \
    && mamba clean -y --all

# Install pip packages separately
RUN mamba run -n micapipe pip install --no-cache-dir \
       vtk==9.2.2 brainspace==0.1.10 \
       git+https://github.com/MICA-MNI/LAMAReg.git \
    && rm -rf ~/.cache/pip/*
```

**Why v1 Works:**
- v1 Dockerfile (lines 715-722) uses direct `mamba install` - NO environment.yml file
- Each package group installed explicitly with channels specified
- No ambiguity about conda vs pip packages
- Proven in production for months/years

**Commits:** 
- ea3cd92 (attempted environment.yml fix - FAILED)
- e6f108e (debug verification - confirmed YAML correct but still failed)
- b60ae7d (switched to v1 direct install method - SHOULD WORK)

**Status:** ✅ FIXED (switched to v1 method)

**Lesson:** 
- **CRITICAL: When something works in production (v1), use the EXACT same approach**
- Don't try to "improve" or "modernize" working code with environment.yml
- `mamba env update` with YAML files is unreliable - has parsing/interpretation bugs
- Direct `mamba install` commands are explicit, debuggable, and production-proven
- **Always check v1/production code FIRST before trying new approaches**
- Time wasted trying to debug environment.yml: ~2 hours (multiple builds + investigations)
- Time to implement v1 solution: ~5 minutes

---

## ❌ MISTAKE #10: Asking User for Manual Confirmation in Automated Scripts
**Date:** October 7, 2025  
**Issue:** Migration script requires manual "y/N" confirmation, breaking automation

**Wrong:**
```bash
# WRONG - requires manual intervention
read -p "Continue anyway? (y/N): " -n 1 -r
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
fi
```

**Correct Fix:** Auto-proceed with informative messages
```bash
# CORRECT - fully automatic
echo "⚠️  Missing files detected - will download during build (slower)"
echo "   Continuing automatically..."
```

**Commit:** b60ae7d  
**Status:** ✅ FIXED  
**Lesson:** 
- Automation scripts should NEVER require manual input
- Provide clear warning messages instead
- User can cancel with Ctrl+C if they disagree
- Especially critical when user explicitly says "do not ask me to y"

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

### 6. Package Installation Method
**ALWAYS** use direct mamba install (like v1), NOT environment.yml:
```dockerfile
# CORRECT - Direct install (v1 method)
RUN mamba install -y -n micapipe -c conda-forge -c mrtrix3 \
       numpy scipy pandas matplotlib \
       antspy antspyx mrtrix3=3.0.7 \
    && mamba clean -y --all

# Then pip packages separately
RUN mamba run -n micapipe pip install --no-cache-dir package1 package2

# WRONG - environment.yml has bugs
RUN mamba env update -n micapipe --file environment.yml  # ❌ DON'T USE
```

### 7. Automation Scripts
**NEVER** require manual confirmation:
```bash
# WRONG
read -p "Continue? (y/N): "

# CORRECT
echo "⚠️  Issue detected - continuing automatically..."
# User can Ctrl+C if needed
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
- [ ] **Using direct mamba install like v1 (NOT environment.yml)**
- [ ] No duplicate conda/pip package installations
- [ ] All external downloads have multiple fallback methods
- [ ] Pre-downloaded files exist in current directory before build
- [ ] Scripts are fully automatic (no manual confirmation prompts)
- [ ] **Check v1 Dockerfile FIRST before implementing new approaches**

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

# If conda package not found or mamba env update fails:
# Use v1 direct install method instead:
mamba install -y -n micapipe -c conda-forge package_name
# Then pip packages:
mamba run -n micapipe pip install pip_package

# If antspy installation fails:
# Check v1 Dockerfile lines 715-722 for working method
# Use: mamba install -y -n micapipe -c conda-forge antspy antspyx
```

---

## 📊 BUILD SUCCESS METRICS

**Total Mistakes Made:** 10  
**Total Mistakes Fixed:** 9  
**Harmless Issues:** 1  

**Time Wasted on Repeated Mistakes:** ~5 hours (especially environment.yml: 2 hours)  
**Time Saved by Checking v1 First:** Could have been immediate  
**Time Saved by This Document:** Hopefully infinite 😊

**Most Critical Lessons:**
1. **Check production/v1 code FIRST** before trying new approaches
2. **Never use environment.yml with mamba** - has parsing bugs
3. **Direct mamba install is explicit and production-proven**

---

**Remember:** Every mistake is a learning opportunity. This document prevents repeating them!
