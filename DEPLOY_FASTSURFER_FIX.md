# Quick Deployment Guide - FastSurfer Aggressive Fix

**Branch**: comprehensive-base-image  
**Latest Commit**: 2c78f3a  
**Date**: October 8, 2025

## What Changed

✅ **FastSurfer moved VERY early** (line 251, right after micapipe environment)  
✅ **Debug output added** to see what `install_env.py` generates  
✅ **Aggressive patching** - catches `torch==2.6.0` with AND without `+cpu` suffix  
✅ **Verification output** shows final versions before environment creation  

## Deploy to Server

```bash
# 1. Pull latest changes
cd ~/micapipe
git pull origin comprehensive-base-image

# 2. Verify you have the fix
git log -1 --oneline
# Should show: 2c78f3a docs: Add detailed explanation of aggressive FastSurfer fix

# 3. Sync to server build location
./migrate_comprehensive_base_to_server.sh

# 4. Go to server build directory
cd /host/cassio/export03/data/enning/downloads

# 5. Start the build (will take 45-90 minutes)
./build_base_image_server.sh 2>&1 | tee rebuild_$(date +%Y%m%d_%H%M%S).log
```

## What to Watch For

The build will now show detailed output around lines 250-320:

### 1. Source YAML Patching (Should see)
```
📦 Patching FastSurfer source YAML to use PyTorch 2.5.1...
✅ Source YAML patched:
  - torch==2.5.1
  - torchvision==0.20.1
```

### 2. Generated YAML (Shows what install_env.py created)
```
📦 Generating CPU-only FastSurfer environment...
✅ Generated CPU YAML - checking PyTorch versions:
--- Pip dependencies ---
```

**LOOK HERE**: If you see `torch==2.6.0` in this output, that's the problem we're fixing!

### 3. Aggressive Patching (Fixes both variations)
```
🔧 Applying patches to CPU YAML (with and without +cpu suffix)...
✅ Patches applied - Final versions:
  - torch==2.5.1+cpu      ← Must show 2.5.1, NOT 2.6.0
  - torchvision==0.20.1+cpu
```

### 4. Environment Creation (The moment of truth)
```
🚀 Creating CPU FastSurfer environment...
Collecting packages...
Installing pip dependencies...
```

**SUCCESS**: If you see "Installing pip dependencies..." complete without error, the fix worked!  
**FAILURE**: If you see "ERROR: No matching distribution found for torch==2.6.0+cpu", check the debug output above.

## If Build Succeeds

🎉 You're done! The base image is built with:
- ✅ FastSurfer with PyTorch 2.5.1
- ✅ DESIGNER with proper dependencies
- ✅ All other tools (FSL, FreeSurfer, ANTs, etc.)

## If Build Still Fails

The debug output will reveal:
1. What version `install_env.py` generated
2. Where in the YAML file torch appears
3. What format it uses (with/without +cpu suffix)

Check the build log and look for the section between these markers:
- `📦 Generating CPU-only FastSurfer environment...`
- `🚀 Creating CPU FastSurfer environment...`

Copy that section and we can see exactly what's happening.

## Quick Status Check

After pulling changes, verify:

```bash
# FastSurfer is early in Dockerfile
grep -n "FastSurfer - MOVED VERY EARLY" Dockerfile.base
# Output: 252:# FastSurfer - MOVED VERY EARLY TO FAIL FAST

# Aggressive patching is present
grep -c "torch==2\.6\.0" Dockerfile.base
# Output: 2 (one with +cpu, one without)

# Debug output is present
grep "Generated CPU YAML - checking" Dockerfile.base
# Should find the echo statement
```

## Expected Timeline

- **Minutes 0-5**: Base system setup, conda installation
- **Minutes 5-10**: micapipe environment creation
- **Minutes 10-15**: FastSurfer (EARLY - this is where we fail fast if there's an issue)
- **Minutes 15-30**: DESIGNER installation
- **Minutes 30-60**: FSL installation (large package)
- **Minutes 60-80**: FreeSurfer installation (large package)
- **Minutes 80-90**: Final tools, cleanup

FastSurfer now fails in minutes 10-15 instead of minutes 60-70, saving 50+ minutes per failed build!

---

**Ready to deploy**: All changes committed and pushed to comprehensive-base-image branch
