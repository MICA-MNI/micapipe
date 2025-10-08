# FastSurfer PyTorch 2.6.0 Issue - Aggressive Fix + Debug Output

**Date**: October 8, 2025  
**Commit**: 67c5798  
**Issue**: `ERROR: No matching distribution found for torch==2.6.0+cpu`

## Problem Analysis

Despite patching the source YAML, `install_env.py` was still generating environment files with `torch==2.6.0+cpu` in the pip dependencies section. The error occurred during `conda env create` when it tried to install pip requirements.

## Root Cause

FastSurfer v2.4.2's `install_env.py` script may have HARDCODED PyTorch versions that override the source YAML, or the generated YAML has multiple sections where torch appears.

## Solution: Three-Layer Defense + Debug Output

### Layer 1: Source YAML Patching (Before install_env.py)
```dockerfile
sed -i 's/torch>=2\.0/torch==2.5.1/g' /opt/FastSurfer/env/fastsurfer.yml
sed -i 's/torchvision>=0\.15/torchvision==0.20.1/g' /opt/FastSurfer/env/fastsurfer.yml
```

### Layer 2: Debug Output (See what was generated)
```dockerfile
echo "✅ Generated CPU YAML - checking PyTorch versions:"
grep -E "torch|torchvision" /opt/FastSurfer/fastsurfer_cpu.yml
echo "--- Pip dependencies ---"
grep -A 50 "pip:" /opt/FastSurfer/fastsurfer_cpu.yml | head -60
```

**Purpose**: Shows us EXACTLY what `install_env.py` generated, so we can see if torch 2.6.0 appears and where.

### Layer 3: Aggressive Patching (Catch ALL variations)
```dockerfile
# Patch BOTH with and without +cpu suffix
sed -i 's/torch==2\.6\.0+cpu/torch==2.5.1+cpu/g' /opt/FastSurfer/fastsurfer_cpu.yml
sed -i 's/torch==2\.6\.0/torch==2.5.1/g' /opt/FastSurfer/fastsurfer_cpu.yml  # ← NEW

sed -i 's/torchvision==0\.21\.0+cpu/torchvision==0.20.1+cpu/g' /opt/FastSurfer/fastsurfer_cpu.yml
sed -i 's/torchvision==0\.21\.0/torchvision==0.20.1/g' /opt/FastSurfer/fastsurfer_cpu.yml  # ← NEW
```

**Key Improvement**: Previous fix only patched `torch==2.6.0\+cpu` (with escaped +). If the YAML contains `torch==2.6.0` WITHOUT the +cpu suffix in some sections, it would be missed.

### Layer 4: Verification (Confirm patches worked)
```dockerfile
echo "✅ Patches applied - Final versions:"
grep -E "torch|torchvision" /opt/FastSurfer/fastsurfer_cpu.yml | head -20
```

## Positioning Change

**Moved from line 680 → line 251** (right after micapipe environment creation)

**Why**:
- Fail fast strategy - catch PyTorch issues early
- Avoid wasting 60+ minutes on FSL/FreeSurfer if FastSurfer will fail
- Debug output appears early in build logs

## What to Look For in Build Logs

When the build runs, you'll see:

```
📦 Patching FastSurfer source YAML to use PyTorch 2.5.1...
✅ Source YAML patched:
  - pip:
    ...
    - torch==2.5.1
    - torchvision==0.20.1

📦 Generating CPU-only FastSurfer environment...
✅ Generated CPU YAML - checking PyTorch versions:
[This will show if install_env.py ignored our patches]

--- Pip dependencies ---
[This will show the full pip section]

🔧 Applying patches to CPU YAML (with and without +cpu suffix)...
✅ Patches applied - Final versions:
[This MUST show torch==2.5.1+cpu, not 2.6.0]

🚀 Creating CPU FastSurfer environment...
```

## Expected Outcomes

### If Source Patches Work
- Generated YAML will already have torch 2.5.1
- Layer 3 patches will have nothing to change
- Build proceeds successfully

### If install_env.py Overrides Source
- Debug output will show torch 2.6.0 in generated YAML
- Layer 3 patches will fix it (both +cpu and non-+cpu variants)
- Verification will confirm torch 2.5.1 is present
- Build proceeds successfully

### If Build Still Fails
Debug output will reveal:
1. Where torch 2.6.0 appears in the YAML
2. What format it uses (2.6.0+cpu vs 2.6.0 vs other)
3. If there are multiple torch entries we're missing

## Manual Verification on Server

After pulling the changes:

```bash
cd ~/micapipe
git pull origin comprehensive-base-image
git log -1 --oneline  # Should show: 67c5798 fix: Move FastSurfer VERY early...

# Verify FastSurfer is now early in Dockerfile
grep -n "FastSurfer - MOVED VERY EARLY" Dockerfile.base
# Should show: line ~252

# Verify aggressive patching
grep "torch==2\.6\.0+cpu\|torch==2\.6\.0" Dockerfile.base
# Should show TWO sed commands (one with +cpu, one without)
```

## Rebuild Command

```bash
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh 2>&1 | tee rebuild_$(date +%Y%m%d_%H%M%S).log
```

Watch for the debug output in the early stages (around 5-10 minutes into the build).

## If This STILL Doesn't Work

The debug output will tell us:
1. Exact format of torch version in generated YAML
2. If it appears in multiple places (conda deps + pip deps)
3. If install_env.py is using a separate requirements.txt

Then we can:
- Patch install_env.py itself to change hardcoded versions
- Or skip install_env.py and manually create the environment YAML
- Or use FastSurfer's official Docker approach with our own patches

## Success Criteria

✅ Build progresses past FastSurfer environment creation  
✅ Debug output shows torch 2.5.1+cpu after patches  
✅ No error about torch==2.6.0+cpu not found  
✅ FastSurfer environment has PyTorch 2.5.1 installed  

---

**Status**: Ready for rebuild with enhanced debugging and aggressive patching
