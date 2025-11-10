# Instructions: Fix DWI Test Data on Server

**Date:** November 10, 2025  
**Issue:** DWI test fails because FOD is 4D, but SynthSeg needs 3D images

## Quick Fix on Server

### 1. Pull Latest Changes
```bash
cd ~/micapipe
git pull origin 156-update-regsynth-to-lamareg
```

### 2. Extract 3D b0 Image from FOD (or find existing b0)

**Option A: If you have raw DWI data with b0 volumes**
```bash
cd /data_/mica3/BIDS_CI/lamareg_test_data/dwi

# Extract b0 from DWI
fslroi <path_to_dwi>/dwi.nii.gz dwi_b0.nii.gz 0 1

# Or if b0 already exists somewhere:
cp <path_to_existing_b0>/b0.nii.gz dwi_b0.nii.gz
```

**Option B: If you have FA map**
```bash
cd /data_/mica3/BIDS_CI/lamareg_test_data/dwi

# Copy FA map as DWI reference
cp <path_to_FA>/FA.nii.gz dwi_FA.nii.gz
```

**Option C: Extract first volume from 4D FOD (quick workaround)**
```bash
cd /data_/mica3/BIDS_CI/lamareg_test_data/dwi

# Extract first SH coefficient (isotropic component)
fslroi dwi_fod.nii.gz dwi_fod_vol0.nii.gz 0 1
```

### 3. Re-run Test
```bash
cd ~/micapipe/tests/lamareg_tests
./test_dwi_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/dwi
```

### 4. What the Test Will Do Now

The updated test automatically checks for 3D images in this order:
1. `dwi_b0.nii.gz` (preferred)
2. `dwi_FA.nii.gz` (good alternative)
3. `dwi_fod_vol0.nii.gz` (extracted from 4D FOD)
4. `dwi_fod.nii.gz` (will warn and fail)

If it finds a 3D image, it will use that instead of the 4D FOD.

## Expected Output (Success)

```bash
==========================================
Checking input files...
==========================================
Moving image: /data_/mica3/BIDS_CI/lamareg_test_data/dwi/T1w_in_dwi_brain.nii.gz
Fixed image: /data_/mica3/BIDS_CI/lamareg_test_data/dwi/dwi_b0.nii.gz

✓ PASS: Input: Moving image
✓ PASS: Input: Fixed image
✓ All input files present - proceeding with LAMAReg registration
Fixed image is 3D - OK for SynthSeg
✓ PASS: Fixed image dimensionality

==========================================
EXECUTING LAMAReg REGISTRATION
==========================================
[timestamp] Starting LAMAReg registration (this may take 10-15 minutes on CPU)...
[timestamp] Moving: T1w_in_dwi_brain.nii.gz
[timestamp] Fixed: dwi_b0.nii.gz
[timestamp] Output prefix: ./test_output/dwi_to_T1w_

<... LAMAReg runs for 10-15 minutes ...>

==========================================
[timestamp] LAMAReg registration completed successfully in 950 seconds
✓ PASS: LAMAReg execution
[timestamp] Output files: 10 created, 0 missing
==========================================
```

## Re-run Setup Script (Alternative)

If you want to regenerate all test data with proper 3D extraction:

```bash
cd ~/micapipe/tests/lamareg_tests
./setup_hcp_test_data.sh
```

The updated setup script will:
- Look for existing b0 images
- Use FA maps if b0 not found
- Extract first volume from FOD if neither available
- Warn clearly about what it's doing

## What Still Needs Fixing

### Main Pipeline Bug
The main pipeline (`functions/02_proc-dwi.sh`) has the same bug:

```bash
# Line 511-513 (CURRENT - WRONG):
Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${fod}" \                      # ← 4D FOD will fail!

# SHOULD BE:
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0_dwi.nii.gz"

Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${dwi_b0}" \                   # ← 3D b0 image
```

This needs to be fixed in the main pipeline code, not just the test!

## Summary

**Problem:** FOD is 4D (SH coefficients), SynthSeg needs 3D anatomical images

**Solution:** Use b0, FA, or extract first FOD volume

**Status:** 
- ✅ Test scripts fixed and pushed
- ✅ Setup script updated
- ✅ Documentation created
- ❌ Main pipeline (02_proc-dwi.sh) still needs fixing

**Next Steps:**
1. Extract/find 3D DWI reference on server
2. Re-run test to verify it works
3. Fix main pipeline to use b0 instead of FOD
4. Test full pipeline on real data
