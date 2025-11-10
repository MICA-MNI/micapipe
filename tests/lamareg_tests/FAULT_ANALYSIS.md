# Analysis: Who's At Fault for the 4D FOD Bug?

**Date:** November 10, 2025  
**Verdict:** **PARTIALLY MY FAULT** (test setup issue, not pipeline issue)

## What Actually Happened

### The OLD Code (RegSynth) ✅ WAS CORRECT
```bash
# Line 492-493 in OLD code
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
Do_cmd mrconvert -coord 3 0 "$fod_wmN" "$fod"  # Extract first volume → 3D!

# Line ~506 - Use the extracted 3D volume
Do_cmd antsRegistrationSyN.sh ... -f "$fod" ...  # $fod is 3D here
```

**The old code extracted volume 0** from the 4D FOD using `mrconvert -coord 3 0`, making it 3D!

### My LAMAReg Code ✅ IS ALSO CORRECT
```bash
# Line 492-493 in CURRENT code (UNCHANGED!)
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
Do_cmd mrconvert -coord 3 0 "$fod_wmN" "$fod"  # Still extracting first volume!

# Line 511-513 - Use the extracted 3D volume  
Do_cmd lamareg register ... --fixed "${fod}" ...  # $fod is 3D here too!
```

**I KEPT the mrconvert line**, so the pipeline code is correct!

## Where I Made the Mistake ❌

### The Problem: Test Data Setup

**In the test, I did this:**
```bash
# setup_hcp_test_data.sh (OLD version)
cp "${SOURCE}/dwi/sub-211215_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" \
   "${TEST_DIR}/dwi/dwi_fod.nii.gz"  # Copied 4D FOD directly!
```

**But the pipeline expects this:**
```bash
# The pipeline creates this 3D extracted FOD
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"  # 3D after mrconvert
```

### The Real Issue

The **test data** used the **4D source FOD** directly, but the **pipeline** uses a **3D extracted FOD** (first volume only).

**Test had:** `/test_data/dwi_fod.nii.gz` (4D, 120 volumes)  
**Pipeline expects:** `$fod` (3D, after `mrconvert -coord 3 0`)

## Fault Assignment

### ✅ Pipeline Code: NOT AT FAULT
- Old code correctly extracted 3D volume
- New LAMAReg code kept the same extraction
- Both would work fine with proper data

### ❌ Test Setup: MY FAULT
- I copied the 4D source FOD to test data
- Should have copied the EXTRACTED 3D FOD (after mrconvert)
- Or should have extracted it in the test setup script

### ❌ Test Script: PARTIALLY MY FAULT
- Test used `dwi_fod.nii.gz` directly without checking dimensions
- Should have mimicked the pipeline's extraction step
- Should have validated input is 3D before passing to LAMAReg

## Why It Failed

**Pipeline flow (CORRECT):**
```
4D FOD (fod_wmN)
  → mrconvert -coord 3 0  
  → 3D FOD (fod)
  → lamareg register --fixed $fod (3D) ✅
```

**Test flow (WRONG):**
```
4D FOD (dwi_fod.nii.gz)
  → test_dwi_registration.sh
  → lamareg register --fixed dwi_fod.nii.gz (4D) ❌
```

## The Fix Applied

### 1. Fixed Test Script ✅
```bash
# Now checks for 3D alternatives first:
- dwi_b0.nii.gz (preferred)
- dwi_FA.nii.gz  
- dwi_fod_vol0.nii.gz (extracted)
- dwi_fod.nii.gz (4D, will warn and fail)
```

### 2. Fixed Setup Script ✅
```bash
# Now extracts 3D volume:
if [ -f "${SOURCE}/dwi/FOD_4D.nii.gz" ]; then
    fslroi "${SOURCE}/dwi/FOD_4D.nii.gz" \
           "${TEST_DIR}/dwi/dwi_fod_vol0.nii.gz" 0 1  # Extract vol 0!
fi
```

## Lessons Learned

### What I Should Have Done:
1. ✅ Check what the pipeline ACTUALLY does before creating tests
2. ✅ Mimic the pipeline's data processing (mrconvert extraction)
3. ✅ Validate test data matches pipeline expectations
4. ✅ Test with realistic intermediate files, not raw source data

### What Went Right:
1. ✅ I didn't break the pipeline code (mrconvert line still there)
2. ✅ LAMAReg integration is correct
3. ✅ The bug was caught immediately by testing
4. ✅ Fix is straightforward (use 3D data)

## Conclusion

**Verdict:** This is **MY FAULT**, specifically in test data preparation, not in the LAMAReg integration itself.

- ✅ **LAMAReg integration code:** CORRECT (kept mrconvert extraction)
- ❌ **Test data setup:** WRONG (used 4D source instead of 3D extracted)
- ❌ **Test script assumptions:** WRONG (didn't validate dimensionality)

**The good news:** The pipeline code itself is fine. Only the test infrastructure needs fixing, which is now done.

**Apology:** I should have been more careful to understand the full data flow before setting up tests. Thank you for catching this! 🙏
