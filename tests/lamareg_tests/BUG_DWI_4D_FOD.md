# Critical Bug: DWI Registration Using 4D FOD with SynthSeg

**Date:** November 10, 2025  
**Severity:** HIGH  
**Status:** BLOCKING DWI REGISTRATION TESTS

## Problem

The DWI registration in `02_proc-dwi.sh` uses a **4D FOD image** as the fixed reference for LAMAReg registration:

```bash
# Line 492
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"

# Line 511-513
Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${fod}" \           # ← 4D FOD image!
```

**This fails because SynthSeg (used internally by LAMAReg) expects 3D anatomical images.**

## Error Message

```
ValueError: operands could not be broadcast together with shapes (4,) (3,)
```

This occurs in SynthSeg when it tries to process a 4D image (FOD with multiple SH coefficients) as a 3D anatomical scan.

## Root Cause

**FOD (Fiber Orientation Distribution) images are 4D:**
- dim1, dim2, dim3: spatial dimensions
- dim4: spherical harmonic (SH) coefficients (typically 45-120 volumes)

**SynthSeg parcellation requires 3D anatomical images:**
- Works with: T1w, T2w, FLAIR, b0, FA maps
- Does NOT work with: 4D images, time series, multi-shell DWI, FOD

## Solution Options

### Option 1: Use b0 Image (Recommended)
The b0 image is a 3D anatomical reference from DWI:

```bash
# Use b0 instead of FOD
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0_dwi.nii.gz"

Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${dwi_b0}" \    # ← 3D b0 image
  --output "$dwi_SyN_output" \
  ...
```

### Option 2: Use FA Map
FA (Fractional Anisotropy) is a 3D scalar map derived from DWI:

```bash
# Use FA map
dwi_FA="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.nii.gz"

Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${dwi_FA}" \    # ← 3D FA map
  --output "$dwi_SyN_output" \
  ...
```

### Option 3: Extract First Volume of FOD
As a workaround, extract the first SH coefficient (which represents the isotropic component):

```bash
# Extract first volume from 4D FOD
dwi_fod_vol0="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_vol0.nii.gz"
Do_cmd fslroi "${fod}" "${dwi_fod_vol0}" 0 1

Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${dwi_fod_vol0}" \    # ← First volume only
  --output "$dwi_SyN_output" \
  ...
```

## Recommended Fix

**Use the b0 image** (Option 1) because:
1. ✅ Already computed in the DWI pipeline
2. ✅ 3D anatomical reference in DWI space
3. ✅ High contrast between brain and background
4. ✅ Standard practice for DWI-to-structural registration
5. ✅ Compatible with SynthSeg parcellation

## Code Changes Needed

### In `02_proc-dwi.sh`

```bash
# BEFORE (line 511-513):
Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${fod}" \

# AFTER:
# Use b0 image for registration (3D reference in DWI space)
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0_dwi.nii.gz"

Do_cmd lamareg register \
  --moving "${T1nativepro_in_dwi_brain}" \
  --fixed "${dwi_b0}" \
```

## Test Data Preparation

For testing, create a 3D DWI reference image:

```bash
# Option 1: Extract b0 from DWI data
fslroi dwi.nii.gz dwi_b0.nii.gz 0 1

# Option 2: Use existing FA map
cp DTI_FA.nii.gz dwi_FA.nii.gz

# Option 3: Extract first FOD volume
fslroi dwi_fod.nii.gz dwi_fod_vol0.nii.gz 0 1
```

Then provide in test data directory:
```
lamareg_test_data/dwi/
  ├── T1w_in_dwi_brain.nii.gz   # Moving image (3D)
  └── dwi_b0.nii.gz              # Fixed image (3D) ← NEW!
```

## Impact Assessment

### Files Affected:
- ✅ `functions/02_proc-dwi.sh` - **NEEDS FIX**
- ✅ `tests/lamareg_tests/test_dwi_registration.sh` - **UPDATED** with workaround
- ✅ Test data setup scripts - **NEEDS UPDATE**

### Registration Quality:
- Using b0 or FA may actually **improve registration quality** compared to FOD
- b0 has better anatomical contrast for structural alignment
- FA highlights white matter structure (good for DWI-specific registration)

### Backward Compatibility:
- Output file names remain the same
- Transformation chains remain the same
- Only the fixed image reference changes

## Testing Status

- ❌ **DWI test currently FAILS** due to 4D FOD input
- ✅ **Test script updated** to detect and warn about 4D images
- ✅ **Test script** will auto-select 3D image if available (b0, FA, mean_b0)
- ⏳ **Main pipeline needs fixing** in `02_proc-dwi.sh`

## Next Steps

1. **Immediate:** Update `02_proc-dwi.sh` to use `dwi_b0` instead of `fod`
2. **Testing:** Prepare test data with 3D b0 image
3. **Validation:** Re-run DWI registration tests
4. **Documentation:** Update pipeline docs to reflect b0 usage

## References

- LAMAReg uses SynthSeg for parcellation-based QC
- SynthSeg paper: "Robust machine learning segmentation" (Billot et al.)
- SynthSeg works with 3D anatomical images only
- Standard DWI registration practice: use b0 or FA as reference

---

**Assignee:** @enningyang  
**Priority:** HIGH (blocks testing and production use)  
**Branch:** 156-update-regsynth-to-lamareg
