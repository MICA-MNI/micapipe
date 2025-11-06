# LAMAReg Two-Stage Registration Update

## Summary
Updated all LAMAReg registration commands to properly handle the two-stage robust registration mode by adding `--secondary-warpfield` and `--inverse-secondary-warpfield` arguments. This prevents automatic warpfield composition, which can cause accuracy loss.

## Why This Change Is Necessary

According to LAMAReg documentation:
- **Robust mode (default)**: LAMAReg performs two-stage registration for improved accuracy
  1. **Stage 1**: Register parcellations (contrast-agnostic, coarse alignment) → produces primary warpfield + affine
  2. **Stage 2**: Fine-tune with direct image registration using Stage 1 as initialization → produces secondary warpfield (refinement)
  3. **Final transform**: Composition of both warpfields (primary_warp ∘ secondary_warp)

- **Warpfield Composition**: 
  - Default: keeps primary and secondary warpfields separate (better accuracy, no composition loss)
  - `--compose` flag: writes a single composed warp (faster application but loses accuracy)
  - **Critical**: Providing exactly one warpfield path without `--compose` is invalid

By specifying both `--secondary-warpfield` and `--inverse-secondary-warpfield`, we ensure:
- No automatic composition (preserves accuracy)
- Explicit control over transformation application order
- Compatibility with ANTs transform chains

## Files Modified

### 1. functions/02_proc-dwi.sh
**Purpose**: DWI preprocessing and registration to T1w

**Changes**:
- Added secondary warpfield paths: `dwi_SyN_warp2`, `dwi_SyN_Invwarp2`
- Added `--secondary-warpfield` and `--inverse-secondary-warpfield` to LAMAReg command
- Updated transformation chains:
  - `trans_T12dwi`: Added `-t ${dwi_SyN_warp2}` at the beginning
  - `trans_dwi2T1`: Added `-t ${dwi_SyN_Invwarp2}` at the end

### 2. functions/02_proc-func.sh
**Purpose**: Functional MRI preprocessing and registration to T1nativepro

**Changes**:
- Added secondary warpfield paths: `SyN_func_warp2`, `SyN_func_Invwarp2`
- Added `--secondary-warpfield` and `--inverse-secondary-warpfield` to LAMAReg command
- Updated transformation chains:
  - `transformsInv`: Added `-t ${SyN_func_warp2}` at the beginning (T1nativepro → func)
  - `transform`: Added `-t ${SyN_func_Invwarp2}` at the end (func → T1nativepro)

### 3. functions/02_proc-flair.sh
**Purpose**: FLAIR image processing and registration to T1nativepro

**Changes**:
- Added secondary warpfield paths: `flair_warpfield2`, `flair_inv_warpfield2`
- Added `--secondary-warpfield` and `--inverse-secondary-warpfield` to LAMAReg command
- Updated `antsApplyTransforms` to include both warpfields: `-t "$flair_warpfield2" -t "$flair_warpfield"`
- Updated JSON metadata to document both warpfields in the transformation command

### 4. functions/03_MPC.sh
**Purpose**: Microstructural profile covariance (MPC) processing - registration to FreeSurfer space

**Changes**:
- Added secondary warpfield paths: `SyN_qMRI2fs_warp2`, `SyN_qMRI2fs_Invwarp2`
- Added `--secondary-warpfield` and `--inverse-secondary-warpfield` to LAMAReg command
- Updated transformation chains:
  - `transformsInv`: Added `-t ${SyN_qMRI2fs_Invwarp2}` at the end (T1_fsnative → qMRI)
  - `transforms`: Added `-t ${SyN_qMRI2fs_warp2}` at the beginning (qMRI → T1_fsnative)

### 5. functions/03_MPC-SWM.sh
**Purpose**: MPC processing for superficial white matter (SWM) - registration to T1nativepro

**Changes**:
- Added secondary warpfield paths: `SyN_qMRI2np_warp2`, `SyN_qMRI2np_Invwarp2`
- Added `--secondary-warpfield` and `--inverse-secondary-warpfield` to LAMAReg command
- Updated transformation chains:
  - `transformsInv`: Added `-t ${SyN_qMRI2np_Invwarp2}` at the end (T1nativepro → qMRI)
  - `transforms`: Added `-t ${SyN_qMRI2np_warp2}` at the beginning (qMRI → T1nativepro)

## Output File Naming Convention

LAMAReg now generates the following files for each registration:

### Primary Stage (Parcellation-based):
- `{prefix}0GenericAffine.mat` - Affine transformation matrix
- `{prefix}1Warp.nii.gz` - Primary warpfield (forward)
- `{prefix}1InverseWarp.nii.gz` - Primary warpfield (inverse)

### Secondary Stage (Direct image registration):
- `{prefix}2Warp.nii.gz` - Secondary warpfield (forward, refinement)
- `{prefix}2InverseWarp.nii.gz` - Secondary warpfield (inverse, refinement)

### Parcellation Outputs:
- `{prefix}_fixed_parc.nii.gz` - Fixed image parcellation (SynthSeg)
- `{prefix}_moving_parc.nii.gz` - Moving image parcellation (SynthSeg)
- `{prefix}_registered_parc.nii.gz` - Registered parcellation
- `{prefix}_dice_scores.csv` - QC metrics (DICE scores)

### Final Warped Image:
- `{prefix}Warped.nii.gz` - Registered moving image

## Transformation Application Order

**Important**: ANTs applies transformations in **reverse order** (right to left)

### Forward Transforms (moving → fixed):
```bash
antsApplyTransforms \
  -t secondary_warp.nii.gz \  # Applied SECOND (refinement)
  -t primary_warp.nii.gz \     # Applied FIRST (coarse alignment)
  -t affine.mat                # Applied BEFORE warps (linear alignment)
```

### Inverse Transforms (fixed → moving):
```bash
antsApplyTransforms \
  -t affine.mat \                        # Applied THIRD (linear, inverted)
  -t primary_inverse_warp.nii.gz \       # Applied SECOND (coarse, inverted)
  -t secondary_inverse_warp.nii.gz       # Applied FIRST (refinement, inverted)
```

## Benefits of This Approach

1. **Preserves Registration Accuracy**: No composition-related accuracy loss
2. **Explicit Control**: Clear visibility of transformation stages
3. **Debugging**: Each stage can be inspected independently
4. **Flexibility**: Can apply transformations to multiple images efficiently
5. **QC Integration**: DICE scores available for each registration stage
6. **Reproducibility**: Deterministic behavior with explicit file paths

## Testing Checklist

Before merging, verify:
- ✅ All 5 modules compile without syntax errors
- ✅ LAMAReg generates both primary and secondary warpfields
- ✅ QC CSV files are created with DICE scores
- ✅ Registered images show improved alignment vs single-stage
- ✅ Transformation application order is correct (forward/inverse)
- ✅ File naming follows LAMAReg convention
- ✅ JSON metadata documents all transformations correctly

## Related Documentation

- **LAMAReg GitHub**: https://github.com/MICA-MNI/LAMAReg
- **Issue**: #156 - Update regSynth to LAMAReg
- **Previous Commit**: 3d7fa0e - "Replace regSynth with LAMAReg as default cross-modality registration"

## Command Example

```bash
# Full LAMAReg two-stage registration
lamareg register \
  --moving moving_image.nii.gz \
  --fixed fixed_image.nii.gz \
  --output registered_output.nii.gz \
  --moving-parc moving_parc.nii.gz \
  --fixed-parc fixed_parc.nii.gz \
  --registered-parc registered_parc.nii.gz \
  --affine transform_0GenericAffine.mat \
  --warpfield transform_1Warp.nii.gz \
  --inverse-warpfield transform_1InverseWarp.nii.gz \
  --secondary-warpfield transform_2Warp.nii.gz \
  --inverse-secondary-warpfield transform_2InverseWarp.nii.gz \
  --qc-csv dice_scores.csv \
  --synthseg-threads 4 \
  --ants-threads 8
```

## Statistics

- **Files Changed**: 5
- **Insertions**: +31
- **Deletions**: -11
- **Net Change**: +20 lines
