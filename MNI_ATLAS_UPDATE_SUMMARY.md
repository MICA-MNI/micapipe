# MNI Atlas Update Summary - Issue #153

## Overview
Updated micapipe to use the MNI NLIN SYM 09a (non-linear symmetric 1x1x1mm) atlas instead of the previous MNI152 0.8mm linear asymmetric atlas. This change improves compatibility with LAMAReg and enhances transformation accuracy.

## Changes Made

### 1. Atlas Files
- **Created**: `MNI152Volumes/mni_icbm152_t1_tal_nlin_sym_09a_brain.nii.gz`
  - Brain-extracted version of the MNI NLIN SYM 09a atlas using the provided mask
  - Generated using Python nibabel: `atlas_data * mask_data`

### 2. functions/01_proc-structural.sh
- **Replaced**: Registration loop that used MNI152 0.8mm and 2mm templates
- **Updated**: Primary registration now uses MNI NLIN SYM 09a 1mm as main template
- **Added**: Legacy 2mm registration for backward compatibility
- **Changed**: Variable names from `MNI152_0.8mm` to `MNI152_1mm_nlinSym`
- **Updated**: Transformation variable references:
  - `mat_MNI152_SyN` → points to 1mm_nlinSym transformations
  - `MNI152_1mm_nlinSym` → new primary template variable
  - `MNI152_08mm` → removed, replaced with `MNI152_1mm_nlinSym`

### 3. functions/utilities.sh
- **Updated**: Export variables for MNI152 transformations
  - `mat_MNI152_SyN` → uses 1mm_nlinSym naming
  - `MNI152_mask` → points to new atlas mask
- **Updated**: JSON transformation descriptions
  - `from-nativepro_brain_to-MNI152_0.8mm` → `from-nativepro_brain_to-MNI152_1mm_nlinSym`
  - `from-MNI152_0p8mm_to-t1nativepro` → `from-MNI152_1mm_nlinSym_to-t1nativepro`
  - `BinaryMask_from` → updated to reflect new atlas

### 4. File Naming Convention
**New naming convention**: `MNI152_1mm_nlinSym` replaces `MNI152_0.8mm`

Transformation files now use:
- `*_from-nativepro_brain_to-MNI152_1mm_nlinSym_mode-image_desc-SyN_*`

## Atlas Information
- **Source**: MNI ICBM 152 Non-linear 6th Generation Symmetric Average Brain Stereotaxic Registration Model
- **Resolution**: 1x1x1mm isotropic
- **Type**: Non-linear symmetric (vs. previous linear asymmetric)
- **Files**:
  - `mni_icbm152_t1_tal_nlin_sym_09a.nii.gz` - Full atlas
  - `mni_icbm152_t1_tal_nlin_sym_09a_mask.nii.gz` - Brain mask
  - `mni_icbm152_t1_tal_nlin_sym_09a_brain.nii.gz` - Brain-extracted atlas (created)

## Compatibility Notes
1. **LAMAReg**: Now uses preferred MNI NLIN SYM 09a template for improved registration accuracy
2. **Backward Compatibility**: 2mm template registration maintained for legacy pipelines
3. **Cerebellum**: Still uses MNI152 0.8mm cerebellar atlas (no cerebellar atlas available for new template yet)

## Testing
- ✅ Syntax validation passed for both modified scripts
- ✅ Atlas files verified and properly formatted
- ✅ Brain-extracted atlas created successfully
- ✅ Transformation variable references updated consistently

## Reference
Fonov, V., Evans, A. C., Botteron, K., Almli, C. R., McKinstry, R. C., Collins, D. L., & Brain Development Cooperative Group. (2011). Unbiased average age-appropriate atlases for pediatric studies. Neuroimage, 54(1), 313-327.