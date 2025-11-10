# LAMAReg Testing Guide for Server

**Last Updated:** November 10, 2025  
**Branch:** 156-update-regsynth-to-lamareg

## Quick Start

### 1. Pull Latest Code
```bash
cd ~/micapipe
git pull origin 156-update-regsynth-to-lamareg
```

### 2. Activate LAMAReg Environment
```bash
conda activate lamareg
# Verify LAMAReg is available
lamareg --version
```

### 3. Prepare Test Data (CRITICAL!)

The DWI test needs a **3D reference image**, not the 4D FOD. Extract it:

```bash
cd /data_/mica3/BIDS_CI/lamareg_test_data/dwi

# Option A: Use mrconvert (MRtrix3 - already installed in LAMAReg env)
mrconvert -coord 3 0 dwi_fod.nii.gz dwi_fod_vol0.nii.gz

# Option B: Use fslroi (FSL 6.0.3 on server)
source /data_/mica1/01_programs/fsl-6-0-3/etc/fslconf/fsl.sh
fslroi dwi_fod.nii.gz dwi_fod_vol0.nii.gz 0 1

# Option C: Use existing b0 or FA if available (best)
# Find and copy: dwi_b0.nii.gz or dwi_FA.nii.gz
```

### 4. Run Tests

```bash
cd ~/micapipe/tests/lamareg_tests

# Option A: Run all tests at once (RECOMMENDED)
./run_all_tests.sh /data_/mica3/BIDS_CI/lamareg_test_data ./test_results 16

# Option B: Run individual tests (if you prefer)
# Test 1: DWI registration (10-15 minutes)
./test_dwi_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/dwi

# Test 2: Functional registration (10-15 minutes)
./test_func_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/func

# Test 3: MPC registration (10-15 minutes)
./test_mpc_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/mpc

# Test 4: MPC-SWM registration (10-15 minutes)
./test_mpc_swm_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/mpc
```

**Note:** The default thread count is 16 (4 for SynthSeg, 12 for ANTs). You can adjust by setting `LAMAREG_THREADS`:
```bash
# Use 8 threads instead
./run_all_tests.sh /data_/mica3/BIDS_CI/lamareg_test_data ./test_results 8

# Or set environment variable for individual tests
export LAMAREG_THREADS=8
./test_dwi_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/dwi
```

---

## Understanding Test Output

### ✅ SUCCESS - LAMAReg Actually Running:
```
==========================================
Checking input files...
==========================================
Moving image: /path/to/moving.nii.gz
Fixed image: /path/to/fixed.nii.gz

✓ PASS: Input: Moving image
✓ PASS: Input: Fixed image
✓ All input files present - proceeding with LAMAReg registration

==========================================
Checking required tools...
==========================================
✓ PASS: LAMAReg installation
LAMAReg version: 1.4.1
✓ PASS: ANTs installation
✓ All required tools are available

==========================================
EXECUTING LAMAReg REGISTRATION
==========================================
[2025-11-10 14:30:15] Starting LAMAReg registration (this may take 10-15 minutes on CPU)...
[2025-11-10 14:30:15] Moving: T1w_in_dwi_brain.nii.gz
[2025-11-10 14:30:15] Fixed: dwi_b0.nii.gz

<... LAMAReg runs for ~10-15 minutes with lots of output ...>

==========================================
[2025-11-10 14:45:30] LAMAReg registration completed successfully in 915 seconds
✓ PASS: LAMAReg execution
[2025-11-10 14:45:30] Output files: 10 created, 0 missing
==========================================

========================================
Test Summary
========================================
Total tests: 5
Passed: 5
Failed: 0

All tests passed!
```

**Key indicators it's working:**
- ✓ Shows "EXECUTING LAMAReg REGISTRATION" banner
- ✓ Takes 10-15 minutes (900+ seconds)
- ✓ Shows "completed successfully in XXX seconds"
- ✓ Shows "10 created, 0 missing" files

### ❌ FAILURE - Missing Input Files:
```
==========================================
Checking input files...
==========================================
Moving image: /path/to/missing.nii.gz
Fixed image: /path/to/missing.nii.gz

✗ FAIL: Input: Moving image - File not found: /path/to/missing.nii.gz
ERROR: Moving image file is missing!
SKIPPING LAMAReg execution - cannot proceed without input files

========================================
Test Summary
========================================
Total tests: 2
Passed: 1
Failed: 1
```

**Key indicators of skipping:**
- ✗ Shows "ERROR: ... file is missing!"
- ✗ Shows "SKIPPING LAMAReg execution"
- ✗ Completes in < 5 seconds

### ❌ FAILURE - 4D Image Used (SynthSeg Error):
```
ValueError: operands could not be broadcast together with shapes (4,) (3,)
Error during processing: Command '['lamareg', 'synthseg', ...]' returned non-zero exit status 1.
✗ FAIL: LAMAReg execution - Registration failed (exit code: 1)
```

**Solution:** Extract 3D volume from 4D FOD (see step 3 above)

---

## Required Test Data Files

### DWI Test (`/data_/mica3/BIDS_CI/lamareg_test_data/dwi/`)
```
T1w_in_dwi_brain.nii.gz    # Moving image (REQUIRED)
dwi_b0.nii.gz              # Fixed - 3D b0 (BEST)
  OR
dwi_FA.nii.gz              # Fixed - 3D FA map (GOOD)
  OR
dwi_fod_vol0.nii.gz        # Fixed - 3D extracted from FOD (OK)
```

### Functional Test (`/data_/mica3/BIDS_CI/lamareg_test_data/func/`)
```
func_brain.nii.gz          # Functional MRI (brain extracted)
T1_nativepro.nii.gz        # T1 native processed
```

### MPC Test (`/data_/mica3/BIDS_CI/lamareg_test_data/mpc/`)
```
microstructural_map.nii.gz # qMRI map (e.g., FA)
T1_fsnative.nii.gz         # FreeSurfer native T1
```

### MPC-SWM Test (`/data_/mica3/BIDS_CI/lamareg_test_data/mpc_swm/`)
```
microstructural_map.nii.gz # qMRI map
T1_nativepro.nii.gz        # T1 native processed
```

---

## Checking Test Results

### Log Files
Each test creates output in `test_output_*/`:
```bash
# View full LAMAReg output
cat test_output_*/test_log.txt

# View test pass/fail summary
cat test_output_*/test_results.txt

# Check output files created
ls -lh test_output_*/dwi_to_T1w_*
```

### Expected Outputs (per test)
```
dwi_to_T1w_0GenericAffine.mat         # Affine transformation
dwi_to_T1w_1Warp.nii.gz                # Primary warpfield
dwi_to_T1w_1InverseWarp.nii.gz        # Inverse primary warp
dwi_to_T1w_2Warp.nii.gz                # Secondary warpfield
dwi_to_T1w_2InverseWarp.nii.gz        # Inverse secondary warp
dwi_to_T1w__moving_parc.nii.gz        # Moving parcellation
dwi_to_T1w__fixed_parc.nii.gz         # Fixed parcellation
dwi_to_T1w__registered_parc.nii.gz    # Registered parcellation
dwi_to_T1w__dice_scores.csv           # QC DICE scores
dwi_to_T1w_Warped.nii.gz              # Warped output image
```

---

## Troubleshooting

### Problem: Test completes too fast (< 5 seconds)
**Cause:** Missing input files  
**Solution:** Check error message, provide missing files

### Problem: "ValueError: operands could not be broadcast (4,) (3,)"
**Cause:** 4D image used instead of 3D  
**Solution:** Extract 3D volume:
```bash
# Best: Use mrconvert (from MRtrix3 - included in LAMAReg conda env)
mrconvert -coord 3 0 dwi_fod.nii.gz dwi_fod_vol0.nii.gz

# Alternative: Use fslroi (FSL 6.0.3 on server)
source /data_/mica1/01_programs/fsl-6-0-3/etc/fslconf/fsl.sh
fslroi dwi_fod.nii.gz dwi_fod_vol0.nii.gz 0 1
```

### Problem: "lamareg command not found"
**Cause:** LAMAReg not in PATH or wrong conda env  
**Solution:** 
```bash
conda activate lamareg
which lamareg
```

### Problem: "Permission denied"
**Cause:** Test scripts not executable  
**Solution:**
```bash
chmod +x test_*.sh
```

---

## Testing Full Pipeline (After Tests Pass)

Once all unit tests pass, test with real data:

```bash
# Activate micapipe environment
source /path/to/micapipe/functions/init.sh

# Run DWI processing on a test subject
micapipe -sub test001 -ses 01 \
         -bids /path/to/bids \
         -out /path/to/derivatives \
         -proc_dwi
```

**Check the output:**
```bash
# Look for LAMAReg warpfields
ls derivatives/micapipe/sub-test001/ses-01/xfm/*LAMAReg*

# Check JSON metadata
cat derivatives/micapipe/sub-test001/ses-01/xfm/*transformations-proc_dwi*.json
# Should show: "registrationMethod": "LAMAReg"

# Verify DICE scores
cat derivatives/micapipe/sub-test001/ses-01/xfm/*dice_scores.csv
```

---

## Summary

**Before testing:**
1. ✅ Pull latest code
2. ✅ Activate LAMAReg conda environment
3. ✅ Extract 3D DWI reference (critical!)

**During testing:**
- Each test takes 10-15 minutes if running correctly
- Watch for "EXECUTING LAMAReg REGISTRATION" banner
- Check execution time is 900+ seconds

**After testing:**
- All 4 tests should pass
- 10 output files per test
- DICE scores CSV created
- Ready for full pipeline testing

---

## Key Files Reference

- `README.md` - Overview and setup
- `BUG_DWI_4D_FOD.md` - Known issue with 4D FOD
- `INSTRUCTIONS_FOR_SERVER.md` - This file (duplicate, can remove)
- `test_common.sh` - Shared test functions
- `test_dwi_registration.sh` - DWI test (standalone)
- `test_func_registration.sh` - Functional test
- `test_mpc_registration.sh` - MPC test
- `test_mpc_swm_registration.sh` - MPC-SWM test
- `setup_hcp_test_data.sh` - Data preparation script
- `run_all_tests.sh` - Run all tests sequentially

---

**Questions?** Check `BUG_DWI_4D_FOD.md` for the 4D FOD issue details.
