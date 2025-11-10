# LAMAReg Test Framework Status

**Date:** November 10, 2025  
**Branch:** 156-update-regsynth-to-lamareg

## Test Files Overview

### ✅ test_dwi_registration.sh
**Status:** FIXED and SIMPLIFIED  
**Last Commit:** 75d072c

**Changes Made:**
- ✅ Removed `set -e` to prevent silent exits
- ✅ Removed all color codes (RED, GREEN, BLUE, NC)
- ✅ Replaced `echo -e` with plain `echo`
- ✅ Changed `((var++))` to `var=$((var + 1))`
- ✅ LAMAReg command properly implemented (direct call, not eval)
- ✅ Tested locally with fake data before pushing

**LAMAReg Integration:**
```bash
lamareg register \
  --moving "$moving_img" \
  --fixed "$fixed_img" \
  --output "$warped" \
  --moving-parc "$moving_parc" \
  --fixed-parc "$fixed_parc" \
  --registered-parc "$reg_parc" \
  --affine "$affine" \
  --warpfield "$warp1" \
  --inverse-warpfield "$invwarp1" \
  --secondary-warpfield "$warp2" \
  --inverse-secondary-warpfield "$invwarp2" \
  --qc-csv "$qc_csv" \
  --synthseg-threads 4 \
  --ants-threads 8 2>&1 | tee -a "$LOG_FILE"
```

---

### ⚠️ test_func_registration.sh
**Status:** USES test_common.sh (potential issues)  
**LAMAReg Integration:** Via `execute_lamareg()` function

**Structure:**
- Sources test_common.sh for shared functions
- Calls `execute_lamareg("func_brain.nii.gz", "T1_nativepro.nii.gz", "func_to_T1nativepro_")`
- Uses shared `test_result()` and `log()` functions

**Potential Issues:**
- Depends on test_common.sh which still has:
  - Color codes with `echo -e`
  - `((var++))` arithmetic
  - May have same silent exit issues

---

### ⚠️ test_mpc_registration.sh
**Status:** USES test_common.sh (potential issues)  
**LAMAReg Integration:** Via `execute_lamareg()` function

**Structure:**
- Sources test_common.sh
- Calls `execute_lamareg("microstructural_map.nii.gz", "T1_fsnative.nii.gz", "qMRI_to_fsnative_")`

**Potential Issues:**
- Same as test_func_registration.sh (depends on test_common.sh)

---

### ⚠️ test_mpc_swm_registration.sh
**Status:** USES test_common.sh (potential issues)  
**LAMAReg Integration:** Via `execute_lamareg()` function

**Structure:**
- Sources test_common.sh
- Calls `execute_lamareg("microstructural_map.nii.gz", "T1_nativepro.nii.gz", "qMRI_to_nativepro_")`

**Potential Issues:**
- Same as above (depends on test_common.sh)

---

## test_common.sh Analysis

**Current Issues:**
```bash
# Line 7-11: Color codes still defined
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Line 49-59: Uses echo -e with colors
test_result() {
    ((TESTS_TOTAL++))  # May fail on some shells
    
    if [ "$result" = "PASS" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $test_name"  # echo -e with color codes
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAIL${NC}: $test_name - $message"  # echo -e with color codes
        ((TESTS_FAILED++))
    fi
}
```

**execute_lamareg() Function (Lines 202-234):**
✅ LAMAReg command is properly implemented
✅ Uses direct call (not eval)
✅ Proper error handling with exit code checking
✅ All required parameters included

---

## Recommendations

### Option 1: Fix test_common.sh (Recommended if other tests need to work)
Apply the same fixes we made to test_dwi_registration.sh:
1. Remove color code variables
2. Replace all `echo -e` with plain `echo`
3. Replace `((var++))` with `var=$((var + 1))`
4. Test locally before pushing

### Option 2: Keep test_dwi_registration.sh independent (Current approach)
- test_dwi_registration.sh is now standalone and working
- Other tests may still have issues due to test_common.sh

---

## LAMAReg Command Verification

Both implementations use the correct LAMAReg syntax:

**Required Parameters:**
- ✅ `--moving` (source image)
- ✅ `--fixed` (target image)
- ✅ `--output` (warped output)
- ✅ `--affine` (affine transformation)
- ✅ `--warpfield` (primary warpfield)
- ✅ `--inverse-warpfield` (inverse primary warp)
- ✅ `--secondary-warpfield` (secondary warpfield)
- ✅ `--inverse-secondary-warpfield` (inverse secondary warp)
- ✅ `--moving-parc` (moving parcellation for QC)
- ✅ `--fixed-parc` (fixed parcellation for QC)
- ✅ `--registered-parc` (output registered parcellation)
- ✅ `--qc-csv` (DICE scores CSV output)
- ✅ `--synthseg-threads` (4 threads)
- ✅ `--ants-threads` (8 threads)

**Output capture:**
- ✅ Both implementations use `2>&1 | tee -a "$LOG_FILE"` to capture output
- ✅ Both check `${PIPESTATUS[0]}` for actual lamareg exit code

---

## Next Steps

1. **If other tests are failing:** Fix test_common.sh with same approach as test_dwi_registration.sh
2. **If other tests are working:** The color code issues may be less severe on the server environment
3. **Test on server:** Pull latest changes and run all tests with real data
4. **Monitor execution:** Verify LAMAReg actually runs (not just prints start message)

---

## Summary

✅ **LAMAReg Integration:** Command syntax is correct in all tests  
⚠️ **Shell Compatibility:** test_common.sh may have same issues as old test_dwi_registration.sh  
✅ **test_dwi_registration.sh:** Fixed and verified locally  
❓ **Other tests:** Need server testing to confirm they work with test_common.sh

**The LAMAReg replacement itself is correctly implemented** - the issues were shell script compatibility, not LAMAReg integration problems.
