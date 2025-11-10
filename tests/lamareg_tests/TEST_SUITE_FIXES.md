# Test Suite Fix: Making Tests Actually Run LAMAReg

**Date:** November 10, 2025  
**Commits:** 75d072c (test_dwi_registration.sh), e9f67ce (test_common.sh)

## Problem Identified

You were absolutely right to be suspicious! The tests were completing **too fast** because they were:

1. **Skipping LAMAReg execution** when input files were missing
2. **Failing silently** without clear warnings
3. **Not showing** whether LAMAReg was actually running or just being skipped

### What Was Happening:

```bash
# Old behavior (silent skip):
validate_inputs()  # File not found
  → return 1
run_lamareg_test() || return 1  # Exits early, never runs LAMAReg
  → Test completes in 2 seconds
  → User thinks it ran but it didn't!
```

## Fixes Applied

### 1. test_dwi_registration.sh (Commit 75d072c)
- ✅ Removed `set -e` (no more silent exits)
- ✅ Removed color codes
- ✅ Replaced `echo -e` with plain `echo`
- ✅ Changed `((var++))` to `var=$((var + 1))`
- ✅ **Standalone test** - doesn't depend on test_common.sh

### 2. test_common.sh (Commit e9f67ce)
- ✅ Removed `set -e`
- ✅ Removed all color codes
- ✅ Replaced `echo -e` with plain `echo`
- ✅ Changed `((var++))` to `var=$((var + 1))`
- ✅ **Added verbose output** showing EXACTLY what's happening

### New Verbose Output:

```bash
==========================================
Checking input files...
==========================================
Moving image: /path/to/moving.nii.gz
Fixed image: /path/to/fixed.nii.gz

✗ FAIL: Input: Moving image - File not found: /path/to/moving.nii.gz
ERROR: Moving image file is missing!
SKIPPING LAMAReg execution - cannot proceed without input files
```

vs when files ARE present:

```bash
==========================================
Checking input files...
==========================================
Moving image: /data_/mica3/BIDS_CI/lamareg_test_data/dwi/T1w_in_dwi_brain.nii.gz
Fixed image: /data_/mica3/BIDS_CI/lamareg_test_data/dwi/dwi_fod.nii.gz

✓ PASS: Input: Moving image
✓ PASS: Input: Fixed image
✓ All input files present - proceeding with LAMAReg registration

==========================================
EXECUTING LAMAReg REGISTRATION
==========================================
[2025-11-10 13:45:23] Starting LAMAReg registration (this may take 10-15 minutes on CPU)...
[2025-11-10 13:45:23] Moving: T1w_in_dwi_brain.nii.gz
[2025-11-10 13:45:23] Fixed: dwi_fod.nii.gz
[2025-11-10 13:45:23] Output prefix: /output/dwi_to_T1w_

<... LAMAReg output for 10-15 minutes ...>

==========================================
[2025-11-10 14:00:45] LAMAReg registration completed successfully in 922 seconds
✓ PASS: LAMAReg execution
[2025-11-10 14:00:45] Output files: 10 created, 0 missing
==========================================
```

## What to Look For When Running Tests

### ✅ Test is ACTUALLY running LAMAReg:
```bash
==========================================
EXECUTING LAMAReg REGISTRATION
==========================================
```
- You'll see LAMAReg output scrolling for 10-15 minutes
- Execution time will be reported (e.g., "completed in 922 seconds")
- Output files will be counted (e.g., "10 created, 0 missing")

### ❌ Test is SKIPPING LAMAReg execution:
```bash
ERROR: Moving image file is missing!
SKIPPING LAMAReg execution - cannot proceed without input files
```
- Test completes in < 5 seconds
- No LAMAReg output shown
- Missing files clearly reported

## Required Input Files

### test_dwi_registration.sh
```bash
TEST_DATA_DIR/
  ├── T1w_in_dwi_brain.nii.gz    # Moving image
  └── dwi_fod.nii.gz              # Fixed image (FOD from DWI)
```

### test_func_registration.sh  
```bash
TEST_DATA_DIR/
  ├── func_brain.nii.gz           # Functional MRI (brain extracted)
  └── T1_nativepro.nii.gz         # T1 native processed
```

### test_mpc_registration.sh
```bash
TEST_DATA_DIR/
  ├── microstructural_map.nii.gz  # qMRI map
  └── T1_fsnative.nii.gz          # FreeSurfer native T1
```

### test_mpc_swm_registration.sh
```bash
TEST_DATA_DIR/
  ├── microstructural_map.nii.gz  # qMRI map
  └── T1_nativepro.nii.gz         # T1 native processed
```

## How to Verify Tests Are Working

### 1. Pull Latest Changes
```bash
cd /Users/enningyang/CodeProj/micapipe
git pull origin 156-update-regsynth-to-lamareg
```

### 2. Run Individual Test
```bash
cd tests/lamareg_tests

# DWI test (standalone, uses own implementation)
./test_dwi_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/dwi

# Functional test (uses test_common.sh)
./test_func_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/func

# MPC test (uses test_common.sh)
./test_mpc_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/mpc

# MPC-SWM test (uses test_common.sh)
./test_mpc_swm_registration.sh /data_/mica3/BIDS_CI/lamareg_test_data/mpc_swm
```

### 3. Check Output Verbosity

**If LAMAReg is running:**
- You'll see "EXECUTING LAMAReg REGISTRATION" banner
- LAMAReg output for ~10-15 minutes
- Execution time displayed (e.g., 900+ seconds)
- File counts shown

**If files are missing:**
- You'll see "ERROR: ... file is missing!"
- Clear "SKIPPING LAMAReg execution" message
- Test completes in seconds

### 4. Check Log Files
```bash
# Each test creates logs in output directory:
cat test_output_*/test_log.txt      # Full LAMAReg output
cat test_output_*/test_results.txt  # Pass/fail summary
```

## Expected Behavior Now

### Scenario 1: Input Files Present ✅
```
Total execution time: 10-15 minutes (actual LAMAReg running)
You'll see:
  - Tool version checks
  - Input validation (PASS)
  - "EXECUTING LAMAReg REGISTRATION" banner
  - Live LAMAReg output
  - Execution time (900+ seconds)
  - Output file verification
  - DICE scores if available
```

### Scenario 2: Input Files Missing ❌
```
Total execution time: < 5 seconds (LAMAReg skipped)
You'll see:
  - Tool version checks
  - Input validation (FAIL)
  - Clear error message
  - "SKIPPING LAMAReg execution" warning
  - Test summary (with failures for missing files)
```

## Testing Checklist

- [ ] Pull latest code (commits 75d072c, e9f67ce)
- [ ] Verify input files exist in test data directories
- [ ] Run test_dwi_registration.sh first (standalone, most stable)
- [ ] Confirm you see "EXECUTING LAMAReg REGISTRATION" banner
- [ ] Confirm LAMAReg actually runs (takes 10+ minutes, shows output)
- [ ] Check execution time is reported (should be 900+ seconds)
- [ ] Verify output files are created and counted
- [ ] Run other tests (func, mpc, mpc_swm)
- [ ] Compare execution times to ensure they're all actually running

## Summary

**Before:** Tests were silently skipping LAMAReg execution when files were missing, completing in seconds but appearing to pass.

**After:** Tests now clearly show:
1. ✅ When LAMAReg is ACTUALLY running (with "EXECUTING" banner and live output)
2. ❌ When LAMAReg is being skipped (with clear error messages)
3. ⏱️ How long LAMAReg took to run (execution time in seconds)
4. 📊 How many output files were created

**You can now trust that if a test completes in < 5 seconds, it's NOT running LAMAReg!**
