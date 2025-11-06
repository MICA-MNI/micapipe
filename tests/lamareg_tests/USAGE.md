# LAMAReg Registration Test Suite

## Quick Start

```bash
# Run all tests
cd tests/lamareg_tests
./run_all_tests.sh /path/to/test/data

# Run individual test
./test_dwi_registration.sh /path/to/test/data

# Validate existing outputs
./validate_outputs.sh /path/to/registration/output
```

## Test Files

| Test Script | Module Tested | Purpose |
|------------|---------------|---------|
| `test_dwi_registration.sh` | `02_proc-dwi.sh` | DWI → T1w registration |
| `test_func_registration.sh` | `02_proc-func.sh` | fMRI → T1nativepro registration |
| `test_flair_registration.sh` | `02_proc-flair.sh` | FLAIR → T1nativepro registration |
| `test_mpc_registration.sh` | `03_MPC.sh` | Microstructural → FreeSurfer registration |
| `test_mpc_swm_registration.sh` | `03_MPC-SWM.sh` | Microstructural → T1nativepro registration |
| `run_all_tests.sh` | All | Master test runner |
| `validate_outputs.sh` | N/A | Output validation utility |
| `test_common.sh` | N/A | Shared test functions |

## What Gets Tested

### 1. Installation & Prerequisites
- ✅ LAMAReg CLI is installed and accessible
- ✅ ANTs tools are available
- ✅ Required dependencies (FSL, Python) present

### 2. Command Syntax Validation
- ✅ LAMAReg commands use correct syntax
- ✅ All required arguments are present
- ✅ Named arguments match LAMAReg specification
- ✅ Secondary warpfield arguments included

### 3. Output File Generation
Checks for presence and validity of:
- `*0GenericAffine.mat` - Affine transformation
- `*1Warp.nii.gz` - Primary warpfield
- `*1InverseWarp.nii.gz` - Primary inverse warpfield
- `*2Warp.nii.gz` - **Secondary warpfield** (robust mode)
- `*2InverseWarp.nii.gz` - **Secondary inverse warpfield**
- `*_fixed_parc.nii.gz` - Fixed image parcellation
- `*_moving_parc.nii.gz` - Moving image parcellation
- `*_registered_parc.nii.gz` - Registered parcellation
- `*_dice_scores.csv` - QC metrics
- `*Warped.nii.gz` - Final registered image

### 4. File Integrity
- ✅ Files exist and are not empty (>100 bytes)
- ✅ NIfTI headers are valid (if FSL available)
- ✅ CSV files are properly formatted
- ✅ Transformation matrices are valid

### 5. Registration Quality (if outputs exist)
- ✅ DICE scores meet minimum thresholds:
  - Global: ≥ 0.70
  - Gray matter: ≥ 0.65
  - White matter: ≥ 0.75
- ✅ No obvious registration failures

### 6. Transformation Chain Validation
- ✅ Forward transforms: `-t warp2 -t warp1 -t affine`
- ✅ Inverse transforms: `-t affine -t invwarp1 -t invwarp2`
- ✅ ANTs can apply transformations

## Test Output Structure

```
test_results_all/
├── master_test_log.txt          # Complete execution log
├── test_summary.txt             # High-level summary
├── test_dwi_registration/
│   ├── test_log.txt
│   ├── test_results.txt
│   └── [registration outputs if run]
├── test_func_registration/
│   └── ...
├── test_flair_registration/
│   └── ...
├── test_mpc_registration/
│   └── ...
└── test_mpc_swm_registration/
    └── ...
```

## Success Criteria

A test **passes** if:
1. ✅ All required tools are installed
2. ✅ Command syntax is validated
3. ✅ (If data exists) All output files generated
4. ✅ (If data exists) DICE scores above thresholds
5. ✅ No errors in execution logs

## Usage Examples

### Example 1: Run Full Test Suite
```bash
# With test data
./run_all_tests.sh /data/micapipe/test_subjects ./results

# Syntax validation only (no data)
./run_all_tests.sh /tmp/dummy_dir
```

### Example 2: Test Specific Module
```bash
# Test DWI registration
./test_dwi_registration.sh /data/test_subject_001/dwi ./dwi_results

# View results
cat ./dwi_results/test_results.txt
```

### Example 3: Validate Existing Pipeline Output
```bash
# Check outputs from real pipeline run
./validate_outputs.sh /data/derivatives/micapipe/sub-001/xfm

# Generate validation report
# Creates: validation_report.txt with file checks and DICE scores
```

### Example 4: Quick Syntax Check (No Data Needed)
```bash
# Tests will validate command syntax even without real data
./test_dwi_registration.sh /tmp/fake_dir

# Output shows:
#  ✓ PASS: LAMAReg installation
#  ✓ PASS: LAMAReg command syntax
#  ✗ FAIL: Input files (expected, no data provided)
```

## Test Data Requirements

### Minimal Test Data (for syntax validation)
- No actual data needed
- Tests validate command structure only

### Full Test Data (for execution testing)
Each test expects specific files:

**DWI Test:**
```
test_data/
├── T1w_in_dwi_brain.nii.gz
└── dwi_fod.nii.gz
```

**Functional Test:**
```
test_data/
├── func_brain.nii.gz
└── T1_nativepro.nii.gz
```

**FLAIR Test:**
```
test_data/
├── FLAIR_preproc.nii.gz
└── T1_nativepro.nii.gz
```

**MPC Tests:**
```
test_data/
├── microstructural_map.nii.gz
├── T1_fsnative.nii.gz        # For MPC
└── T1_nativepro.nii.gz       # For MPC-SWM
```

## Interpreting Results

### All Tests Pass
```
✓ PASS: LAMAReg installation
✓ PASS: ANTs installation
✓ PASS: LAMAReg command syntax
✓ PASS: Transformation chain
...
========================================
Total tests: 15
Passed: 15
Failed: 0

✓ ALL TESTS PASSED!
```
**Action**: Ready for production use

### Syntax Validation Pass, No Data
```
✓ PASS: LAMAReg installation
✓ PASS: LAMAReg command syntax
✗ FAIL: Input: Moving image - File not found
...
Total tests: 10
Passed: 8
Failed: 2
```
**Action**: Normal for syntax-only testing. Commands are correct.

### Registration Quality Issues
```
✓ PASS: All output files generated
✗ FAIL: DICE scores quality (avg=0.58) - Below threshold (0.70)
```
**Action**: Check input data quality, review LAMAReg verbose output

## Troubleshooting

### "lamareg command not found"
```bash
# Install LAMAReg
pip install lamareg
# or
conda install -c conda-forge lamareg

# Verify
which lamareg
lamareg --help
```

### "antsApplyTransforms not found"
```bash
# Install ANTs
# macOS:
brew install ants
# Linux:
sudo apt-get install ants
```

### Low DICE Scores
- Check input image quality and contrast
- Verify images are from same subject
- Review LAMAReg verbose output for warnings
- Try adjusting thread counts

### Tests Take Too Long
- Reduce image resolution for testing
- Use subset of data
- Run individual tests instead of full suite

## Extending the Tests

### Add New Test
1. Copy existing test script
2. Update TEST_NAME, TEST_MODULE, and file paths
3. Implement `run_lamareg_test()` function
4. Add to `TESTS` array in `run_all_tests.sh`

### Modify Thresholds
Edit `MIN_DICE_*` variables in individual test scripts or `test_common.sh`

### Add Custom Validation
Extend `validate_outputs()` in `test_common.sh`

## CI/CD Integration

```yaml
# Example GitHub Actions workflow
- name: Run LAMAReg Tests
  run: |
    cd tests/lamareg_tests
    ./run_all_tests.sh $TEST_DATA_PATH
    
- name: Upload Test Results
  uses: actions/upload-artifact@v3
  with:
    name: test-results
    path: tests/lamareg_tests/test_results_all/
```

## Notes

- Tests designed to run with or without real data
- Syntax validation works without inputs
- Full validation requires actual test images
- All scripts are self-contained and portable
- DICE score thresholds are configurable
- Tests run independently (can run in parallel)

## Support

For issues with:
- **LAMAReg**: https://github.com/MICA-MNI/LAMAReg/issues
- **micapipe**: https://github.com/MICA-MNI/micapipe/issues
- **These tests**: Create issue in micapipe repo with "test" label
