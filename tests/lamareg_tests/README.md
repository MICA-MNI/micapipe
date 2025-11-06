# LAMAReg Registration Unit Tests

This directory contains unit tests for validating LAMAReg registration functionality in micapipe.

## Test Structure

Each test file is designed to test a specific registration module independently:

1. **test_dwi_registration.sh** - DWI to T1w registration (02_proc-dwi.sh)
2. **test_func_registration.sh** - Functional MRI to T1nativepro registration (02_proc-func.sh)
3. **test_flair_registration.sh** - FLAIR to T1nativepro registration (02_proc-flair.sh)
4. **test_mpc_registration.sh** - Microstructural to FreeSurfer registration (03_MPC.sh)
5. **test_mpc_swm_registration.sh** - Microstructural to T1nativepro registration (03_MPC-SWM.sh)

## Running Tests

### Individual Test
```bash
cd tests/lamareg_tests
./test_dwi_registration.sh /path/to/test/data
```

### All Tests
```bash
cd tests/lamareg_tests
./run_all_tests.sh /path/to/test/data
```

### Validation Only (Check Outputs)
```bash
cd tests/lamareg_tests
./validate_outputs.sh /path/to/output/directory
```

## Test Data Requirements

Each test expects the following structure:
```
test_data/
├── sub-001/
│   ├── anat/
│   │   ├── sub-001_T1w.nii.gz
│   │   └── sub-001_FLAIR.nii.gz
│   ├── dwi/
│   │   ├── sub-001_dwi.nii.gz
│   │   └── sub-001_dwi.bval
│   │   └── sub-001_dwi.bvec
│   └── func/
│       └── sub-001_task-rest_bold.nii.gz
```

## What Each Test Validates

### 1. Command Syntax
- LAMAReg command is properly formatted
- All required arguments are present
- Named arguments match LAMAReg CLI specification

### 2. Output Files
- Primary warpfield: `*1Warp.nii.gz`
- Secondary warpfield: `*2Warp.nii.gz`
- Inverse warpfields: `*1InverseWarp.nii.gz`, `*2InverseWarp.nii.gz`
- Affine: `*0GenericAffine.mat`
- Parcellations: `*_fixed_parc.nii.gz`, `*_moving_parc.nii.gz`, `*_registered_parc.nii.gz`
- QC CSV: `*_dice_scores.csv`
- Warped output: `*Warped.nii.gz`

### 3. File Integrity
- All output files exist and are not empty
- NIfTI headers are valid
- Transformation matrices are properly formatted
- CSV files contain valid DICE scores

### 4. Registration Quality
- DICE scores are above acceptable thresholds
- Warped images show proper alignment
- No obvious artifacts or failures

### 5. Transformation Chains
- Forward and inverse transforms are correctly ordered
- ANTs transformation application works correctly
- Multiple transform composition works as expected

## Expected Output

Each test produces:
1. **test_results.txt** - Summary of pass/fail for each check
2. **test_log.txt** - Detailed execution log
3. **qc_report.html** - Visual QC report (if available)
4. **dice_scores.csv** - Aggregated DICE scores across all tests

## Success Criteria

A test passes if:
- ✅ All required output files are generated
- ✅ File sizes are reasonable (not empty, not suspiciously small)
- ✅ DICE scores > 0.7 for brain regions
- ✅ No error messages in logs
- ✅ Registered images show visual alignment

## Troubleshooting

### Test Fails: Missing Output Files
- Check that LAMAReg is installed: `which lamareg`
- Verify input data exists and is valid
- Review test logs for specific errors

### Test Fails: Low DICE Scores
- Check input image quality
- Verify modality contrast is sufficient
- Review LAMAReg verbose output

### Test Fails: Command Syntax Error
- Ensure LAMAReg version is compatible (>= 0.2.0)
- Check that all required arguments are provided
- Verify file paths are absolute and exist

## Test Maintenance

- Update tests when LAMAReg CLI changes
- Add new tests for new registration workflows
- Keep expected DICE thresholds up to date
- Document any known limitations or edge cases
