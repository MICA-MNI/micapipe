#!/bin/bash
#
# Common functions for LAMAReg registration tests
# This file should be sourced by individual test scripts
#

# NO set -e to avoid silent exits

# Test counter
TESTS_PASSED=0
TESTS_FAILED=0
TESTS_TOTAL=0

# Required output files
REQUIRED_OUTPUTS=(
    "0GenericAffine.mat"
    "1Warp.nii.gz"
    "1InverseWarp.nii.gz"
    "2Warp.nii.gz"
    "2InverseWarp.nii.gz"
    "_fixed_parc.nii.gz"
    "_moving_parc.nii.gz"
    "_registered_parc.nii.gz"
    "_dice_scores.csv"
    "Warped.nii.gz"
)

# Default thresholds
MIN_DICE_GLOBAL=${MIN_DICE_GLOBAL:-0.70}
MIN_DICE_GM=${MIN_DICE_GM:-0.65}
MIN_DICE_WM=${MIN_DICE_WM:-0.75}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to test result
test_result() {
    local test_name="$1"
    local result="$2"
    local message="$3"
    
    TESTS_TOTAL=$((TESTS_TOTAL + 1))
    
    if [ "$result" = "PASS" ]; then
        echo "✓ PASS: $test_name"
        echo "✓ PASS: $test_name" >> "$RESULTS_FILE"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo "✗ FAIL: $test_name - $message"
        echo "✗ FAIL: $test_name - $message" >> "$RESULTS_FILE"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

# Function to check file exists and is not empty
check_file() {
    local filepath="$1"
    local filename=$(basename "$filepath")
    
    if [ ! -f "$filepath" ]; then
        test_result "File exists: $filename" "FAIL" "File not found"
        return 1
    fi
    
    local filesize=$(stat -f%z "$filepath" 2>/dev/null || stat -c%s "$filepath" 2>/dev/null)
    if [ "$filesize" -lt 100 ]; then
        test_result "File size: $filename" "FAIL" "File too small ($filesize bytes)"
        return 1
    fi
    
    test_result "File exists: $filename" "PASS" ""
    return 0
}

# Function to validate NIfTI file
validate_nifti() {
    local filepath="$1"
    local filename=$(basename "$filepath")
    
    if command -v fslinfo &> /dev/null; then
        if fslinfo "$filepath" &> /dev/null; then
            test_result "NIfTI valid: $filename" "PASS" ""
            return 0
        else
            test_result "NIfTI valid: $filename" "FAIL" "Invalid NIfTI header"
            return 1
        fi
    else
        log "Warning: fslinfo not available, skipping NIfTI validation"
        return 0
    fi
}

# Function to check DICE scores
check_dice_scores() {
    local csv_file="$1"
    
    if [ ! -f "$csv_file" ]; then
        test_result "DICE scores CSV" "FAIL" "File not found"
        return 1
    fi
    
    # Check if CSV has content
    local line_count=$(wc -l < "$csv_file")
    if [ "$line_count" -lt 2 ]; then
        test_result "DICE scores content" "FAIL" "CSV has insufficient data"
        return 1
    fi
    
    test_result "DICE scores CSV" "PASS" ""
    
    # Parse and validate DICE scores
    if command -v python3 &> /dev/null; then
        local avg_dice=$(python3 -c "
import csv
import sys
try:
    with open('$csv_file', 'r') as f:
        reader = csv.DictReader(f)
        scores = [float(row['dice']) for row in reader if 'dice' in row]
        if scores:
            print(sum(scores) / len(scores))
        else:
            print('0')
except:
    print('0')
" 2>/dev/null)
        
        if (( $(echo "$avg_dice >= $MIN_DICE_GLOBAL" | bc -l) )); then
            test_result "DICE scores quality (avg=$avg_dice)" "PASS" ""
        else
            test_result "DICE scores quality (avg=$avg_dice)" "FAIL" "Below threshold ($MIN_DICE_GLOBAL)"
        fi
    fi
}

# Function to setup output paths
setup_output_paths() {
    local prefix="$1"
    
    AFFINE="${prefix}0GenericAffine.mat"
    WARP1="${prefix}1Warp.nii.gz"
    INVWARP1="${prefix}1InverseWarp.nii.gz"
    WARP2="${prefix}2Warp.nii.gz"
    INVWARP2="${prefix}2InverseWarp.nii.gz"
    MOVING_PARC="${prefix}_moving_parc.nii.gz"
    FIXED_PARC="${prefix}_fixed_parc.nii.gz"
    REG_PARC="${prefix}_registered_parc.nii.gz"
    QC_CSV="${prefix}_dice_scores.csv"
    WARPED="${prefix}Warped.nii.gz"
}

# Function to validate inputs
validate_inputs() {
    local moving_img="$1"
    local fixed_img="$2"
    
    echo ""
    echo "=========================================="
    echo "Checking input files..."
    echo "=========================================="
    echo "Moving image: $moving_img"
    echo "Fixed image: $fixed_img"
    echo ""
    
    if [ ! -f "$moving_img" ]; then
        test_result "Input: Moving image" "FAIL" "File not found: $moving_img"
        echo "ERROR: Moving image file is missing!"
        echo "SKIPPING LAMAReg execution - cannot proceed without input files"
        echo ""
        return 1
    fi
    test_result "Input: Moving image" "PASS" ""
    
    if [ ! -f "$fixed_img" ]; then
        test_result "Input: Fixed image" "FAIL" "File not found: $fixed_img"
        echo "ERROR: Fixed image file is missing!"
        echo "SKIPPING LAMAReg execution - cannot proceed without input files"
        echo ""
        return 1
    fi
    test_result "Input: Fixed image" "PASS" ""
    
    echo "✓ All input files present - proceeding with LAMAReg registration"
    echo ""
    return 0
}

# Function to check required tools
check_required_tools() {
    echo ""
    echo "=========================================="
    echo "Checking required tools..."
    echo "=========================================="
    
    # Check LAMAReg
    if ! command -v lamareg &> /dev/null; then
        test_result "LAMAReg installation" "FAIL" "lamareg command not found"
        echo "ERROR: lamareg is not installed or not in PATH"
        echo "Install LAMAReg or add it to your PATH before running tests"
        return 1
    fi
    
    # Get LAMAReg version
    local lamareg_version=$(lamareg --version 2>&1 | head -1 || echo "unknown")
    test_result "LAMAReg installation" "PASS" ""
    echo "LAMAReg version: $lamareg_version"
    
    # Check ANTs
    if ! command -v antsApplyTransforms &> /dev/null; then
        test_result "ANTs installation" "FAIL" "antsApplyTransforms not found"
        return 1
    fi
    test_result "ANTs installation" "PASS" ""
    
    echo "✓ All required tools are available"
    echo ""
    return 0
}

# Function to execute LAMAReg registration
execute_lamareg() {
    local moving_img="$1"
    local fixed_img="$2"
    local output_prefix="$3"
    
    echo ""
    echo "=========================================="
    echo "EXECUTING LAMAReg REGISTRATION"
    echo "=========================================="
    log "Starting LAMAReg registration (this may take 10-15 minutes on CPU)..."
    log "Moving: $(basename $moving_img)"
    log "Fixed: $(basename $fixed_img)"
    log "Output prefix: $output_prefix"
    echo ""
    
    local start_time=$(date +%s)
    
    lamareg register \
      --moving "$moving_img" \
      --fixed "$fixed_img" \
      --output "${WARPED}" \
      --moving-parc "${MOVING_PARC}" \
      --fixed-parc "${FIXED_PARC}" \
      --registered-parc "${REG_PARC}" \
      --affine "${AFFINE}" \
      --warpfield "${WARP1}" \
      --inverse-warpfield "${INVWARP1}" \
      --secondary-warpfield "${WARP2}" \
      --inverse-secondary-warpfield "${INVWARP2}" \
      --qc-csv "${QC_CSV}" \
      --synthseg-threads 4 \
      --ants-threads 8 2>&1 | tee -a "$LOG_FILE"
    
    local exit_code=${PIPESTATUS[0]}
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    echo ""
    echo "=========================================="
    if [ $exit_code -eq 0 ]; then
        log "LAMAReg registration completed successfully in ${duration} seconds"
        test_result "LAMAReg execution" "PASS" ""
        
        # Check output files were created
        local files_created=0
        local files_missing=0
        for file in "$WARPED" "$AFFINE" "$WARP1" "$INVWARP1" "$WARP2" "$INVWARP2" "$MOVING_PARC" "$FIXED_PARC" "$REG_PARC" "$QC_CSV"; do
            if [ -f "$file" ]; then
                files_created=$((files_created + 1))
            else
                files_missing=$((files_missing + 1))
                log "WARNING: Expected output file not created: $(basename $file)"
            fi
        done
        
        log "Output files: $files_created created, $files_missing missing"
        echo "=========================================="
        echo ""
        return 0
    else
        log "LAMAReg registration FAILED with exit code $exit_code after ${duration} seconds"
        test_result "LAMAReg execution" "FAIL" "Registration failed (exit code: $exit_code)"
        echo "=========================================="
        echo ""
        return 1
    fi
}

# Function to setup output paths for LAMAReg
setup_output_paths() {
    local prefix="$1"
    AFFINE="${prefix}0GenericAffine.mat"
    WARP1="${prefix}1Warp.nii.gz"
    INVWARP1="${prefix}1InverseWarp.nii.gz"
    WARP2="${prefix}2Warp.nii.gz"
    INVWARP2="${prefix}2InverseWarp.nii.gz"
    MOVING_PARC="${prefix}_moving_parc.nii.gz"
    FIXED_PARC="${prefix}_fixed_parc.nii.gz"
    REG_PARC="${prefix}_registered_parc.nii.gz"
    QC_CSV="${prefix}_dice_scores.csv"
    WARPED="${prefix}Warped.nii.gz"
}

# Function to validate all output files
validate_outputs() {
    local prefix="$1"
    
    setup_output_paths "$prefix"
    
    echo ""
    echo "Validating output files..."
    
    for suffix in "${REQUIRED_OUTPUTS[@]}"; do
        filepath="${prefix}${suffix}"
        if [ -f "$filepath" ]; then
            check_file "$filepath"
            
            # Validate NIfTI files
            if [[ "$filepath" == *.nii.gz ]]; then
                validate_nifti "$filepath"
            fi
        fi
    done
    
    # Check DICE scores
    if [ -f "$QC_CSV" ]; then
        check_dice_scores "$QC_CSV"
    fi
}

# Function to print test summary
print_summary() {
    echo ""
    echo "========================================"
    echo "Test Summary"
    echo "========================================"
    echo "Total tests: ${TESTS_TOTAL}"
    echo "Passed: ${TESTS_PASSED}"
    echo "Failed: ${TESTS_FAILED}"
    echo ""
    
    echo "" >> "$RESULTS_FILE"
    echo "========================================" >> "$RESULTS_FILE"
    echo "Summary:" >> "$RESULTS_FILE"
    echo "  Total tests: $TESTS_TOTAL" >> "$RESULTS_FILE"
    echo "  Passed: $TESTS_PASSED" >> "$RESULTS_FILE"
    echo "  Failed: $TESTS_FAILED" >> "$RESULTS_FILE"
    echo "" >> "$RESULTS_FILE"
    
    if [ $TESTS_FAILED -eq 0 ]; then
        echo "All tests passed!"
        echo "Result: ALL TESTS PASSED" >> "$RESULTS_FILE"
        return 0
    else
        echo "Some tests failed. Check $RESULTS_FILE for details."
        echo "Result: SOME TESTS FAILED" >> "$RESULTS_FILE"
        return 1
    fi
}

# Main test framework
main() {
    echo "========================================"
    echo "$TEST_NAME"
    echo "========================================"
    
    # Parse arguments
    if [ $# -lt 1 ]; then
        echo "Error: Test data directory required"
        echo "Usage: $0 <test_data_dir> [output_dir]"
        exit 1
    fi
    
    TEST_DATA_DIR="$1"
    OUTPUT_DIR="${2:-$(pwd)/test_output_$(basename "$0" .sh)}"
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    LOG_FILE="$OUTPUT_DIR/test_log.txt"
    RESULTS_FILE="$OUTPUT_DIR/test_results.txt"
    
    # Initialize results file
    echo "$TEST_NAME" > "$RESULTS_FILE"
    echo "Test Date: $(date)" >> "$RESULTS_FILE"
    echo "Module: $TEST_MODULE" >> "$RESULTS_FILE"
    echo "========================================" >> "$RESULTS_FILE"
    echo "" >> "$RESULTS_FILE"
    
    echo ""
    echo "Starting tests..."
    echo ""
    
    # Check required tools
    check_required_tools
    
    # Run test-specific function
    if type run_lamareg_test &>/dev/null; then
        run_lamareg_test
    fi
    
    # Validate outputs if they exist
    if [ -n "$OUTPUT_PREFIX" ]; then
        validate_outputs "$OUTPUT_DIR/$OUTPUT_PREFIX"
    fi
    
    # Print summary
    print_summary
    
    # Exit with appropriate code
    [ $TESTS_FAILED -eq 0 ] && exit 0 || exit 1
}
