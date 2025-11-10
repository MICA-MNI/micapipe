#!/bin/bash
#
# Simplified LAMAReg DWI Registration Test
# Tests: 02_proc-dwi.sh - T1w to DWI registration
#
# Usage: ./test_dwi_registration.sh <test_data_dir> [output_dir]
#

# NO set -e to avoid silent exits

# Test configuration
TEST_NAME="DWI Registration Test"
TEST_MODULE="02_proc-dwi.sh"

# Counters
TESTS_TOTAL=0
TESTS_PASSED=0
TESTS_FAILED=0

# Get directories
TEST_DATA_DIR="${1:-}"
OUTPUT_DIR="${2:-./test_output}"

if [ -z "$TEST_DATA_DIR" ]; then
    echo "Usage: $0 <test_data_dir> [output_dir]"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Set up log files
LOG_FILE="$OUTPUT_DIR/test.log"
RESULTS_FILE="$OUTPUT_DIR/results.txt"

# Clear previous results
> "$LOG_FILE"
> "$RESULTS_FILE"

echo "========================================"
echo "$TEST_NAME"
echo "========================================"
echo "Test data: $TEST_DATA_DIR"
echo "Output: $OUTPUT_DIR"
echo "Log: $LOG_FILE"
echo ""

# Simple log function
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Simple test_result function
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

# Function to simulate LAMAReg registration test
test_lamareg_registration() {
    log "Starting LAMAReg registration test..."
    
    # Check if lamareg is available
    if ! command -v lamareg &> /dev/null; then
        test_result "LAMAReg installation" "FAIL" "lamareg command not found"
        echo "Error: lamareg is not installed or not in PATH"
        return 1
    fi
    test_result "LAMAReg installation" "PASS" ""
    
    # Define test file paths
    local moving_img="$TEST_DATA_DIR/T1w_in_dwi_brain.nii.gz"
    local fixed_img="$TEST_DATA_DIR/dwi_fod.nii.gz"
    local output_prefix="$OUTPUT_DIR/dwi_to_T1w_"
    
    # Check input files exist
    if [ ! -f "$moving_img" ]; then
        test_result "Input: Moving image" "FAIL" "File not found: $moving_img"
        log "File not found: $moving_img"
        return 1
    fi
    test_result "Input: Moving image" "PASS" ""
    
    if [ ! -f "$fixed_img" ]; then
        test_result "Input: Fixed image" "FAIL" "File not found: $fixed_img"
        return 1
    fi
    test_result "Input: Fixed image" "PASS" ""
    
    # Define output paths
    local affine="${output_prefix}0GenericAffine.mat"
    local warp1="${output_prefix}1Warp.nii.gz"
    local invwarp1="${output_prefix}1InverseWarp.nii.gz"
    local warp2="${output_prefix}2Warp.nii.gz"
    local invwarp2="${output_prefix}2InverseWarp.nii.gz"
    local moving_parc="${output_prefix}_moving_parc.nii.gz"
    local fixed_parc="${output_prefix}_fixed_parc.nii.gz"
    local reg_parc="${output_prefix}_registered_parc.nii.gz"
    local qc_csv="${output_prefix}_dice_scores.csv"
    local warped="${output_prefix}Warped.nii.gz"
    
    # Run LAMAReg registration
    log "Executing LAMAReg registration (this may take 10-15 minutes on CPU)..."
    log "Moving: $(basename $moving_img)"
    log "Fixed: $(basename $fixed_img)"
    
    # Execute LAMAReg directly (not through eval to avoid quoting issues)
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
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ]; then
        log "LAMAReg registration completed successfully"
        test_result "LAMAReg execution" "PASS" ""
    else
        log "LAMAReg registration failed with exit code $exit_code"
        test_result "LAMAReg execution" "FAIL" "Registration failed (exit code: $exit_code)"
        return 1
    fi
    
    # Check output files
    for file in "$warped" "$affine" "$warp1" "$invwarp1" "$moving_parc" "$fixed_parc" "$reg_parc" "$qc_csv"; do
        if [ -f "$file" ]; then
            test_result "Output file: $(basename $file)" "PASS" ""
        else
            test_result "Output file: $(basename $file)" "FAIL" "File not created"
        fi
    done
    
    return 0
}

echo ""
echo "Starting tests..."
echo ""

# Run tests
test_lamareg_registration

# Print summary
echo ""
echo "========================================"
echo "Test Summary"
echo "========================================"
echo "Total tests: $TESTS_TOTAL"
echo "Passed: $TESTS_PASSED"
echo "Failed: $TESTS_FAILED"
echo ""

if [ $TESTS_FAILED -eq 0 ]; then
    echo "All tests PASSED!"
    exit 0
else
    echo "Some tests FAILED. Check $RESULTS_FILE for details."
    exit 1
fi
