#!/bin/bash
#
# Common functions for LAMAReg registration tests
# This file should be sourced by individual test scripts
#

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

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
    
    ((TESTS_TOTAL++))
    
    if [ "$result" = "PASS" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $test_name"
        echo "✓ PASS: $test_name" >> "$RESULTS_FILE"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAIL${NC}: $test_name - $message"
        echo "✗ FAIL: $test_name - $message" >> "$RESULTS_FILE"
        ((TESTS_FAILED++))
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
    
    if [ ! -f "$moving_img" ]; then
        test_result "Input: Moving image" "FAIL" "File not found: $moving_img"
        echo -e "${YELLOW}Note: Place test data in expected location or skip execution test${NC}"
        return 1
    fi
    test_result "Input: Moving image" "PASS" ""
    
    if [ ! -f "$fixed_img" ]; then
        test_result "Input: Fixed image" "FAIL" "File not found: $fixed_img"
        return 1
    fi
    test_result "Input: Fixed image" "PASS" ""
    
    return 0
}

# Function to check required tools
check_required_tools() {
    # Check LAMAReg
    if ! command -v lamareg &> /dev/null; then
        test_result "LAMAReg installation" "FAIL" "lamareg command not found"
        echo -e "${RED}Error: lamareg is not installed or not in PATH${NC}"
        return 1
    fi
    test_result "LAMAReg installation" "PASS" ""
    
    # Check ANTs
    if ! command -v antsApplyTransforms &> /dev/null; then
        test_result "ANTs installation" "FAIL" "antsApplyTransforms not found"
        return 1
    fi
    test_result "ANTs installation" "PASS" ""
    
    return 0
}

# Function to validate all output files
validate_outputs() {
    local prefix="$1"
    
    setup_output_paths "$prefix"
    
    echo ""
    echo -e "${BLUE}Validating output files...${NC}"
    
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
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Test Summary${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo -e "Total tests: ${TESTS_TOTAL}"
    echo -e "${GREEN}Passed: ${TESTS_PASSED}${NC}"
    echo -e "${RED}Failed: ${TESTS_FAILED}${NC}"
    echo ""
    
    echo "" >> "$RESULTS_FILE"
    echo "========================================" >> "$RESULTS_FILE"
    echo "Summary:" >> "$RESULTS_FILE"
    echo "  Total tests: $TESTS_TOTAL" >> "$RESULTS_FILE"
    echo "  Passed: $TESTS_PASSED" >> "$RESULTS_FILE"
    echo "  Failed: $TESTS_FAILED" >> "$RESULTS_FILE"
    echo "" >> "$RESULTS_FILE"
    
    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}All tests passed!${NC}"
        echo "Result: ALL TESTS PASSED" >> "$RESULTS_FILE"
        return 0
    else
        echo -e "${RED}Some tests failed. Check $RESULTS_FILE for details.${NC}"
        echo "Result: SOME TESTS FAILED" >> "$RESULTS_FILE"
        return 1
    fi
}

# Main test framework
main() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}${TEST_NAME}${NC}"
    echo -e "${BLUE}========================================${NC}"
    
    # Parse arguments
    if [ $# -lt 1 ]; then
        echo -e "${RED}Error: Test data directory required${NC}"
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
    echo -e "${BLUE}Starting tests...${NC}"
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
