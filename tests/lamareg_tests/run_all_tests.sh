#!/bin/bash
#
# Run All LAMAReg Registration Tests
# Usage: ./run_all_tests.sh <test_data_dir> [output_dir]
#

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
NC='\033[0m'

echo -e "${MAGENTA}========================================${NC}"
echo -e "${MAGENTA}  LAMAReg Registration Test Suite${NC}"
echo -e "${MAGENTA}========================================${NC}"
echo ""

# Parse arguments
if [ $# -lt 1 ]; then
    echo -e "${RED}Error: Test data directory required${NC}"
    echo "Usage: $0 <test_data_dir> [output_dir]"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/test/data ./test_results"
    exit 1
fi

TEST_DATA_DIR="$1"
OUTPUT_DIR="${2:-$(pwd)/test_results_all}"

# Create output directory
mkdir -p "$OUTPUT_DIR"
MASTER_LOG="$OUTPUT_DIR/master_test_log.txt"
SUMMARY_FILE="$OUTPUT_DIR/test_summary.txt"

echo "LAMAReg Registration Test Suite" > "$SUMMARY_FILE"
echo "Test Date: $(date)" >> "$SUMMARY_FILE"
echo "Test Data: $TEST_DATA_DIR" >> "$SUMMARY_FILE"
echo "Output Dir: $OUTPUT_DIR" >> "$SUMMARY_FILE"
echo "========================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# Test scripts and their corresponding data subdirectories
TESTS=(
    "test_dwi_registration.sh"
    "test_func_registration.sh"
    "test_flair_registration.sh"
    "test_mpc_registration.sh"
    "test_mpc_swm_registration.sh"
)

# Test names
TEST_NAMES=(
    "DWI Registration"
    "Functional MRI Registration"
    "FLAIR Registration"
    "MPC Registration"
    "MPC-SWM Registration"
)

# Data subdirectories for each test
TEST_SUBDIRS=(
    "dwi"
    "func"
    "flair"
    "mpc"
    "mpc"
)

# Results tracking
TOTAL_TESTS=${#TESTS[@]}
PASSED_TESTS=0
FAILED_TESTS=0
SKIPPED_TESTS=0

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Run each test
for i in "${!TESTS[@]}"; do
    test_script="${TESTS[$i]}"
    test_name="${TEST_NAMES[$i]}"
    test_subdir="${TEST_SUBDIRS[$i]}"
    test_data_path="$TEST_DATA_DIR/$test_subdir"
    test_output_dir="$OUTPUT_DIR/$(basename "$test_script" .sh)"
    
    # Check if test data directory exists
    if [ ! -d "$test_data_path" ]; then
        echo ""
        echo -e "${YELLOW}========================================${NC}"
        echo -e "${YELLOW}Skipping: $test_name${NC}"
        echo -e "${YELLOW}========================================${NC}"
        echo -e "${YELLOW}Data directory not found: $test_data_path${NC}"
        echo ""
        echo "⊘ SKIPPED: $test_name (no data)" >> "$SUMMARY_FILE"
        ((SKIPPED_TESTS++))
        continue
    fi
    
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Running: $test_name${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
    
    # Make script executable
    chmod +x "$SCRIPT_DIR/$test_script"
    
    # Run test (disable exit on error temporarily to capture exit code)
    set +e
    "$SCRIPT_DIR/$test_script" "$test_data_path" "$test_output_dir" 2>&1 | tee -a "$MASTER_LOG"
    test_exit_code=$?
    set -e
    
    # Check result
    if [ $test_exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ PASSED: $test_name${NC}"
        echo "✓ PASSED: $test_name" >> "$SUMMARY_FILE"
        ((PASSED_TESTS++))
    else
        echo -e "${RED}✗ FAILED: $test_name${NC}"
        echo "✗ FAILED: $test_name" >> "$SUMMARY_FILE"
        ((FAILED_TESTS++))
    fi
done

# Final summary
echo ""
echo -e "${MAGENTA}========================================${NC}"
echo -e "${MAGENTA}  Test Suite Summary${NC}"
echo -e "${MAGENTA}========================================${NC}"
echo -e "Total tests: ${TOTAL_TESTS}"
echo -e "${GREEN}Passed: ${PASSED_TESTS}${NC}"
echo -e "${RED}Failed: ${FAILED_TESTS}${NC}"
echo -e "${YELLOW}Skipped: ${SKIPPED_TESTS}${NC}"
echo ""
echo -e "Detailed results: ${OUTPUT_DIR}"
echo -e "Master log: ${MASTER_LOG}"
echo -e "Summary: ${SUMMARY_FILE}"
echo ""

# Write final summary
echo "" >> "$SUMMARY_FILE"
echo "========================================" >> "$SUMMARY_FILE"
echo "Final Summary:" >> "$SUMMARY_FILE"
echo "  Total tests: $TOTAL_TESTS" >> "$SUMMARY_FILE"
echo "  Passed: $PASSED_TESTS" >> "$SUMMARY_FILE"
echo "  Failed: $FAILED_TESTS" >> "$SUMMARY_FILE"
echo "  Skipped: $SKIPPED_TESTS" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}════════════════════════════════════════${NC}"
    echo -e "${GREEN}  ✓ ALL TESTS PASSED!${NC}"
    echo -e "${GREEN}════════════════════════════════════════${NC}"
    echo "Result: ALL TESTS PASSED" >> "$SUMMARY_FILE"
    exit 0
else
    echo -e "${RED}════════════════════════════════════════${NC}"
    echo -e "${RED}  ✗ SOME TESTS FAILED${NC}"
    echo -e "${RED}════════════════════════════════════════${NC}"
    echo "Result: SOME TESTS FAILED" >> "$SUMMARY_FILE"
    exit 1
fi
