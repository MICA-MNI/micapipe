#!/bin/bash
#
# Run All LAMAReg Registration Tests
# Usage: ./run_all_tests.sh <test_data_dir> [output_dir] [threads]
#

# Colors (removed -e flag for compatibility)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
NC='\033[0m'

echo "${MAGENTA}========================================${NC}"
echo "${MAGENTA}  LAMAReg Registration Test Suite${NC}"
echo "${MAGENTA}========================================${NC}"
echo ""

# Parse arguments
if [ $# -lt 1 ]; then
    echo "${RED}Error: Test data directory required${NC}"
    echo "Usage: $0 <test_data_dir> [output_dir] [threads]"
    echo ""
    echo "Arguments:"
    echo "  test_data_dir   Base directory containing test data subdirectories"
    echo "  output_dir      Output directory (default: ./test_results_all)"
    echo "  threads         Number of threads for LAMAReg (default: 16)"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/test/data ./test_results 16"
    exit 1
fi

TEST_DATA_DIR="$1"
OUTPUT_DIR="${2:-$(pwd)/test_results_all}"
THREADS="${3:-16}"

export LAMAREG_THREADS="$THREADS"

# Create output directory
mkdir -p "$OUTPUT_DIR"
MASTER_LOG="$OUTPUT_DIR/master_test_log.txt"
SUMMARY_FILE="$OUTPUT_DIR/test_summary.txt"

echo "LAMAReg Registration Test Suite" > "$SUMMARY_FILE"
echo "Test Date: $(date)" >> "$SUMMARY_FILE"
echo "Test Data: $TEST_DATA_DIR" >> "$SUMMARY_FILE"
echo "Output Dir: $OUTPUT_DIR" >> "$SUMMARY_FILE"
echo "Threads: $THREADS" >> "$SUMMARY_FILE"
echo "========================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

echo "Configuration:"
echo "  Test data: $TEST_DATA_DIR"
echo "  Output: $OUTPUT_DIR"
echo "  Threads: $THREADS"
echo ""

# Test scripts and their corresponding data subdirectories
TESTS=(
    "test_dwi_registration.sh:dwi"
    "test_func_registration.sh:func"
    "test_mpc_registration.sh:mpc"
    "test_mpc_swm_registration.sh:mpc"
)

# Test names
TEST_NAMES=(
    "DWI Registration"
    "Functional MRI Registration"
    "MPC Registration"
    "MPC-SWM Registration"
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
    # Parse test:subdir format
    IFS=':' read -r test_script test_subdir <<< "${TESTS[$i]}"
    test_name="${TEST_NAMES[$i]}"
    test_data_path="$TEST_DATA_DIR/$test_subdir"
    test_output_dir="$OUTPUT_DIR/$(basename "$test_script" .sh)"
    
    # Check if test data directory exists
    if [ ! -d "$test_data_path" ]; then
        echo ""
        echo "${YELLOW}========================================${NC}"
        echo "${YELLOW}Skipping: $test_name${NC}"
        echo "${YELLOW}========================================${NC}"
        echo "${YELLOW}Data directory not found: $test_data_path${NC}"
        echo ""
        echo "⊘ SKIPPED: $test_name (no data)" >> "$SUMMARY_FILE"
        SKIPPED_TESTS=$((SKIPPED_TESTS + 1))
        continue
    fi
    
    echo ""
    echo "${BLUE}========================================${NC}"
    echo "${BLUE}Running: $test_name${NC}"
    echo "${BLUE}Data: $test_data_path${NC}"
    echo "${BLUE}Output: $test_output_dir${NC}"
    echo "${BLUE}========================================${NC}"
    echo ""
    
    # Make script executable
    chmod +x "$SCRIPT_DIR/$test_script"
    
    # Run test and capture exit code
    "$SCRIPT_DIR/$test_script" "$test_data_path" "$test_output_dir" 2>&1 | tee -a "$MASTER_LOG"
    test_exit_code=$?
    
    # Check result
    if [ $test_exit_code -eq 0 ]; then
        echo "${GREEN}✓ PASSED: $test_name${NC}"
        echo "✓ PASSED: $test_name" >> "$SUMMARY_FILE"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        echo "${RED}✗ FAILED: $test_name (exit code: $test_exit_code)${NC}"
        echo "✗ FAILED: $test_name (exit code: $test_exit_code)" >> "$SUMMARY_FILE"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
done

# Final summary
echo ""
echo "${MAGENTA}========================================${NC}"
echo "${MAGENTA}  Test Suite Summary${NC}"
echo "${MAGENTA}========================================${NC}"
echo "Total tests: ${TOTAL_TESTS}"
echo "${GREEN}Passed: ${PASSED_TESTS}${NC}"
echo "${RED}Failed: ${FAILED_TESTS}${NC}"
echo "${YELLOW}Skipped: ${SKIPPED_TESTS}${NC}"
echo ""
echo "Detailed results: ${OUTPUT_DIR}"
echo "Master log: ${MASTER_LOG}"
echo "Summary: ${SUMMARY_FILE}"
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

if [ $FAILED_TESTS -eq 0 ] && [ $SKIPPED_TESTS -lt $TOTAL_TESTS ]; then
    echo "${GREEN}════════════════════════════════════════${NC}"
    echo "${GREEN}  ✓ ALL TESTS PASSED!${NC}"
    echo "${GREEN}════════════════════════════════════════${NC}"
    echo "Result: ALL TESTS PASSED" >> "$SUMMARY_FILE"
    exit 0
else
    echo "${RED}════════════════════════════════════════${NC}"
    echo "${RED}  ✗ SOME TESTS FAILED OR SKIPPED${NC}"
    echo "${RED}════════════════════════════════════════${NC}"
    echo "Result: SOME TESTS FAILED OR SKIPPED" >> "$SUMMARY_FILE"
    exit 1
fi
