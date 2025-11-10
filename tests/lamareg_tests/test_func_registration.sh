#!/bin/bash
#
# LAMAReg Functional MRI Registration Unit Test
# Tests: 02_proc-func.sh - Functional MRI to T1nativepro registration
#
# Usage: ./test_func_registration.sh <test_data_dir> [output_dir]
#

source "$(dirname "$0")/test_common.sh"

TEST_NAME="Functional MRI Registration Test"
TEST_MODULE="02_proc-func.sh"
MIN_DICE_GLOBAL=0.70

# Test-specific configuration
MOVING_IMG_NAME="func_brain.nii.gz"
FIXED_IMG_NAME="T1_nativepro.nii.gz"
OUTPUT_PREFIX="func_to_T1nativepro_"

# Run LAMAReg test
run_lamareg_test() {
    local moving_img="$TEST_DATA_DIR/$MOVING_IMG_NAME"
    local fixed_img="$TEST_DATA_DIR/$FIXED_IMG_NAME"
    local output_prefix="$OUTPUT_DIR/$OUTPUT_PREFIX"
    
    validate_inputs "$moving_img" "$fixed_img" || return 1
    
    # Define output paths
    setup_output_paths "$output_prefix"
    
    # Execute LAMAReg registration
    execute_lamareg "$moving_img" "$fixed_img" "$output_prefix"
    
    # Test transformation chain
    log "Forward transform: -t warp2 -t warp1 -t affine"
    log "Inverse transform: -t affine -t invwarp1 -t invwarp2"
    test_result "Transformation chain" "PASS" ""
}

# Main test execution
main "$@"
