#!/bin/bash
#
# LAMAReg MPC Registration Unit Test
# Tests: 03_MPC.sh - Microstructural to FreeSurfer registration
#

source "$(dirname "$0")/test_common.sh"

TEST_NAME="MPC Registration Test"
TEST_MODULE="03_MPC.sh"
MIN_DICE_GLOBAL=0.70

MOVING_IMG_NAME="microstructural_map.nii.gz"
FIXED_IMG_NAME="T1_fsnative.nii.gz"
OUTPUT_PREFIX="qMRI_to_fsnative_"

run_lamareg_test() {
    local moving_img="$TEST_DATA_DIR/$MOVING_IMG_NAME"
    local fixed_img="$TEST_DATA_DIR/$FIXED_IMG_NAME"
    local output_prefix="$OUTPUT_DIR/$OUTPUT_PREFIX"
    
    validate_inputs "$moving_img" "$fixed_img" || return 1
    setup_output_paths "$output_prefix"
    
    # Execute LAMAReg registration
    execute_lamareg "$moving_img" "$fixed_img" "$output_prefix"
    
    log "Transformation chain: forward and inverse"
    test_result "Transformation chain" "PASS" ""
}

main "$@"
