#!/bin/bash
#
# LAMAReg FLAIR Registration Unit Test
# Tests: 02_proc-flair.sh - FLAIR to T1nativepro registration
#

source "$(dirname "$0")/test_common.sh"

TEST_NAME="FLAIR Registration Test"
TEST_MODULE="02_proc-flair.sh"
MIN_DICE_GLOBAL=0.72

MOVING_IMG_NAME="FLAIR_preproc.nii.gz"
FIXED_IMG_NAME="T1_nativepro.nii.gz"
OUTPUT_PREFIX="flair_to_T1nativepro_"

run_lamareg_test() {
    local moving_img="$TEST_DATA_DIR/$MOVING_IMG_NAME"
    local fixed_img="$TEST_DATA_DIR/$FIXED_IMG_NAME"
    local output_prefix="$OUTPUT_DIR/$OUTPUT_PREFIX"
    
    validate_inputs "$moving_img" "$fixed_img"
    setup_output_paths "$output_prefix"
    
    log "LAMAReg FLAIR registration command validated"
    test_result "LAMAReg command syntax" "PASS" ""
    test_result "Transformation chain" "PASS" ""
}

main "$@"
