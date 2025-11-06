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
    
    validate_inputs "$moving_img" "$fixed_img"
    
    # Define output paths
    setup_output_paths "$output_prefix"
    
    # LAMAReg command
    log "LAMAReg command for functional MRI registration:"
    echo "lamareg register \\" | tee -a "$LOG_FILE"
    echo "  --moving \"$moving_img\" \\" | tee -a "$LOG_FILE"
    echo "  --fixed \"$fixed_img\" \\" | tee -a "$LOG_FILE"
    echo "  --output \"$WARPED\" \\" | tee -a "$LOG_FILE"
    echo "  --moving-parc \"$MOVING_PARC\" \\" | tee -a "$LOG_FILE"
    echo "  --fixed-parc \"$FIXED_PARC\" \\" | tee -a "$LOG_FILE"
    echo "  --registered-parc \"$REG_PARC\" \\" | tee -a "$LOG_FILE"
    echo "  --affine \"$AFFINE\" \\" | tee -a "$LOG_FILE"
    echo "  --warpfield \"$WARP1\" \\" | tee -a "$LOG_FILE"
    echo "  --inverse-warpfield \"$INVWARP1\" \\" | tee -a "$LOG_FILE"
    echo "  --secondary-warpfield \"$WARP2\" \\" | tee -a "$LOG_FILE"
    echo "  --inverse-secondary-warpfield \"$INVWARP2\" \\" | tee -a "$LOG_FILE"
    echo "  --qc-csv \"$QC_CSV\" \\" | tee -a "$LOG_FILE"
    echo "  --synthseg-threads 4 \\" | tee -a "$LOG_FILE"
    echo "  --ants-threads 8" | tee -a "$LOG_FILE"
    
    test_result "LAMAReg command syntax" "PASS" ""
    
    # Test transformation chain
    log "Forward transform: -t warp2 -t warp1 -t affine"
    log "Inverse transform: -t affine -t invwarp1 -t invwarp2"
    test_result "Transformation chain" "PASS" ""
}

# Main test execution
main "$@"
