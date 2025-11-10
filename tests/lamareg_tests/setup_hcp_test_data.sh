#!/bin/bash
#
# Setup LAMAReg Test Data from HCP Subject
# 
# Usage: ./setup_hcp_test_data.sh
#

set -e

# Source directories
DERIV_DIR="/data_/mica3/BIDS_HCP-57sub/derivatives/micapipe_v0.2.0"
SUB="sub-211215"
SOURCE="${DERIV_DIR}/${SUB}"

# Target directory for test data
TEST_DIR="/data_/mica3/BIDS_CI/lamareg_test_data"

echo "========================================"
echo "Setting up LAMAReg test data"
echo "========================================"
echo "Source: ${SOURCE}"
echo "Target: ${TEST_DIR}"
echo ""

# Create test data directories
mkdir -p ${TEST_DIR}/{dwi,func,flair,mpc}

echo "Step 1: DWI Registration Test Data"
echo "-----------------------------------"
# For DWI test - need T1w in DWI space and a 3D DWI reference
# IMPORTANT: LAMAReg/SynthSeg requires 3D images, not 4D FOD!
if [ -f "${SOURCE}/dwi/${SUB}_space-dwi_desc-T1w_nativepro_SyN.nii.gz" ]; then
    echo "  ✓ Found T1w in DWI space"
    cp "${SOURCE}/dwi/${SUB}_space-dwi_desc-T1w_nativepro_SyN.nii.gz" \
       "${TEST_DIR}/dwi/T1w_in_dwi_brain.nii.gz"
else
    echo "  ✗ Missing: T1w in DWI space"
fi

# Try to find a 3D DWI reference image (in priority order)
# 1. b0 image (best for registration)
# 2. FA map (good white matter contrast)
# 3. Extract first volume from FOD (fallback)

if [ -f "${SOURCE}/dwi/${SUB}_space-dwi_desc-b0_dwi.nii.gz" ]; then
    echo "  ✓ Found b0 image (3D - ideal for LAMAReg)"
    cp "${SOURCE}/dwi/${SUB}_space-dwi_desc-b0_dwi.nii.gz" \
       "${TEST_DIR}/dwi/dwi_b0.nii.gz"
elif [ -f "${SOURCE}/dwi/${SUB}_space-dwi_model-DTI_map-FA.nii.gz" ]; then
    echo "  ✓ Found FA map (3D - good for LAMAReg)"
    cp "${SOURCE}/dwi/${SUB}_space-dwi_model-DTI_map-FA.nii.gz" \
       "${TEST_DIR}/dwi/dwi_FA.nii.gz"
elif [ -f "${SOURCE}/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" ]; then
    echo "  ⚠ Found 4D FOD - extracting first volume as 3D reference"
    echo "    (Note: Using full b0 or FA would be better)"
    # Extract first volume from 4D FOD to create 3D reference
    if command -v fslroi &> /dev/null; then
        fslroi "${SOURCE}/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" \
               "${TEST_DIR}/dwi/dwi_fod_vol0.nii.gz" 0 1
        echo "  ✓ Created 3D reference from FOD (first volume)"
    else
        echo "  ✗ fslroi not available - cannot extract 3D from 4D FOD"
        # Copy the 4D FOD anyway but warn it will fail
        cp "${SOURCE}/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" \
           "${TEST_DIR}/dwi/dwi_fod.nii.gz"
        echo "  ✗ WARNING: Copied 4D FOD - this WILL FAIL with SynthSeg!"
    fi
else
    echo "  ✗ Missing: No DWI reference image found (b0, FA, or FOD)"
fi

echo ""
echo "Step 2: Functional MRI Test Data"
echo "-----------------------------------"
# For functional test - need func brain and T1 nativepro
# HCP data has nested directory structure
FUNC_BRAIN=$(find ${SOURCE}/func -name "*_space-func_desc-*_brain.nii.gz" 2>/dev/null | head -1)
if [ -n "$FUNC_BRAIN" ] && [ -f "$FUNC_BRAIN" ]; then
    echo "  ✓ Found functional MRI brain"
    cp "$FUNC_BRAIN" "${TEST_DIR}/func/func_brain.nii.gz"
else
    echo "  ✗ Missing: Functional MRI brain"
fi

if [ -f "${SOURCE}/anat/${SUB}_space-nativepro_T1w.nii.gz" ]; then
    echo "  ✓ Found T1 nativepro"
    cp "${SOURCE}/anat/${SUB}_space-nativepro_T1w.nii.gz" \
       "${TEST_DIR}/func/T1_nativepro.nii.gz"
else
    echo "  ✗ Missing: T1 nativepro"
fi

echo ""
echo "Step 3: FLAIR Test Data (if available)"
echo "-----------------------------------"
# For FLAIR test - may not always be present
if [ -f "${SOURCE}/anat/${SUB}_space-flair_FLAIR.nii.gz" ]; then
    echo "  ✓ Found FLAIR preprocessed"
    cp "${SOURCE}/anat/${SUB}_space-flair_FLAIR.nii.gz" \
       "${TEST_DIR}/flair/FLAIR_preproc.nii.gz"
    cp "${SOURCE}/anat/${SUB}_space-nativepro_T1w.nii.gz" \
       "${TEST_DIR}/flair/T1_nativepro.nii.gz"
else
    echo "  ⚠ FLAIR not available (this is optional)"
fi

echo ""
echo "Step 4: MPC Test Data"
echo "-----------------------------------"
# For MPC test - need microstructural map and T1 images
# Use FA map as microstructural map
if [ -f "${SOURCE}/dwi/${SUB}_space-dwi_model-DTI_map-FA.nii.gz" ]; then
    echo "  ✓ Found FA map (microstructural)"
    cp "${SOURCE}/dwi/${SUB}_space-dwi_model-DTI_map-FA.nii.gz" \
       "${TEST_DIR}/mpc/microstructural_map.nii.gz"
else
    echo "  ✗ Missing: FA map"
fi

# T1 nativepro for MPC-SWM
if [ -f "${SOURCE}/anat/${SUB}_space-nativepro_T1w.nii.gz" ]; then
    echo "  ✓ Found T1 nativepro (for MPC-SWM)"
    cp "${SOURCE}/anat/${SUB}_space-nativepro_T1w.nii.gz" \
       "${TEST_DIR}/mpc/T1_nativepro.nii.gz"
fi

# FreeSurfer T1 for MPC
if [ -f "${SOURCE}/anat/${SUB}_space-fsnative_T1w.nii.gz" ]; then
    echo "  ✓ Found T1 in FreeSurfer native space"
    cp "${SOURCE}/anat/${SUB}_space-fsnative_T1w.nii.gz" \
       "${TEST_DIR}/mpc/T1_fsnative.nii.gz"
else
    echo "  ✗ Missing: T1 in FreeSurfer native space"
fi

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo ""

# Count files created
TOTAL_FILES=$(find ${TEST_DIR} -name "*.nii.gz" | wc -l)
echo "Total files created: $TOTAL_FILES"
echo ""

echo "Files by modality:"
for dir in dwi func flair mpc; do
    count=$(ls ${TEST_DIR}/${dir}/*.nii.gz 2>/dev/null | wc -l)
    echo "  ${dir}: ${count} files"
done

echo ""
echo "Test data ready at: ${TEST_DIR}"
echo ""
echo "Next steps:"
echo "  1. Navigate to micapipe test directory:"
echo "     cd /path/to/micapipe/tests/lamareg_tests"
echo ""
echo "  2. Run individual tests:"
echo "     ./test_dwi_registration.sh ${TEST_DIR}/dwi"
echo "     ./test_func_registration.sh ${TEST_DIR}/func"
echo "     ./test_mpc_registration.sh ${TEST_DIR}/mpc"
echo "     ./test_mpc_swm_registration.sh ${TEST_DIR}/mpc"
echo ""
echo "  3. Or run all tests:"
echo "     ./run_all_tests.sh ${TEST_DIR}"
echo ""
echo "========================================"
