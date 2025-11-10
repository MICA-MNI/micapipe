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
# For DWI test - need T1w in DWI space and FOD
if [ -f "${SOURCE}/dwi/${SUB}_space-dwi_desc-T1w_nativepro-brain.nii.gz" ]; then
    echo "  ✓ Found T1w in DWI space (brain)"
    cp "${SOURCE}/dwi/${SUB}_space-dwi_desc-T1w_nativepro-brain.nii.gz" \
       "${TEST_DIR}/dwi/T1w_in_dwi_brain.nii.gz"
else
    echo "  ✗ Missing: T1w in DWI space"
fi

if [ -f "${SOURCE}/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" ]; then
    echo "  ✓ Found DWI FOD (white matter normalized)"
    cp "${SOURCE}/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" \
       "${TEST_DIR}/dwi/dwi_fod.nii.gz"
else
    echo "  ✗ Missing: DWI FOD"
fi

echo ""
echo "Step 2: Functional MRI Test Data"
echo "-----------------------------------"
# For functional test - need func brain and T1 nativepro
FUNC_FILE=$(ls ${SOURCE}/func/${SUB}_task-*_space-func_desc-brain_bold.nii.gz 2>/dev/null | head -1)
if [ -f "$FUNC_FILE" ]; then
    echo "  ✓ Found functional MRI brain"
    cp "$FUNC_FILE" "${TEST_DIR}/func/func_brain.nii.gz"
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
FS_ORIG="${SOURCE}/anat/surfaces/freesurfer/${SUB}/mri/orig.mgz"
if [ -f "$FS_ORIG" ]; then
    echo "  ✓ Found FreeSurfer orig.mgz, converting..."
    mri_convert "$FS_ORIG" "${TEST_DIR}/mpc/T1_fsnative.nii.gz"
else
    echo "  ✗ Missing: FreeSurfer orig.mgz"
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
