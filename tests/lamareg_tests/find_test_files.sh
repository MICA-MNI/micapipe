#!/bin/bash
#
# Diagnostic script to find available test files
#

SOURCE="/data_/mica3/BIDS_HCP-57sub/derivatives/micapipe_v0.2.0/sub-211215"

echo "========================================="
echo "Finding available test files"
echo "========================================="
echo ""

echo "DWI Directory Files:"
echo "-------------------"
ls -1 ${SOURCE}/dwi/*.nii.gz 2>/dev/null | head -20

echo ""
echo "Functional Directory Files:"
echo "-------------------"
find ${SOURCE}/func -name "*.nii.gz" 2>/dev/null | head -20

echo ""
echo "Anat Directory Files:"
echo "-------------------"
ls -1 ${SOURCE}/anat/*.nii.gz 2>/dev/null | head -20

echo ""
echo "FreeSurfer Directory:"
echo "-------------------"
ls -1 ${SOURCE}/anat/surfaces/freesurfer/sub-211215/mri/*.mgz 2>/dev/null | head -10

