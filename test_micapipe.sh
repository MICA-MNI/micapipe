#!/bin/bash
#
# Test micapipe Docker/Singularity image with sample data
#
# This script runs micapipe processing tests on the server and outputs
# results to a dedicated test directory.
#
# Usage:
#   ./test_micapipe.sh [docker|singularity] [image_path]
#
# Examples:
#   ./test_micapipe.sh singularity                    # Use default SIF path
#   ./test_micapipe.sh singularity /path/to/image.sif # Use specific SIF
#   ./test_micapipe.sh docker micapipe:latest         # Test Docker image
#

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================
VERSION="v0.2.3"
CONTAINER_TYPE="${1:-singularity}"  # docker or singularity

# Default image paths
DEFAULT_SINGULARITY_IMG="/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif"
DEFAULT_DOCKER_IMG="ghcr.io/mica-mni/micapipe:latest"

# Get image path from argument or use default
if [ "$CONTAINER_TYPE" == "singularity" ]; then
    CONTAINER_IMG="${2:-$DEFAULT_SINGULARITY_IMG}"
elif [ "$CONTAINER_TYPE" == "docker" ]; then
    CONTAINER_IMG="${2:-$DEFAULT_DOCKER_IMG}"
else
    echo "❌ ERROR: Container type must be 'docker' or 'singularity'"
    echo "Usage: $0 [docker|singularity] [image_path]"
    exit 1
fi

# Test data paths
BIDS_DIR="/data/mica3/BIDS_CI/rawdata"
FS_LICENSE="/data_/mica3/BIDS_CI/license_fc.txt"
TMP_DIR="/tmp"

# Output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
TEST_BASE="/host/cassio/export03/data/enning/test_results"
TEST_OUTPUT="${TEST_BASE}/${CONTAINER_TYPE}_${VERSION}_${TIMESTAMP}"

# Test subjects/sessions
TEST_SUBJECT="sub-mri3T"
TEST_SESSION="ses-01"

# ============================================================================
# Header
# ============================================================================
echo "============================================="
echo "🧪 MICAPIPE TEST SUITE"
echo "============================================="
echo ""
echo "📦 Container Type: ${CONTAINER_TYPE}"
echo "🖼️  Container Image: ${CONTAINER_IMG}"
echo "📊 Test Subject:    ${TEST_SUBJECT}"
echo "📊 Test Session:    ${TEST_SESSION}"
echo "📍 Output Base:     ${TEST_OUTPUT}"
echo ""

# ============================================================================
# Pre-flight Checks
# ============================================================================
echo "🔍 Running pre-flight checks..."
echo ""

# Check if container runtime is available
if [ "$CONTAINER_TYPE" == "singularity" ]; then
    if ! command -v singularity &> /dev/null; then
        echo "❌ ERROR: singularity command not found"
        exit 1
    fi
    echo "✅ Singularity found: $(singularity --version)"
    
    # Check if SIF exists
    if [ ! -f "$CONTAINER_IMG" ]; then
        echo "❌ ERROR: Singularity image not found: ${CONTAINER_IMG}"
        echo ""
        echo "Available SIF files:"
        find /host/cassio/export03/data/enning/singularity -name "*.sif" -exec ls -lh {} \; 2>/dev/null || echo "  (none found)"
        exit 1
    fi
    echo "✅ SIF image found: ${CONTAINER_IMG}"
    
elif [ "$CONTAINER_TYPE" == "docker" ]; then
    if ! command -v docker &> /dev/null; then
        echo "❌ ERROR: docker command not found"
        exit 1
    fi
    echo "✅ Docker found: $(docker --version)"
    
    # Check if Docker image exists
    if ! docker image inspect "$CONTAINER_IMG" &> /dev/null; then
        echo "❌ ERROR: Docker image not found: ${CONTAINER_IMG}"
        echo ""
        echo "Available micapipe images:"
        docker images | grep micapipe || echo "  (none found)"
        exit 1
    fi
    echo "✅ Docker image found: ${CONTAINER_IMG}"
fi

# Check if test data exists
if [ ! -d "$BIDS_DIR" ]; then
    echo "❌ ERROR: BIDS test data not found: ${BIDS_DIR}"
    exit 1
fi
echo "✅ BIDS data found: ${BIDS_DIR}"

if [ ! -f "$FS_LICENSE" ]; then
    echo "❌ ERROR: FreeSurfer license not found: ${FS_LICENSE}"
    exit 1
fi
echo "✅ FreeSurfer license found: ${FS_LICENSE}"

# Check if test subject/session exists
SUBJECT_DIR="${BIDS_DIR}/${TEST_SUBJECT}/${TEST_SESSION}"
if [ ! -d "$SUBJECT_DIR" ]; then
    echo "❌ ERROR: Test subject/session not found: ${SUBJECT_DIR}"
    exit 1
fi
echo "✅ Test subject/session found: ${SUBJECT_DIR}"

# Create output directory
echo ""
echo "📁 Creating test output directory..."
mkdir -p "${TEST_OUTPUT}" || { echo "❌ Failed to create output directory"; exit 1; }
chmod 777 "${TEST_OUTPUT}" 2>/dev/null || true
echo "✅ Output directory created: ${TEST_OUTPUT}"

# ============================================================================
# Test Functions
# ============================================================================

# Function to run micapipe test
run_micapipe_test() {
    local recon_type=$1
    local recon_flag=""
    local test_name=""
    
    if [[ "$recon_type" == "freesurfer" ]]; then
        test_name="freesurfer"
        recon_flag="-freesurfer"
        outdir="${TEST_OUTPUT}_freesurfer"
    else
        test_name="fastsurfer"
        recon_flag=""  # Default is FastSurfer
        outdir="${TEST_OUTPUT}_fastsurfer"
    fi
    
    echo ""
    echo "============================================="
    echo "🚀 Running Test: ${test_name}"
    echo "============================================="
    echo "📍 Output: ${outdir}"
    echo ""
    
    # Create output directory
    mkdir -p "${outdir}" || { echo "❌ Failed to create ${outdir}"; exit 1; }
    chmod 777 "${outdir}" 2>/dev/null || true
    
    # Build container command
    if [[ "$CONTAINER_TYPE" == "docker" ]]; then
        CONTAINER_CMD="docker run -ti --rm \
            -v ${BIDS_DIR}:/bids \
            -v ${outdir}:/out \
            -v ${TMP_DIR}:/tmp \
            -v ${FS_LICENSE}:/opt/licence.txt \
            ${CONTAINER_IMG}"
    elif [[ "$CONTAINER_TYPE" == "singularity" ]]; then
        CONTAINER_CMD="singularity run --writable-tmpfs --containall \
            -B ${BIDS_DIR}:/bids \
            -B ${outdir}:/out \
            -B ${TMP_DIR}:/tmp \
            -B ${FS_LICENSE}:/opt/licence.txt \
            ${CONTAINER_IMG}"
    fi
    
    # Record start time
    START_TIME=$(date +%s)
    
    # Run micapipe
    echo "⏱️  Test started at: $(date)"
    echo ""
    
    ${CONTAINER_CMD} \
        -bids /bids \
        -out /out \
        -fs_licence /opt/licence.txt \
        -threads 15 \
        -sub ${TEST_SUBJECT} \
        -ses ${TEST_SESSION} \
        -proc_structural \
        -proc_surf \
        -post_structural \
        -proc_dwi \
        -GD \
        -proc_func \
        -MPC \
        -MPC_SWM \
        -SC \
        -SWM \
        -QC_subj \
        -proc_flair \
        -atlas economo,aparc \
        -dwi_rpe /bids/${TEST_SUBJECT}/${TEST_SESSION}/dwi/${TEST_SUBJECT}_${TEST_SESSION}_acq-b0_dir-PA_epi.nii.gz \
        -dwi_upsample \
        -func_pe /bids/${TEST_SUBJECT}/${TEST_SESSION}/fmap/${TEST_SUBJECT}_${TEST_SESSION}_acq-fmri_dir-AP_epi.nii.gz \
        -func_rpe /bids/${TEST_SUBJECT}/${TEST_SESSION}/fmap/${TEST_SUBJECT}_${TEST_SESSION}_acq-fmri_dir-PA_epi.nii.gz \
        -mpc_acq T1map \
        -regSynth \
        -tracts 10000 \
        -microstructural_img /bids/${TEST_SUBJECT}/${TEST_SESSION}/anat/${TEST_SUBJECT}_${TEST_SESSION}_acq-T1_T1map.nii.gz \
        -microstructural_reg /bids/${TEST_SUBJECT}/${TEST_SESSION}/anat/${TEST_SUBJECT}_${TEST_SESSION}_acq-inv1_T1map.nii.gz \
        ${recon_flag}
    
    TEST_EXIT_CODE=$?
    
    # Record end time
    END_TIME=$(date +%s)
    ELAPSED_TIME=$((END_TIME - START_TIME))
    ELAPSED_MIN=$((ELAPSED_TIME / 60))
    ELAPSED_SEC=$((ELAPSED_TIME % 60))
    
    echo ""
    echo "============================================="
    if [ ${TEST_EXIT_CODE} -eq 0 ]; then
        echo "✅ ${test_name} TEST PASSED"
        echo "============================================="
        echo "⏱️  Test time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
        echo "📍 Results: ${outdir}"
        echo ""
        
        # Check for key output files
        echo "🔍 Checking output files..."
        
        if [ -d "${outdir}/micapipe_v${VERSION}" ]; then
            echo "✅ micapipe output directory exists"
            
            # Count key output types
            NUM_LOGS=$(find "${outdir}" -name "*.log" 2>/dev/null | wc -l)
            NUM_NIFTIS=$(find "${outdir}" -name "*.nii.gz" 2>/dev/null | wc -l)
            NUM_SURFACES=$(find "${outdir}" -name "*.surf.gii" 2>/dev/null | wc -l)
            
            echo "   - Log files: ${NUM_LOGS}"
            echo "   - NIfTI files: ${NUM_NIFTIS}"
            echo "   - Surface files: ${NUM_SURFACES}"
        else
            echo "⚠️  Warning: Expected output directory not found"
        fi
        
        return 0
    else
        echo "❌ ${test_name} TEST FAILED"
        echo "============================================="
        echo "⏱️  Test time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
        echo "❌ Exit code: ${TEST_EXIT_CODE}"
        echo "📍 Partial results may be in: ${outdir}"
        echo ""
        echo "Check logs for errors:"
        echo "   find ${outdir} -name '*.log' -exec tail -50 {} \\;"
        return 1
    fi
}

# ============================================================================
# Run Tests
# ============================================================================
echo ""
echo "============================================="
echo "🧪 STARTING TEST SUITE"
echo "============================================="
echo ""

# Test 1: FastSurfer (default)
run_micapipe_test "fastsurfer"
FASTSURFER_RESULT=$?

echo ""
echo "---"
echo ""

# Test 2: FreeSurfer
run_micapipe_test "freesurfer"
FREESURFER_RESULT=$?

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "============================================="
echo "📊 TEST SUITE SUMMARY"
echo "============================================="
echo ""

if [ ${FASTSURFER_RESULT} -eq 0 ]; then
    echo "✅ FastSurfer Test: PASSED"
else
    echo "❌ FastSurfer Test: FAILED"
fi

if [ ${FREESURFER_RESULT} -eq 0 ]; then
    echo "✅ FreeSurfer Test: PASSED"
else
    echo "❌ FreeSurfer Test: FAILED"
fi

echo ""
echo "📍 Test Results Location:"
echo "   ${TEST_OUTPUT}_fastsurfer"
echo "   ${TEST_OUTPUT}_freesurfer"
echo ""

# Overall result
if [ ${FASTSURFER_RESULT} -eq 0 ] && [ ${FREESURFER_RESULT} -eq 0 ]; then
    echo "🎉 ALL TESTS PASSED!"
    exit 0
else
    echo "❌ SOME TESTS FAILED"
    exit 1
fi
