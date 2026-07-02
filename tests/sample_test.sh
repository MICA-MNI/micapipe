#!/bin/bash
#
# This script runs micapipe tests.

version=0.2.3
container=singularity

# ---------------------------------------------------------------------------
# CI I/O root — one shared, non-personal location on the runner's scratch
# volume. Everything the CI run reads and writes lives under here (SIF, input
# BIDS, FreeSurfer licence, scratch, outputs), so nothing is tied to a
# personal home dir. Override the root, or any individual path, via the env
# vars below.
# ---------------------------------------------------------------------------
CI_ROOT="${MICAPIPE_CI_ROOT:-/export03/data/action-runner}"

container_img="${MICAPIPE_CI_SIF:-${CI_ROOT}/micapipe_v1_beta.sif}"
bids="${MICAPIPE_CI_BIDS:-${CI_ROOT}/rawdata}"
fs_lic="${MICAPIPE_CI_LICENSE:-${CI_ROOT}/license_fc.txt}"
tmp="${MICAPIPE_CI_TMP:-${CI_ROOT}/tmp}"
mkdir -p "${tmp}" 2>/dev/null || true

# Timestamped output root (never overwrites a previous run).
outdir_base="${MICAPIPE_CI_OUTDIR:-${CI_ROOT}/output}/${container}_${version}"
if ! mkdir -p "${outdir_base}" 2>/dev/null; then
    echo "Warning: ${outdir_base} is not writable. Falling back to /tmp."
    outdir_base="/tmp/micapipe_ci/${container}_${version}"
    mkdir -p "${outdir_base}" || { echo "Failed to create fallback directory"; exit 1; }
fi

timestamp=$(date +%Y%m%d_%H%M%S)
outdir="${outdir_base}_${timestamp}"

function run_test(){
    recon=$1

    if [[ "$recon" == "freesurfer" ]]; then
        out="${outdir}_freesurfer"
        recon="-freesurfer"
    else
        out="${outdir}"
        recon=""
    fi

    mkdir -p "${out}" || { echo "Failed to create ${out}"; exit 1; }
    chmod 777 "${out}" || echo "Warning: could not set permissions for ${out}"

    if [[ "$container" == "docker" ]]; then
      command="docker run -ti --rm -v ${bids}:/bids -v ${out}:/out -v ${tmp}:/tmp -v ${fs_lic}:/opt/licence.txt ${container_img}"
    elif [[ "$container" == "singularity" ]]; then
      command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmp}:/tmp -B ${fs_lic}:/opt/licence.txt ${container_img}"
    fi

    echo "Start running sample test."

    # Run the test for session 01 as an example
    for i in 01; do
      ses=ses-${i}
      sub=sub-mri3T
      ${command} \
      -bids /bids -out /out -fs_licence /opt/licence.txt -threads 15 -sub ${sub} -ses ${ses} \
      -proc_structural -proc_surf -post_structural -proc_dwi -GD -proc_func -MPC -MPC_SWM -SC -SWM -QC_subj -proc_flair \
      -atlas economo,aparc \
      -dwi_rpe /bids/${sub}/${ses}/dwi/${sub}_${ses}_acq-b0_dir-PA_epi.nii.gz -dwi_upsample \
      -func_pe /bids/${sub}/${ses}/fmap/${sub}_${ses}_acq-fmri_dir-AP_epi.nii.gz \
      -func_rpe /bids/${sub}/${ses}/fmap/${sub}_${ses}_acq-fmri_dir-PA_epi.nii.gz \
      -mpc_acq T1map -tracts 10000 \
      -microstructural_img /bids/${sub}/${ses}/anat/${sub}_${ses}_acq-T1_T1map.nii.gz ${recon}
    done
}

run_test "fastsurfer"
run_test "freesurfer"
