#!/bin/bash
#
# This script runs micapipe tests.

version=0.2.3
container=singularity
container_img=/data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v0.2.3.sif

bids=/data/mica3/BIDS_CI/rawdata
fs_lic=/data_/mica3/BIDS_CI/license_fc.txt
tmp=/tmp

# Define a primary output directory
outdir_base="/data/mica1/03_projects/enning/BIDS_CI/${container}_${version}"

# Test if outdir_base is writable; if not, fall back to /tmp
if ! mkdir -p "${outdir_base}" 2>/dev/null; then
    echo "Warning: ${outdir_base} is not writable. Falling back to /tmp."
    outdir_base="/export02/local/singularity_tmp/${container}_${version}"
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
