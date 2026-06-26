#!/bin/bash
#
# This script runs micapipe on the same subject on:
# SINGLE session 3T
# multi-session 3T
# multi-session 7T

# input arguments
version=$1
container=$2
container_img=$3

# Local variables
bids=/data/mica3/BIDS_CI/rawdata
fs_lic=/data_/mica3/BIDS_CI/license_fc.txt
tmp=/tmp
outdir=/data/mica3/BIDS_CI/${container}_${version}
threads=15
tracts=100000

# ------------------------------------------------------------------------------
function run_test(){
    recon=$1

    # Conditional statement for freesurfer/fastsurfer run
    if [[ "$recon" == "freesurfer" ]]; then
        out="${outdir}_freesurfer"
        recon="-freesurfer"
    else
        recon=""
        out="${outdir}"
    fi

    # Create output directory
    mkdir ${out}
    chmod 777 ${out}

    # Create command string
    if [[ "$container" == "docker" ]]; then
      command="docker run -ti --rm -v ${bids}:/bids -v ${out}:/out -v ${tmp}:/tmp -v ${fs_lic}:/opt/licence.txt ${container_img}"
    elif [[ "$container" == "singularity" ]]; then
      command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmp}:/tmp -B ${fs_lic}:/opt/licence.txt ${container_img}"
    fi

    # ------------------------------------------------------------------------------
    # SINGLE session one shot workflow
    sub=sub-mri

    ${command} \
    -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} \
    -proc_structural -proc_surf -post_structural -proc_dwi -GD -proc_func -SC -SWM -QC_subj -proc_flair \
    -atlas economo,aparc -flairScanStr T2w -mainScanStr task-music_bold -NSR -noFIX \
    -dwi_main /bids/${sub}/dwi/${sub}_acq-b1000_dir-PA_dwi.nii.gz \
    -tracts ${tracts} ${recon}

    # Run fMRI additional acquisition
    ${command} \
    -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} \
    -proc_func -mainScanStr fingertap -NSR -noFIX ${recon}

    # Session 01 and 02
    for i in 01 02; do
    ses=ses-${i}

    # 3T multi session one shot workflow -fastsurfer
    sub=sub-mri3T
    ${command} \
    -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
    -proc_structural -proc_surf -post_structural -proc_dwi -GD -proc_func -MPC -MPC_SWM -SC -SWM -QC_subj -proc_flair \
    -atlas economo,aparc \
    -dwi_rpe /bids/${sub}/${ses}/dwi/${sub}_${ses}_acq-b0_dir-PA_epi.nii.gz -dwi_upsample \
    -func_pe /bids/${sub}/${ses}/fmap/${sub}_${ses}_acq-fmri_dir-AP_epi.nii.gz \
    -func_rpe /bids/${sub}/${ses}/fmap/${sub}_${ses}_acq-fmri_dir-PA_epi.nii.gz \
    -mpc_acq T1map -tracts ${tracts} \
    -microstructural_img /bids/${sub}/${ses}/anat/${sub}_${ses}_acq-T1_T1map.nii.gz ${recon}

    # 7T multi session one shot workflow
    sub=sub-mri7T
    ${command} \
    -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
    -proc_structural -proc_surf -post_structural -proc_dwi -GD -proc_func -MPC -MPC_SWM -SC -SWM -QC_subj \
    -uni -T1wStr acq-uni_T1map -atlas economo,aparc \
    -dwi_rpe /bids/${sub}/${ses}/dwi/${sub}_${ses}_acq-b0_dir-PA_epi.nii.gz -dwi_upsample \
    -mainScanStr task-rest_echo-1_bold,task-rest_echo-2_bold,task-rest_echo-3_bold \
    -func_pe /bids/${sub}/${ses}/fmap/${sub}_${ses}_acq-fmri_dir-AP_epi.nii.gz \
    -func_rpe /bids/${sub}/${ses}/fmap/${sub}_${ses}_acq-fmri_dir-PA_epi.nii.gz \
    -mpc_acq T1map -tracts ${tracts} \
    -microstructural_img /bids/${sub}/${ses}/anat/${sub}_${ses}_acq-T1_T1map.nii.gz ${recon}

    done

    # ------------------------------------------------------------------------------
    # Run Group QC
    ${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -QC
}

# Run test with fastsurfer
run_test "fastsurfer"

# Run test with freesurfer
run_test "freesurfer"