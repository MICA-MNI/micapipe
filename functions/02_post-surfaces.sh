#!/bin/bash
#
# DWI structural processing with bash:
#
# Preprocessing workflow for diffusion MRI.
#
# This workflow makes use of MRtrix3
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out Directory
#   $4 : Temporal directory (default /tmp)
#
# ONLY for scripting and debugging:
#TEST=ON
# source utilities
source $MICAPIPE/functions/utilities.sh

BIDS=$1
id=$2
out=$3
tmp=$4

# Assigns variables names
bids_variables $BIDS $id $out
# print the names on the terminal
bids_print.variables-post

# Check inputs: Nativepro T1
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi

# Check inputs: freesurfer space T1
if [ ! -f ${T1freesurfr} ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi

# Sets wb_command to only use one thread
if [[ -z $OMP_NUM_THREADS ]]; then OMP_NUM_THREADS=4; Info "wb_command will use $OMP_NUM_THREADS threads"; fi

#------------------------------------------------------------------------------#
Title "Running MICA POST-structural processing"

#	Timer
aloita=$(date +%s)

# Check tmp dir: temporary directory
random_str=$RANDOM
if [ -z ${tmp} ]; then tmp=/tmp/${random_str}_post_structural_${id}; fi
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Freesurface SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Compute warp of native structural to Freesurfer and apply to 5TT and first
Info "Native surfaces to conte69-64k vertices (both hemispheres)"
if [[ ! -f  ${dir_conte69}/${id}_rh_midthickness_32k_fs_LR_fsspace_cras_corrected.surf.gii ]] ; then
    for hemisphere in l r; do
      HEMI=`echo $hemisphere | tr [:lower:] [:upper:]`
        # Build the conte69-32k sphere and midthickness surface
        Do_cmd wb_shortcuts -freesurfer-resample-prep \
            ${dir_freesurfer}/surf/${hemisphere}h.white \
            ${dir_freesurfer}/surf/${hemisphere}h.pial \
            ${dir_freesurfer}/surf/${hemisphere}h.sphere.reg \
            ${util_surface}/fs_LR-deformed_to-fsaverage.${HEMI}.sphere.32k_fs_LR.surf.gii \
            ${dir_freesurfer}/surf/${hemisphere}h.midthickness.surf.gii \
            ${dir_conte69}/${id}_${hemisphere}h_midthickness_32k_fs_LR.surf.gii \
            ${dir_conte69}/${id}_${hemisphere}h_sphereReg.surf.gii
        # Resample white and pial surfaces to conte69-32k
        for surface in pial white; do
            Do_cmd mris_convert ${dir_freesurfer}/surf/${hemisphere}h.${surface} ${dir_conte69}/${hemisphere}h.${surface}.surf.gii
            Do_cmd wb_command -surface-resample \
                ${dir_conte69}/${hemisphere}h.${surface}.surf.gii \
                ${dir_conte69}/${id}_${hemisphere}h_sphereReg.surf.gii \
                ${util_surface}/fs_LR-deformed_to-fsaverage.${HEMI}.sphere.32k_fs_LR.surf.gii \
                BARYCENTRIC \
                ${dir_conte69}/${id}_${hemisphere}h_${surface}_32k_fs_LR.surf.gii
        done
    done
fi


# -----------------------------------------------------------------------------------------------
# Clean temporal directory and temporal fsaverage5
Do_cmd rm -rfv $tmp  ${dir_surf}/fsaverage5

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Post-structural processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\t\t\tlogs:${dir_logs}/post_structural.txt"
echo "${id}, post_structural, CONTE69, $(date), `printf "%0.3f\n" ${eri}`" >> ${out}/brain-proc.csv
