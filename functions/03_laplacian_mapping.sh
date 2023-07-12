#!/bin/bash
#
# Cortical morphology metrics processing:
#
# Generates vertexwise (native, fsa5, and conte69) outputs for:
#   Cortical Thickness
#   Mean Curvature
#
# This workflow makes use of freesurfer outputs and custom python scripts
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micapipe/tree/master/parcellations
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out Directory
#
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
PROC=$8
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE"/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfRecon "${post_struct_json}" | awk -F '"' '{print $4}')
set_surface_directory "${recon}"

#------------------------------------------------------------------------------#
Title "Subcortical White Matter Mapping (SWM)"

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_post-morpho_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Data location
dataDir="${dir_subjsurf}/surf"

#------------------------------------------------------------------------------#
# Lapplacian mapping and surface generation

mkdir -p "${tmp}/"


# Apply transformation from surface space to nativepro space
Do_cmd mri_convert "${dir_subjsurf}/mri/aparc+aseg.mgz" "${tmp}/aparc+aseg.nii.gz"

mat_fsnative_affine=${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_T1w_
T1_fsnative_affine=${mat_fsnative_affine}0GenericAffine.mat
Do_cmd antsApplyTransforms -d 3 -i "${tmp}/aparc+aseg.nii.gz" -r "$T1nativepro" -n MultiLabel -t "$T1_fsnative_affine" -o "${tmp}/aparc+aseg_space-nativepro.nii.gz" -u int

# Solve a Laplace field
Do_cmd python "$MICAPIPE"/functions/laplace_solver.py "${tmp}/aparc+aseg_space-nativepro.nii.gz" "${tmp}/wm-laplace.nii.gz"
# Shift a given surface along the Laplace field
mkdir -p "${dir_conte69}/wm-equipotentials"
for depths in 0.05 0.1 0.2 0.3; do
    Do_cmd python "$MICAPIPE"/functions/laplace_surf_interp.py "${dir_conte69}/${idBIDS}_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii" "${tmp}/wm-laplace.nii.gz" "${tmp}/""${dir_conte69}/wm-equipotentials/{idBIDS}_hemi-L_surfdepth-${depths}" $depths
    Do_cmd python "$MICAPIPE"/functions/laplace_surf_interp.py "${dir_conte69}/${idBIDS}_hemi-R_space-nativepro_surf-fsLR-32k_label-white.surf.gii" "${tmp}/wm-laplace.nii.gz" "${tmp}/""${dir_conte69}/wm-equipotentials/{idBIDS}_hemi-R_surfdepth-${depths}" $depths
done

exit
#------------------------------------------------------------------------------#
# End if module has been processed
cleanup "$tmp" "$nocleanup" "$here"
