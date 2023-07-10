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

# Make output directory
outDir="${dir_maps}"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_subjsurf}/surf"

#------------------------------------------------------------------------------#
# ImageMath "${tmp}/ventricles.nii.gz" Binarize "${dir_surf}/mri/aparc+aseg.mgz" 4 43 31 63

mkdir -p "${tmp}/"
mkdir "${proc_struct}/swm/pericort/"

# Get maps
Do_cmd mri_binarize --i "${dir_surf}/mri/aparc+aseg.mgz" --match 4 --match 43 --match 31 --match 63 \
                        --o "${tmp}/ventricles.nii.gz"
Do_cmd mri_binarize --i "${dir_surf}/mri/aparc+aseg.mgz" --min 1000 --max 2999 \
                        --o "${tmp}/cortex1.nii.gz"
Do_cmd mri_binarize --i "${dir_surf}/mri/aparc+aseg.mgz" --match 54 --match 18 \
                        --o "${tmp}/cortex2.nii.gz"
Do_cmd mri_binarize --i "${dir_surf}/mri/aparc+aseg.mgz" --match 2   --match 4   --match 11  --match 12  \
                        --match 26  --match 17  --match 31  --match 10 --match 5  --match 28  --match 13  --match 30 \
                        --match 41  --match 43  --match 50  --match 51  --match 58  --match 53  --match 63  --match 49 \
                        --match 44 --match 60  --match 52  --match 62 --match 77  --match 255 --match 254 --match 253 \
                        --match 252 --match 251 --match 72  --match 80 \
                        --o "${tmp}/wm.nii.gz"

Do_cmd ImageMath 3 "${tmp}/cortex.nii.gz" + "${tmp}/cortex1.nii.gz" "${tmp}/cortex2.nii.gz"

#------------------------------------------------------------------------------#
# Apply transformation from surface space to nativepro space
mat_fsnative_affine=${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_T1w_
T1_fsnative_affine=${mat_fsnative_affine}0GenericAffine.mat

for i in cortex wm ventricles; do
Do_cmd antsApplyTransforms -d 3 -i "${tmp}/${i}.nii.gz" -r "$T1nativepro" -n GenericLabel -t "$T1_fsnative_affine" -o "${tmp}/space-nativepro_${i}.nii.gz" -u int
done
#------------------------------------------------------------------------------#
# Lapplacian mapping and surface generation
exit
#------------------------------------------------------------------------------#
# End if module has been processed
cleanup "$tmp" "$nocleanup" "$here"
