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

# Test inputs: Nativepro T1
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi

# Test inputs: freesurfer-orig
if [ ! -f ${T15ttgen} ]; then Error "Subject $id doesn't have T1_tt5"; exit; fi

# Test inputs: 5TT
if [ ! -f ${T1fast_all} ]; then Error "Subject $id doesn't have T1_fast segmentation"; exit; fi

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
T1_fsspace=${proc_struct}/${id}_t1w_${res}mm_fsspace.nii.gz
T1natpro_str=$(basename $T1nativepro .nii.gz)

Do_cmd bbregister --mov $T1nativepro --s $id --reg ${dir_warp}/${T1natpro_str}_t1w2fs.lta --init-fsl --t1 --o ${T1_fsspace}
Do_cmd lta_convert --inlta ${dir_warp}/${T1natpro_str}_t1w2fs.lta --outlta ${dir_warp}/${T1natpro_str}_fs2t1w.lta --invert

Do_cmd mri_vol2vol --mov ${T15ttgen} --targ ${T1_fsspace} --o ${T15ttgen/nativepro/fsspace} --lta ${dir_warp}/${T1natpro_str}_t1w2fs.lta
Do_cmd mri_vol2vol --mov ${T1fast_seg} --targ ${T1_fsspace} --o ${tmp}/first_fs.nii.gz --lta ${dir_warp}/${T1natpro_str}_t1w2fs.lta --interp nearest

## Deal with FreeSurfer c_ras offset
MatrixX=$(mri_info --cras ${dir_freesurfer}/mri/orig.mgz | cut -f1 -d' ')
MatrixY=$(mri_info --cras ${dir_freesurfer}/mri/orig.mgz | cut -f2 -d' ')
MatrixZ=$(mri_info --cras ${dir_freesurfer}/mri/orig.mgz | cut -f3 -d' ')
echo "1 0 0 ""$MatrixX" > ${dir_freesurfer}/mri/c_ras.mat
echo "0 1 0 ""$MatrixY" >> ${dir_freesurfer}/mri/c_ras.mat
echo "0 0 1 ""$MatrixZ" >> ${dir_freesurfer}/mri/c_ras.mat
echo "0 0 0 1" >> ${dir_freesurfer}/mri/c_ras.mat

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
        for surface in pial white midthickness; do
            # Apply affine transformation to pial, white, and midthickness surfaces to bring them in line with native scans.
            Do_cmd wb_command -surface-apply-affine ${dir_conte69}/${id}_${hemisphere}h_${surface}_32k_fs_LR.surf.gii \
                ${dir_freesurfer}/mri/c_ras.mat \
                ${dir_conte69}/${id}_${hemisphere}h_${surface}_32k_fs_LR_fsspace_cras_corrected.surf.gii
        done
    done
fi

## Glasser and Vos de Wael areas to volume.
for area in glasser_360 vosdewael_400 vosdewael_300 vosdewael_200 vosdewael_100 schaefer_100 schaefer_200 schaefer_300 schaefer_400 schaefer_500 schaefer_600 schaefer_800 schaefer_1000 aparc_aseg economo; do
    if [[ ! -e ${dir_volum}/${id}_${area}_subcortex_${res}_native.nii.gz ]]; then
        # Fill the ribbon.
        for side in lh rh; do
            Do_cmd wb_command -label-to-volume-mapping \
                ${util_parcelations}/${area}_conte69_${side}.label.gii \
                ${dir_conte69}/*_${side}_midthickness_32k_fs_LR_fsspace_cras_corrected.surf.gii \
                ${T1_fsspace} \
                ${tmp}/${area}_${side}.nii.gz \
                -ribbon-constrained \
                ${dir_conte69}/*_${side}_white_32k_fs_LR_fsspace_cras_corrected.surf.gii \
                ${dir_conte69}/*_${side}_pial_32k_fs_LR_fsspace_cras_corrected.surf.gii
        done
        # Combine left and right hemispheres by adding them EXCLUDING any overlapping voxels.
        [[ $area == glasser_360 ]] && wb_command -volume-math '((r>0)*180+r+l)*!(r>0&&l>0)' ${dir_volum}/${id}_${area}_${res}_native.nii.gz -var l ${tmp}/${area}_lh.nii.gz -var r ${tmp}/${area}_rh.nii.gz \
            || wb_command -volume-math '(r+l)*!(r>0&&l>0)' ${dir_volum}/${id}_${area}_${res}_native.nii.gz -var l ${tmp}/${area}_lh.nii.gz -var r ${tmp}/${area}_rh.nii.gz

        # Create a subcortical atlas as well. Assign overlap to the cortex.
        maxval=$(fsl5.0-fslstats ${dir_volum}/${id}_${area}_${res}_native.nii.gz -p 100 | sed 's:\..*::')
        wb_command -volume-math "(h+(s+${maxval})*(h==0&&s!=0))" ${dir_volum}/${id}_${area}_subcortex_${res}_native.nii.gz -var h ${dir_volum}/${id}_${area}_${res}_native.nii.gz -var s ${tmp}/first_fs.nii.gz
    fi
done


# Clean temporal directory
Do_cmd rm -rfv $tmp

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Post-structural processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\t\t\tlogs:${dir_logs}/post_structural.txt"

# echo "${id}, post_structural, DONE, $(date), `printf "%0.3f\n" ${eri}`" >> ${out}/brain-proc.csv
