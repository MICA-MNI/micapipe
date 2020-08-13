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

# Test inputs: Nativepro T1
if [ ! -f ${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz ]
# Test inputs: freesurfer-orig
if [ ! -f ${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz ]
# Test inputs: 5TT
fiveTT
if [ "$N" -lt 1 ]; then Error "Subject $id doesn't have T1 on: \n\t\t\t${subject_bids}/anat"; exit; fi


#------------------------------------------------------------------------------#
Title "Running MICA post-structural processing"


# Subjects dir for freesurfer
export SUBJECTS_DIR=/
# Sets wb_command to only use one thread
export OMP_NUM_THREADS=1

# Set basic parameters.
subject=$1
nativeT1w=$2
fiveTT=$3
mmResolution=$4
fsDirectory=$5
conte69Directory=$6
surfToVolumeDirectory=$7
volumeDirectory=$8
warpDirectory=$9
parcDirectory=${10}
templateDirectory=${11}
subcorticalSeg=${12}

processingDirectory=$(setProcessingDirectory SETDEFAULT)

# TRAP: Remove temporary directory
function cleanup() {
    find $processingDirectory -maxdepth 1 -type f -name "*.nii.gz" -delete
    find $processingDirectory -empty -delete
}
trap cleanup EXIT

# Make directories
for x in "$surfToVolumeDirectory" "$conte69Directory" "$processingDirectory" "$warpDirectory"; do
    [[ ! -d "$x" ]] && mkdir -p "$x"
done

# Compute warp of native structural to Freesurfer and apply to 5TT and first
export SUBJECTS_DIR="$fsDirectory"/../
Do_cmd bbregister --mov "$nativeT1w" --s "$subject" --reg "$warpDirectory"/$(basename $nativeT1w .nii.gz)_t1w2fs.lta --init-fsl --t1 --o "$volumeDirectory"/"$subject"_t1w_"$mmResolution"mm_fsspace.nii.gz
Do_cmd lta_convert --inlta "$warpDirectory"/$(basename $nativeT1w .nii.gz)_t1w2fs.lta --outlta "$warpDirectory"/$(basename $nativeT1w .nii.gz)_fs2t1w.lta --invert
Do_cmd mrconvert $fiveTT ${fiveTT/.mif.gz/.nii.gz}
Do_cmd mri_vol2vol --mov ${fiveTT/.mif.gz/.nii.gz} --targ "$volumeDirectory"/"$subject"_t1w_"$mmResolution"mm_fsspace.nii.gz --o $(echo $fiveTT | sed 's:.mif.gz:.nii.gz:' | sed 's:native:fsspace:') --lta "$warpDirectory"/$(basename $nativeT1w .nii.gz)_t1w2fs.lta
Do_cmd mri_vol2vol --mov ${subcorticalSeg} --targ "$volumeDirectory"/"$subject"_t1w_"$mmResolution"mm_fsspace.nii.gz --o $processingDirectory/first_fs.nii.gz --lta "$warpDirectory"/$(basename $nativeT1w .nii.gz)_t1w2fs.lta --interp nearest

## Deal with FreeSurfer c_ras offset
MatrixX=$(mri_info --cras "$fsDirectory"/mri/orig.mgz | cut -f1 -d' ')
MatrixY=$(mri_info --cras "$fsDirectory"/mri/orig.mgz | cut -f2 -d' ')
MatrixZ=$(mri_info --cras "$fsDirectory"/mri/orig.mgz | cut -f3 -d' ')
echo "1 0 0 ""$MatrixX" > "$fsDirectory"/mri/c_ras.mat
echo "0 1 0 ""$MatrixY" >> "$fsDirectory"/mri/c_ras.mat
echo "0 0 1 ""$MatrixZ" >> "$fsDirectory"/mri/c_ras.mat
echo "0 0 0 1" >> "$fsDirectory"/mri/c_ras.mat

if [[ ! -f  "$conte69Directory"/"$subject"_rh_midthickness_32k_fs_LR_fsspace_cras_corrected.surf.gii ]] ; then
    for hemisphere in l r; do
        # Build the conte69-32k sphere and midthickness surface
        Do_cmd wb_shortcuts -freesurfer-resample-prep \
            "$fsDirectory"/surf/"$hemisphere"h.white \
            "$fsDirectory"/surf/"$hemisphere"h.pial \
            "$fsDirectory"/surf/"$hemisphere"h.sphere.reg \
            "$templateDirectory"/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
            "$fsDirectory"/surf/"$hemisphere"h.midthickness.surf.gii \
            "$conte69Directory"/"$subject"_"$hemisphere"h_midthickness_32k_fs_LR.surf.gii \
            "$conte69Directory"/"$subject"_"$hemisphere"h_sphereReg.surf.gii
        # Resample white and pial surfaces to conte69-32k
        for surface in pial white; do
            Do_cmd mris_convert "$fsDirectory"/surf/"$hemisphere"h."$surface" "$conte69Directory"/"$hemisphere"h."$surface".surf.gii
            Do_cmd wb_command -surface-resample \
                "$conte69Directory"/"$hemisphere"h."$surface".surf.gii \
                "$conte69Directory"/"$subject"_"$hemisphere"h_sphereReg.surf.gii \
                "$templateDirectory"/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
                BARYCENTRIC \
                "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR.surf.gii
        done
        for surface in pial white midthickness; do
            # Apply affine transformation to pial, white, and midthickness surfaces to bring them in line with native scans.
            Do_cmd wb_command -surface-apply-affine "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR.surf.gii \
                "$fsDirectory"/mri/c_ras.mat \
                "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR_fsspace_cras_corrected.surf.gii
        done
    done
fi

## Glasser and Vos de Wael areas to volume.
for area in glasser_360 vosdewael_400 vosdewael_300 vosdewael_200 vosdewael_100 schaefer_100 schaefer_200 schaefer_300 schaefer_400 schaefer_500 schaefer_600 schaefer_800 schaefer_1000 aparc_aseg economo; do
    if [[ ! -e "$surfToVolumeDirectory"/"$subject"_"$area"_subcortex_"$mmResolution"_native.nii.gz ]]; then
        # Fill the ribbon.
        for side in lh rh; do
            Do_cmd wb_command -label-to-volume-mapping \
                "$parcDirectory"/"$area"_conte69_"$side".label.gii \
                "$conte69Directory"/*_"$side"_midthickness_32k_fs_LR_fsspace_cras_corrected.surf.gii \
                "$volumeDirectory"/"$subject"_t1w_"$mmResolution"mm_fsspace.nii.gz \
                "$processingDirectory"/"$area"_"$side".nii.gz \
                -ribbon-constrained \
                "$conte69Directory"/*_"$side"_white_32k_fs_LR_fsspace_cras_corrected.surf.gii \
                "$conte69Directory"/*_"$side"_pial_32k_fs_LR_fsspace_cras_corrected.surf.gii
        done
        # Combine left and right hemispheres by adding them EXCLUDING any overlapping voxels.
        [[ $area == glasser_360 ]] && wb_command -volume-math '((r>0)*180+r+l)*!(r>0&&l>0)' "$surfToVolumeDirectory"/"$subject"_"$area"_"$mmResolution"_native.nii.gz -var l "$processingDirectory"/"$area"_lh.nii.gz -var r "$processingDirectory"/"$area"_rh.nii.gz \
            || wb_command -volume-math '(r+l)*!(r>0&&l>0)' "$surfToVolumeDirectory"/"$subject"_"$area"_"$mmResolution"_native.nii.gz -var l "$processingDirectory"/"$area"_lh.nii.gz -var r "$processingDirectory"/"$area"_rh.nii.gz

        # Create a subcortical atlas as well. Assign overlap to the cortex.
        maxval=$(fsl5.0-fslstats "$surfToVolumeDirectory"/"$subject"_"$area"_"$mmResolution"_native.nii.gz -p 100 | sed 's:\..*::')
        wb_command -volume-math "(h+(s+${maxval})*(h==0&&s!=0))" "$surfToVolumeDirectory"/"$subject"_"$area"_subcortex_"$mmResolution"_native.nii.gz -var h "$surfToVolumeDirectory"/"$subject"_"$area"_"$mmResolution"_native.nii.gz -var s $processingDirectory/first_fs.nii.gz
    fi
done

# Notification of completition
Title "Post-structural processing ended:\n\t\t\tlogs:${dir_logs}/post_structural.txt"
echo "$"
