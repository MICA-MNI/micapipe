#!/bin/bash
#
# DWI structural processing with bash:
#
# Preprocessing workflow for diffusion MRI.
#
# This workflow makes use of freesurfer, FSL workbench and MRtrix3
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
# TEST=ON
# source utilities
source $MICAPIPE/functions/utilities.sh

BIDS=$1
id=$2
out=$3
tmp=$4

# CORES should be defined as GLOBAL variable
if [[ -z $CORES ]]; then CORES=6; Info "ANTs will use $CORES CORES"; fi

# Sets wb_command to only use one thread
if [[ -z $OMP_NUM_THREADS ]]; then OMP_NUM_THREADS=4; Info "wb_command will use $OMP_NUM_THREADS threads"; fi

# Assigns variables names
bids_variables $BIDS $id $out
# print the names on the terminal
bids_print.variables-post

# Check inputs: Nativepro T1
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi

# Check inputs: freesurfer space T1
if [ ! -f ${T1freesurfr} ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi


#------------------------------------------------------------------------------#
Title "Running MICA POST-structural processing"

#	Timer
aloita=$(date +%s)
here=`pwd`

# Check tmp dir: temporary directory
random_str=$RANDOM
if [ -z ${tmp} ]; then tmp=/tmp/${random_str}_post_structural_${id}; fi
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Freesurface SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

#------------------------------------------------------------------------------#
# Compute affine matrix from Freesurfer space to nativepro
T1_in_fs=${tmp}/T1.nii.gz
T1_fsspace=${proc_struct}/${id}_t1w_${res}mm_fsspace.nii.gz
mat_fsspace_affine=${dir_warp}/${id}_t1w_${res}mm_fsspace_to_nativepro_
T1_fsspace_affine=${mat_fsspace_affine}0GenericAffine.mat

Do_cmd mrconvert $T1freesurfr $T1_in_fs
Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro -m $T1_in_fs -o $mat_fsspace_affine -t a -n $CORES -p d
Do_cmd antsApplyTransforms -d 3 -i $T1nativepro -r $T1_in_fs -t [${T1_fsspace_affine},1] -o $T1_fsspace -v -u int

#------------------------------------------------------------------------------#
# Cerebellar parcellation from MNI152_0.8mm to T1-Nativepro
Info "Cerebellum parcellation to T1-nativepro Volume"
# Variables
T1str_nat=${id}_t1w_${res}mm_nativepro
T1str_fs=${id}_t1w_${res}mm_fsspace
atlas_cerebellum=${util_MNIvolumes}/MNI152_T1_0.8mm_cerebellum.nii.gz       # cerebellar lobules atlas
mat_MNI152_SyN=${dir_warp}/${T1str_nat}_brain_to_0.8mm_MNI152_SyN_brain_    # transformation strings nativepro to MNI152_0.8mm
T1_MNI152_InvWarp=${mat_MNI152_SyN}1InverseWarp.nii.gz                      # Inversewarp - nativepro to MNI152_0.8mm
T1_MNI152_affine=${mat_MNI152_SyN}0GenericAffine.mat                        # Affine matrix - nativepro to MNI152_0.8mm
T1_seg_cerebellum=${dir_volum}/${T1str_nat}_cerebellum.nii.gz               # Cerebellar output

# Apply inverse transfrom from MNI152-cerebellum to T1-nativepro
Do_cmd antsApplyTransforms -d 3 -i $atlas_cerebellum \
              -r ${T1nativepro} \
              -n GenericLabel -t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp} \
              -o ${T1_seg_cerebellum} -v -u int

Info "Subcortical parcellation to T1-nativepro Volume"
Do_cmd cp ${T1fast_all} ${dir_volum}/${T1str_nat}_subcortical.nii.gz

#------------------------------------------------------------------------------#
# Create parcellation on nativepro space
Info "fsaverage5 annnot parcellations to T1-nativepro Volume"
Do_cmd cp -vR ${util_surface}/fsaverage5 ${dir_surf}
cd $util_parcelations
for parc in lh.*.annot; do
   parc_annot=${parc/lh./}
   parc_str=`echo ${parc_annot} | awk -F '_fsa5' '{print $1}'`
 	 for hemi in lh rh; do
   		Info "Running surface $hemi $parc_annot to $subject"
   		Do_cmd mri_surf2surf --hemi $hemi \
   		  --srcsubject fsaverage5 \
   		  --trgsubject ${id} \
   		  --sval-annot ${hemi}.${parc_annot} \
   		  --tval ${dir_freesurfer}/label/${hemi}.${parc_annot}
 	 done
   fs_mgz=${tmp}/${parc_str}.mgz
   fs_tmp=${tmp}/${parc_str}_in_T1.mgz
   fs_nii=${tmp}/${T1str_fs}_${parc_str}.nii.gz                   # labels in fsspace tmp dir
   labels_nativepro=${dir_volum}/${T1str_nat}_${parc_str}.nii.gz  # lables in nativepro

   # Register the annot surface parcelation to the T1-freesurfer volume
   Do_cmd mri_aparc2aseg --s ${id} --o ${fs_mgz} --annot ${parc_annot/.annot/} --new-ribbon
   Do_cmd mri_label2vol --seg ${fs_mgz} --temp ${dir_freesurfer}/mri/T1.mgz --o $fs_tmp --regheader ${dir_freesurfer}/mri/aseg.mgz
   Do_cmd mrconvert $fs_tmp $fs_nii -force      # mgz to nifti_gz
   Do_cmd fslreorient2std $fs_nii $fs_nii       # reorient to standard
   Do_cmd fslmaths $fs_nii -thr 1000 $fs_nii    # threshold the labels
   # Register parcellation to nativepro
   Do_cmd antsApplyTransforms -d 3 -i $fs_nii -r $T1nativepro -n GenericLabel -t $T1_fsspace_affine -o $labels_nativepro -v -u int
done
cd $here

#------------------------------------------------------------------------------#
# Compute warp of native structural to Freesurfer and apply to 5TT and first
Info "Native surfaces to conte69-64k vertices (both hemispheres)"
if [[ ! -f  ${dir_conte69}/${id}_rh_midthickness_32k_fs_LR.surf.gii ]] ; then
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
echo "${id}, post_structural, TEST, $(date), `printf "%0.3f\n" ${eri}`" >> ${out}/brain-proc.csv
