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
#
umask 003
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
atlas=$8
PROC=$9
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Manage manual inputs: Parcellations
cd "$util_parcelations"
if [[ "$atlas" == "DEFAULT" ]]; then
  atlas_parc=($(ls lh.*annot))
  N="${#atlas_parc[*]}"
  Info "Selected parcellations: DEFAULT, N=${N}"
else
  IFS=',' read -ra atlas_parc <<< "$atlas"
  for i in "${!atlas_parc[@]}"; do atlas_parc[i]=$(ls lh."${atlas_parc[$i]}"_mics.annot 2>/dev/null); done
  atlas_parc=("${atlas_parc[@]}")
  N="${#atlas_parc[*]}"
  Info "Selected parcellations: $atlas, N=${N}"
fi
cd "$here"

# Check inputs: Nativepro T1
if [ "$N" -eq 0 ]; then
  Error "Provided -atlas do not match with any on MICAPIPE, try one of the following list:
\t\taparc-a2009s,aparc,economo,glasser-360
\t\tschaefer-1000,schaefer-100,schaefer-200
\t\tschaefer-300,schaefer-400,schaefer-500
\t\tschaefer-600,schaefer-700,schaefer-800
\t\tschaefer-900,vosdewael-100,vosdewael-200
\t\tvosdewael-300,vosdewael-400"; exit
fi

# Check inputs: Nativepro T1
if [ ! -f "${proc_struct}/${idBIDS}"_space-nativepro_t1w.nii.gz ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi
if [ ! -f "$T1fast_seg" ]; then Error "Subject $id doesn't have FAST: run -proc_structural"; exit; fi
# Check inputs: freesurfer space T1
if [ ! -f "$T1freesurfr" ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${idBIDS}/mri/T1.mgz"; exit; fi

#------------------------------------------------------------------------------#
Title "POST-structural processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
# print the names on the terminal
bids_print.variables-post

# GLOBAL variables for this script
Info "Saving temporal dir: $nocleanup"
Info "ANTs will use $threads threads"
Info "wb_command will use $OMP_NUM_THREADS threads"

#	Timer
aloita=$(date +%s)
Nsteps=0

# Create script specific temp directory
tmp=${tmpDir}/${RANDOM}_micapipe_post-struct_${idBIDS}
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Freesurface SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

#------------------------------------------------------------------------------#
# Compute affine matrix from Freesurfer space to nativepro
T1_in_fs=${tmp}/T1.nii.gz
T1_fsnative=${proc_struct}/${idBIDS}_space-fsnative_t1w.nii.gz
mat_fsnative_affine=${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_t1w_
T1_fsnative_affine=${mat_fsnative_affine}0GenericAffine.mat

if [[ ! -f "$T1_fsnative" ]] || [[ ! -f "$T1_fsnative_affine" ]]; then
    Do_cmd mrconvert "$T1freesurfr" "$T1_in_fs"
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1nativepro" -m "$T1_in_fs" -o "$mat_fsnative_affine" -t a -n "$threads" -p d
    Do_cmd antsApplyTransforms -d 3 -i "$T1nativepro" -r "$T1_in_fs" -t ["${T1_fsnative_affine}",1] -o "$T1_fsnative" -v -u int
    if [[ -f ${T1_fsnative} ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a T1 on FreeSurfer space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Cerebellar parcellation from MNI152_0.8mm to T1-Nativepro
Info "Cerebellum parcellation to T1-nativepro Volume"
# Variables
T1str_nat="${idBIDS}_space-nativepro_t1w_atlas"
T1str_fs="${idBIDS}_space-fsnative_t1w"
atlas_cerebellum="${util_MNIvolumes}/MNI152_T1_0.8mm_cerebellum.nii.gz"       # cerebellar lobules atlas
mat_MNI152_SyN="${dir_warp}/${idBIDS}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_"    # transformation strings nativepro to MNI152_0.8mm
T1_MNI152_InvWarp="${mat_MNI152_SyN}1InverseWarp.nii.gz"                      # Inversewarp - nativepro to MNI152_0.8mm
T1_MNI152_affine="${mat_MNI152_SyN}0GenericAffine.mat"                        # Affine matrix - nativepro to MNI152_0.8mm
T1_seg_cerebellum="${dir_volum}/${T1str_nat}-cerebellum.nii.gz"               # Cerebellar output
T1_seg_subcortex="${dir_volum}/${T1str_nat}-subcortical.nii.gz"

# Apply inverse transfrom from MNI152-cerebellum to T1-nativepro
if [[ ! -f "$T1_seg_cerebellum" ]]; then
    Do_cmd antsApplyTransforms -d 3 -i "$atlas_cerebellum" \
                -r "$T1nativepro" \
                -n GenericLabel -t ["$T1_MNI152_affine",1] -t "$T1_MNI152_InvWarp" \
                -o "$T1_seg_cerebellum" -v -u int
    if [[ -f "$T1_seg_cerebellum" ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a Cerebellum parcellation on T1-nativepro"; ((Nsteps++))
fi

Info "Subcortical parcellation to T1-nativepro Volume"
if [[ ! -f "$T1_seg_subcortex" ]]; then
    Do_cmd cp "$T1fast_seg" "$T1_seg_subcortex"
    if [[ -f "$T1_seg_subcortex" ]]; then ((Nsteps++)); fi
    Do_cmd rm -rf "$proc_struct/first/"
else
    Info "Subject ${id} has a Subcortical parcellation on T1-nativepro"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Create parcellation on nativepro space
Info "fsaverage5 annnot parcellations to T1-nativepro Volume"
cd "$util_parcelations"
for parc in "${atlas_parc[@]}"; do
    parc_annot="${parc/lh./}"
    parc_str=$(echo "${parc_annot}" | awk -F '_mics' '{print $1}')
    if [[ ! -f "${dir_volum}/${T1str_nat}-${parc_str}.nii.gz" ]]; then
        for hemi in lh rh; do
        Info "Running surface $hemi $parc_annot to $subject"
        Do_cmd mri_surf2surf --hemi "$hemi" \
               		  --srcsubject fsaverage5 \
               		  --trgsubject "$idBIDS" \
               		  --sval-annot "${hemi}.${parc_annot}" \
               		  --tval "${dir_freesurfer}/label/${hemi}.${parc_annot}"
        done
        fs_mgz="${tmp}/${parc_str}.mgz"
        fs_tmp="${tmp}/${parc_str}_in_T1.mgz"
        fs_nii="${tmp}/${T1str_fs}_${parc_str}.nii.gz"                   # labels in fsnative tmp dir
        labels_nativepro="${dir_volum}/${T1str_nat}-${parc_str}.nii.gz"  # lables in nativepro

        # Register the annot surface parcelation to the T1-freesurfer volume
        Do_cmd mri_aparc2aseg --s "$idBIDS" --o "$fs_mgz" --annot "${parc_annot/.annot/}" --new-ribbon
        Do_cmd mri_label2vol --seg "$fs_mgz" --temp "$T1freesurfr" --o "$fs_tmp" --regheader "${dir_freesurfer}/mri/aseg.mgz"
        Do_cmd mrconvert "$fs_tmp" "$fs_nii" -force      # mgz to nifti_gz
        Do_cmd fslreorient2std "$fs_nii" "$fs_nii"       # reorient to standard
        Do_cmd fslmaths "$fs_nii" -thr 1000 "$fs_nii"    # threshold the labels

        # Register parcellation to nativepro
        Do_cmd antsApplyTransforms -d 3 -i "$fs_nii" -r "$T1nativepro" -n GenericLabel -t "$T1_fsnative_affine" -o "$labels_nativepro" -v -u int
        if [[ -f "$labels_nativepro" ]]; then ((Nsteps++)); fi
    else
        Info "Subject ${id} has a ${parc_str} segmentation on T1-nativepro space"
        ((Nsteps++))
    fi
done
Do_cmd rm -rf ${dir_warp}/*Warped.nii.gz 2>/dev/null

#------------------------------------------------------------------------------#
# Compute warp of native structural to Freesurfer and apply to 5TT and first
if [[ ! -f "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-rh_midthickness.surf.gii" ]]; then
    for hemisphere in l r; do
      Info "Native surfaces to conte69-64k vertices (${hemisphere}h hemisphere)"
      HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
        # Build the conte69-32k sphere and midthickness surface
        Do_cmd wb_shortcuts -freesurfer-resample-prep \
            "${dir_freesurfer}/surf/${hemisphere}h.white" \
            "${dir_freesurfer}/surf/${hemisphere}h.pial" \
            "${dir_freesurfer}/surf/${hemisphere}h.sphere.reg" \
            "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
            "${dir_freesurfer}/surf/${hemisphere}h.midthickness.surf.gii" \
            "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemisphere}h_midthickness.surf.gii" \
            "${dir_conte69}/${idBIDS}_${hemisphere}h_sphereReg.surf.gii"
        # Resample white and pial surfaces to conte69-32k
        for surface in pial white; do
            Do_cmd mris_convert "${dir_freesurfer}/surf/${hemisphere}h.${surface}" "${tmp}/${hemisphere}h.${surface}.surf.gii"
            Do_cmd wb_command -surface-resample \
                "${tmp}/${hemisphere}h.${surface}.surf.gii" \
                "${dir_conte69}/${idBIDS}_${hemisphere}h_sphereReg.surf.gii" \
                "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
                BARYCENTRIC \
                "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemisphere}h_${surface}.surf.gii"
            if [[ -f "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemisphere}h_${surface}.surf.gii" ]]; then ((Nsteps++)); fi
        done
    done
else
    Info "Subject ${idBIDS} has surfaces on conte69"; Nsteps=$((Nsteps+4))
fi

# -----------------------------------------------------------------------------------------------
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
N=$((N+7)) # total number of steps
if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "Post-structural processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\tSteps completed : $(printf "%02d" $Nsteps)/$(printf "%02d" $N)
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/post_structural_*.txt)"
# Print QC stamp
micapipe_procStatus "${id}" "${SES/ses-/}" "post_structural" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "post_structural" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
cleanup "$tmp" "$nocleanup" "$here"
