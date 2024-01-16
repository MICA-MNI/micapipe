#!/bin/bash
#
# Superficial White Matter generation
#
# Generates vertexwise (native, fsa5, and conte69) outputs for:
#   all the qMRI (nii.gz) in /maps
#
# This workflow makes use of freesurfer outputs and custom python scripts
#
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
    #MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
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
Title "Superficial White Matter Mapping (SWM)"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_swm_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# create a json QC file
module_json="${dir_QC}/${idBIDS}_module-SWM.json"

# Create the working directory
mkdir -p "${tmp}/"

#------------------------------------------------------------------------------#
# Laplacian surface generation
Nwm=$(ls "${dir_conte69}/${idBIDS}_hemi-"*_surf-fsnative_label-swm*.surf.gii 2>/dev/null | wc -l)
if [[ "$Nwm" -lt 6 ]]; then ((N++))
    # Import the surface segmentation to NIFTI
    T1fs_seg="${tmp}/aparc+aseg.nii.gz"
    Do_cmd mri_convert "${dir_subjsurf}/mri/aparc+aseg.mgz" "${T1fs_seg}"

    # Move the segmentation to T1_nativepro space
    mat_fsnative_affine="${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_T1w_"
    T1_fsnative_affine="${mat_fsnative_affine}0GenericAffine.mat"
    T1nativepro_seg="${tmp}/aparc+aseg_space-nativepro.nii.gz"
    Do_cmd antsApplyTransforms -d 3 -i "${T1fs_seg}" -r "${T1nativepro}" -t "${T1_fsnative_affine}" -o "${T1nativepro_seg}" -n GenericLabel -v -u int

    # Generate the laplacian field
    WM_laplace=${tmp}/wm-laplace.nii.gz
    Do_cmd python "$MICAPIPE/functions/laplace_solver.py" "${T1nativepro_seg}" "${WM_laplace}"

    # Create the surfaces at 01 02 and 03 deeps
    for HEMI in L R; do
      # Prepare the white matter surface
      Do_cmd cp "${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-white.surf.gii ${tmp}/${HEMI}_wm.surf.gii"
      # Run SWM
      Do_cmd python "${MICAPIPE}"/functions/surface_generator.py "${tmp}/${HEMI}_wm.surf.gii" "${WM_laplace}" "${dir_conte69}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-swm" "1,2,3"
    done
    Nwm=$(ls "${dir_conte69}/${idBIDS}_hemi-"*_surf-fsnative_label-swm*.surf.gii 2>/dev/null | wc -l)
    if [[ "$Nwm" -ge 6 ]]; then ((Nsteps++)); fi
else
    Info "Subject ${idBIDS} has SWM surfaces"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# Mapping files
maps=(${dir_maps}/*nii*)
for map in ${maps[*]}; do
  map_id=$(echo "${map/.nii.gz/}" | awk -F 'map-' '{print $2}')
  # Map to surface: swm
      for HEMI in L R; do
          for i in $(ls "${dir_conte69}/${idBIDS}_hemi-L"_surf-fsnative_label-swm*mm.surf.gii); do
              label=$(echo "${i/.surf.gii/}" | awk -F 'label-' '{print $2}')
              Info "Mapping ${map_id} SWM-${label} to fsLR-32k, fsLR-5k and fsaverage5"
              surf_fsnative="${dir_conte69}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-${label}.surf.gii"
              # MAPPING metric to surfaces
              map_to-surfaces "${map}" "${surf_fsnative}" "${dir_maps}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-${label}_${map_id}.func.gii" "${HEMI}" "${label}_${map_id}" "${dir_maps}"
          done
      done
done

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status "SWM"
micapipe_procStatus "${id}" "${SES/ses-/}" "SWM" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "SWM" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
