#!/bin/bash
#
# Microstructural imaging processing:
#
# Preprocessing workflow for qT1.
# Generates microstructural profiles and mpc matrices on specified parcellations
#
# This workflow makes use of freesurfer and custom python scripts
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
input_im=$8
mpc_reg=$9
mpc_str=${10}
synth_reg=${11}
reg_nonlinear=${12}
PROC=${13}
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE"/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check dependencies Status: POST_STRUCTURAL
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfRecon "${post_struct_json}" | awk -F '"' '{print $4}')
set_surface_directory "${recon}"

# Variables naming for multiple acquisitions
if [[ "${mpc_str}" == DEFAULT ]]; then
  mpc_str="qMRI"
  mpc_p="acq-qMRI"
else
  mpc_p="acq-${mpc_str}"
fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-MPC-SWM-${mpc_str}.json"
micapipe_check_json_status "${module_json}" "MPC-SWM"

# Check microstructural image input flag and set parameters accordingly
if [[ "$input_im" == "DEFAULT" ]]; then microImage="$bids_T1map"; else microImage="${input_im}"; fi
Note "Microstructural image =" "$microImage"

# Check microstructural image to registrer
if [[ "$mpc_reg" == "DEFAULT" ]]; then regImage="${bids_inv1}"; else regImage="${mpc_reg}"; fi
Note "Microstructural image for registration =" "$regImage"

# Exit if microImage or Registration image do not exists
if [ ! -f "${microImage}" ]; then Error "Image for MPC-SWM was not found or the path is wrong!!!"; exit; fi
if [ ! -f "${regImage}" ]; then Error "Image for MPC-SWM registration was not found or the path is wrong!!!"; exit; fi

#------------------------------------------------------------------------------#
Title "Microstructural Profiles Covariance - SWM\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-post
Note "Saving temporal dir : " "${nocleanup}"
Note "Parallel processing : " "${threads} threads"
Note "tmp dir   : " "${tmpDir}"
Note "recon     : " "${recon}"
Note "synth_reg : " "${synth_reg}"
Note "reg_nonlinear : " "${reg_nonlinear}"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
#tmp="${tmpDir}/${RANDOM}_micapipe_${mpc_str}-MPC_${id}"
tmp="${tmpDir}/${RANDOM}_micapipe_mpc-swm_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Freesurface SUBJECTs directory
export SUBJECTS_DIR="$dir_surf"
outDir="${subject_dir}/mpc-swm/${mpc_p}"
Note "acqMRI:" "${mpc_str}"

#------------------------------------------------------------------------------#
# Registration between both images
qMRI_warped=${proc_struct}/${idBIDS}_space-nativepro_${mpc_str}.nii.gz
# This is only for the json metada
qMRI_reference="${T1nativepro}"

# Affine transformations
str_qMRI2np_xfm="${dir_warp}/${idBIDS}_from-${mpc_str}_to-nativepro_"
mat_qMRI2np_xfm="${str_qMRI2np_xfm}0GenericAffine.mat"

# SyN_transformations
SyN_qMRI2np_warp="${str_qMRI2np_xfm}1Warp.nii.gz"
SyN_qMRI2np_Invwarp="${str_qMRI2np_xfm}1InverseWarp.nii.gz"

# Apply transformations
if [[ ${reg_nonlinear}  == "TRUE" ]]; then
    # SyN from T1_nativepro to t1-nativepro
    reg="s"
    transformsInv="-t [${mat_qMRI2np_xfm},1] -t ${SyN_qMRI2np_Invwarp}" # T1nativepro to qMRI
    transforms="-t ${SyN_qMRI2np_warp} -t ${mat_qMRI2np_xfm}"  # qMRI to T1nativepro T1nativepro to qMRI
else
    reg="a"
    transformsInv="-t [${mat_qMRI2np_xfm},1]"  # T1nativepro to qMRI
    transforms="-t ${mat_qMRI2np_xfm}"   # qMRI to T1nativepro
fi

synthseg_native() {
  mri_img=$1
  mri_str=$2
  mri_synth="${tmp}/${mri_str}_synthsegGM.nii.gz"
  Do_cmd mri_synthseg --i "${mri_img}" --o "${tmp}/${mri_str}_synthseg.nii.gz" --robust --threads "$threads" --cpu
  Do_cmd fslmaths "${tmp}/${mri_str}_synthseg.nii.gz" -uthr 42 -thr 42 -bin -mul -39 -add "${tmp}/${mri_str}_synthseg.nii.gz" "${mri_synth}"
}

# Calculate the registrations
if [[ ! -f "$qMRI_warped" ]] || [[ ! -f "$mat_qMRI2np_xfm" ]]; then ((N++))
    img_fixed="${T1nativepro}"
    img_moving="${regImage}"

    # Registration with synthseg
    if [[ "${synth_reg}" == "TRUE" ]]; then
      Info "Running label based affine registrations"
      synthseg_native "${T1nativepro}" "T1w"
      synthseg_native "${regImage}" "qT1"
      img_fixed="${tmp}/T1w_synthsegGM.nii.gz"
      img_moving="${tmp}/qT1_synthsegGM.nii.gz"
    fi

    # Registrations from t1-fsnative to qMRI
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$img_fixed" -m "$img_moving" -o "$str_qMRI2np_xfm" -t "${reg}" -n "$threads" -p d -i ["${img_fixed}","${img_moving}",0]
    rm "${dir_warp}/${idBIDS}"*_Warped.nii.gz

    # Check if transformation file exist
    if [ ! -f "${mat_qMRI2np_xfm}" ]; then Error "Registration between ${mpc_str} and T1nativepro FAILED. Check you inputs!"; cleanup "$tmp" "$nocleanup" "$here"; exit; fi

    # Apply transformations: from qMRI to T1-fsnative
    Do_cmd antsApplyTransforms -d 3 -i "$microImage" -r "$T1nativepro" "${transforms}" -o "$qMRI_warped" -v -u int
    if [[ -f ${qMRI_warped} ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a ${mpc_str} on Surface space"; ((Nsteps++)); ((N++))
fi

# Convert the ANTs transformation file for wb_command
wb_affine="${tmp}/${idBIDS}_from-fsnative_to_qMRI_wb.mat"
Do_cmd c3d_affine_tool -itk "${mat_qMRI2np_xfm}" -o "${wb_affine} -inv"

#------------------------------------------------------------------------------#
# Laplacian surface generation
num_surfs=16 #15, 10
thickness=0.2 #0.2, 0.3

# Check if the directory exists and change the permissions
[[ ! -d "$outDir" ]] && mkdir -p "$outDir" && chmod -R 770 "$outDir"
# Create the MPC-SWM module json file
json_mpc "$microImage" "${outDir}/${idBIDS}_MPC-SWM-${mpc_str}.json"

# Define 16 surfaces in 3mm separated by 0.2mm step (3 mm=0.2mm step * 15 surfaces)
deepths=($(seq 0 "${thickness}" 3))

# Replace the blanck space with comma
printf -v deepths_comma '%s,' "${deepths[@]}"
deepths_comma=$(echo "${deepths_comma%,}")

# Run the SWM generation ig the fsLR-5k MPC does not exist
MPC_fsLR5k="${outDir}/${idBIDS}_surf-fsLR-5k_desc-MPC.shape.gii"
if [[ ! -f "${MPC_fsLR5k}" ]]; then ((N++))
    # Import the surface segmentation to NIFTI
    T1fs_seg="${tmp}/aparc+aseg.nii.gz"
    Do_cmd mri_convert "${dir_subjsurf}/mri/aparc+aseg.mgz" "${T1fs_seg}"

    # Move the segmentation to T1_nativepro space
    # tranformations where calculated on -post_structural
    mat_fsnative_affine="${dir_warp}/${idBIDS}_from-fsnative_to_nativepro_T1w_"
    T1_fsnative_affine="${mat_fsnative_affine}0GenericAffine.mat"
    T1nativepro_seg="${tmp}/aparc+aseg_space-nativepro.nii.gz"
    Do_cmd antsApplyTransforms -d 3 -i "${T1fs_seg}" -r "${T1nativepro}" -t "${T1_fsnative_affine}" -o "${T1nativepro_seg}" -n GenericLabel -v -u int

    # Generate the laplacian field
    WM_laplace=${tmp}/wm-laplace.nii.gz
    Info "Generating laplacian field"
    Do_cmd python "${MICAPIPE}/functions/laplace_solver.py" "${T1nativepro_seg}" "${WM_laplace}"

    # Create the surfaces by depths
    Info "Creating superficial white matter surfaces"
    for HEMI in L R; do
      # Prepare the white matter surface
      Do_cmd cp "${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-white.surf.gii ${tmp}/${HEMI}_wm.surf.gii"
      # Run SWM
      Do_cmd python "${MICAPIPE}"/functions/surface_generator.py "${tmp}/${HEMI}_wm.surf.gii" "${WM_laplace}" "${tmp}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-swm" "${deepths_comma}"

      # Copy white matter as surface 0.0
      cp "${tmp}/${HEMI}_wm.surf.gii" "${tmp}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-swm0.0mm.surf.gii"

      # find all laplacian surfaces and list by creation time
      for i in ${!deepths[@]} ; do
          mm="${deepths[$i]}mm"
          label="MPC-"$((i+1))
          # SWM surface for each deepth in NATIVEPRO
          surf_swm="${tmp}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-swm${mm}.surf.gii"
          # SWM surface for each deepth in qMRI space
          out_surf="${tmp}/${idBIDS}_hemi-${HEMI}_space-qMRI_surf-fsnative_label-swm${mm}.surf.gii"
          # SWM map for each deepth in qMRI space
          out_feat="${outDir}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-${label}.func.gii"

          Info "Creating ${HEMI} swm ${mm}"
          # Apply transformation to register surface to nativepro
          wb_command -surface-apply-affine "${surf_swm}" "${wb_affine}" "${out_surf}"
          # Apply Non-linear Warpfield to register surface to nativepro
          if [[ ${reg_nonlinear}  == "TRUE" ]]; then Do_cmd wb_command -surface-apply-warpfield "${out_surf}" "${SyN_qMRI2np_Invwarp}" "${out_surf}"; fi
          # Sample intensity and resample to other surfaces fsaverage5, fsLR-32k and fsLR-5k
          map_to-surfaces "${microImage}" "${out_surf}" "${out_feat}" "${HEMI}" "${label}" "${outDir}"
       done
    done
    Nwm=$(ls "${tmp}/${idBIDS}_hemi-"*_space-qMRI_surf-fsnative_label-swm*.surf.gii 2>/dev/null | wc -l)
    if [[ "$Nwm" -ge $((num_surfs*2)) ]]; then ((Nsteps++)); fi
else
    Info "Subject ${idBIDS} has SWM surfaces"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# Create MPC connectomes and Intensity profiles per parcellations
parcellations=($(find "$dir_volum" -name "*atlas*" ! -name "*cerebellum*" ! -name "*subcortical*"))
for seg in "${parcellations[@]}"; do
    parc=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    parc_annot="${parc}_mics.annot"
    MPC_int="${outDir}/${idBIDS}_atlas-${parc}_desc-intensity_profiles.shape.gii"
    if [[ ! -f "$MPC_int" ]]; then ((N++))
        Info "Running MPC on $parc"
        Do_cmd python "$MICAPIPE"/functions/surf2mpc.py "$out" "$id" "$SES" "$num_surfs" "$parc_annot" "$dir_subjsurf" "${mpc_p}" "mpc-swm"
        if [[ -f "$MPC_int" ]]; then ((Nsteps++)); fi
    else Info "Subject ${id} has MPC connectome and intensity profile on ${parc}"; ((Nsteps++)); ((N++)); fi
done

#------------------------------------------------------------------------------#
# Create vertex-wise MPC connectome and directory cleanup
if [[ ! -f "${MPC_fsLR5k}" ]]; then ((N++))
  Info "Running MPC vertex-wise on fsLR-5k"
  Do_cmd python "$MICAPIPE"/functions/build_mpc-vertex.py "$out" "$id" "$SES" "${mpc_p}" "mpc-swm" "$num_surfs"
  if [[ -f "${MPC_fsLR5k}" ]]; then ((Nsteps++)); fi
else Info "Subject ${id} has MPC vertex-wise on fsLR-5k"; ((Nsteps++)); ((N++)); fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status "MPC-SWM"
micapipe_procStatus "${id}" "${SES/ses-/}" "MPC-SWM-${mpc_str}" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "MPC-SWM-${mpc_str}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
