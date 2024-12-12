#!/bin/bash
#
# T1w Structural processing with bash:
#
# Preprocessing workflow for structural T1w.
#
# This workflow makes use of ANTS, FSL and AFNI
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
T1wStr=$8
UNI=${9}
MF=${10}
PROC=${11}
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

# Manage manual inputs: T1w images
if [[ "$T1wStr" != "DEFAULT" ]]; then
  IFS=',' read -ra bids_T1wStr <<< "$T1wStr"
  for i in "${!bids_T1wStr[@]}"; do bids_T1wStr[i]=$(ls "${subject_bids}/anat/${idBIDS}_${bids_T1wStr[$i]}.nii"* 2>/dev/null); done
  bids_T1ws=("${bids_T1wStr[@]}")
fi

# End script if no T1 are found
Nimgs="${#bids_T1ws[*]}"  # total number of T1w
if [ "$Nimgs" -lt 1 ]; then Error "Subject $id doesn't have T1 on: \n\t\t\t${subject_bids}/anat"; exit; fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_structural.json"
micapipe_check_json_status "${module_json}" "proc_structural"

# If UNi is selected and multiple t1Str (3) are included the script will assing possitional values:
# 1:UNI, 2:INV1, 3:INV2
if [[ "${UNI}" == "TRUE" ]]; then
  if [ "$Nimgs" -gt 1 ]; then
    bids_inv1=${bids_T1ws[1]}
    bids_inv2=${bids_T1ws[2]}
    bids_T1ws=(${bids_T1ws[0]})
  fi
  # Look for Inversion-1 and Inversion-2 if UNI is selected
  if [ ! -f "${bids_inv1}" ]; then Error "Subject $id doesn't have INV1 on: \n\t${bids_inv1}"; exit; fi
  if [ ! -f "${bids_inv2}" ]; then Error "Subject $id doesn't have INV2 on: \n\t${bids_inv2}"; exit; fi
fi

#------------------------------------------------------------------------------#
Title "Structural processing\n\t\tmicapipe $Version, $PROC"
micapipe_software
# print the names on the terminal
bids_print.variables
bids_print.variables-structural

# GLOBAL variables for this script
Info "ANTs will use $threads threads"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_proc_struc-vol_${id}"
Do_cmd mkdir -p "$tmp" "${proc_struct}"/first
Note "Saving temporal dir:" "$nocleanup"
Note "\t\ttmp:" "${tmp}"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# BIDS T1w processing
Nimgs=${#bids_T1ws[@]} # total number of T1w
n=$((Nimgs - 1))
T1str_nat="${idBIDS}_space-nativepro_T1w"
T1reo="${tmp}/${T1str_nat}_reo.nii.gz"
T1n4="${tmp}/${T1str_nat}_n4.nii.gz"
T1nativepro="${proc_struct}/${T1str_nat}.nii.gz"
T1nativepro_brain="${proc_struct}/${idBIDS}_space-nativepro_T1w_brain.nii.gz"
T1nativepro_mask="${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_mask.nii.gz"
T1nativepro_first="${proc_struct}/first/${T1str_nat}.nii.gz"
T1nativepro_5tt="${T1nativepro/.nii.gz/_5tt.nii.gz}"
procstruct_json="${proc_struct}/${T1str_nat}.json"

# Creates the T1w_nativepro for structural processing
if [ ! -f "${proc_struct}/${T1str_nat}".nii.gz ] || [ ! -f "${proc_struct}/${T1str_nat}"_brain.nii.gz ]; then ((N++))
    # Reorient  to LPI with AFNI
    # LPI is the standard 'neuroscience' orientation, where the x-axis is
    # Left-to-Right, the y-axis is Posterior-to-Anterior, and the z-axis is Inferior-to-Superior.
    for ((k=0; k<=n; k++)); do
      nom=$(echo "${bids_T1ws[k]}" | awk -F '/' '{print $NF}')
      t1_reo=${tmp}/${nom/sub/reo}
      t1_obl=${tmp}/${nom/sub/obl}
      # Deoblique
      Do_cmd 3dresample -orient LPI -prefix "$t1_obl" -inset "${bids_T1ws[k]}"
      # Reorient to standart
      Do_cmd fslreorient2std "$t1_obl" "$t1_reo"
    done

    # If multiple T1w were provided, Register and average to the first T1w
    if [ "$Nimgs" -gt 1 ]; then
      reo_T1ws=(${tmp}/reo-*)
      ref=${reo_T1ws[0]} # reference to registration
      ref_run=$(echo "${reo_T1ws[0]}" | awk -F 'run-' '{print $2}'| sed 's:.nii.gz::g')
      t1ref="run-${ref_run}"
      # Loop over each T1
      for ((i=1; i<=n; i++)); do
          run=$(echo "${reo_T1ws[i]}" | awk -F 'run-' '{print $2}'| sed 's:.nii.gz::g')
          T1mat_str="${tmp}/${idBIDS}_T1w_from-run-${run}_to_${t1ref}_"
          T1mat="${T1mat_str}0GenericAffine.mat"
          T1run_2_T1="${tmp}/${id}_T1w_from-run-${run}_to_${t1ref}.nii.gz"
          Info "Registering T1w_run-${run} to ${t1ref}"
          Do_cmd antsRegistrationSyN.sh -d 3 -m "${reo_T1ws[i]}" -f "$ref"  -o "$T1mat_str" -t a -n "$threads" -p d
          Do_cmd antsApplyTransforms -d 3 -i "${reo_T1ws[i]}" -r "$ref" -t "$T1mat" -o "$T1run_2_T1" -v -u int
      done
      # Calculate the mean over all T1w registered to the 1st T1w run
      t1s_reg=$(find "${tmp}"/*_to_"${t1ref}".nii.gz)
      t1s_add=$(echo "-add $(echo ${t1s_reg} | sed 's: : -add :g')")
      Do_cmd fslmaths "$ref" "$t1s_add" -div "$Nimgs" "$T1reo" -odt float

    # If only one T1w is provided
    elif [ "$Nimgs" -eq 1 ]; then
      Do_cmd fslchfiletype NIFTI_GZ "$t1_reo"
      Do_cmd cp -v "$t1_reo"* "$T1reo"
    fi

    # Differential workflow if is mp2rage or mprage
    if  [[ "${UNI}" == "TRUE" ]]; then
      Info "Removing background noise from mp2rage UNI-T1map"
      # UNI mp2rage denoising
      Do_cmd "${MICAPIPE}"/functions/mp2rage_denoise.py "$T1reo" "${bids_inv1[0]}" "${bids_inv2[0]}" "$T1n4" --mf "${MF}"
    else
      mv "$T1reo" "$T1n4"
    fi

    Info "T1w_nativepro biasfield correction and intensity rescaling"
    # Intensity Non-uniform correction - N4
    Do_cmd N4BiasFieldCorrection  -d 3 -i "$T1n4" -r -o "$T1n4" -v

    # Rescale intensity [100,0]
    Do_cmd ImageMath 3 "$T1nativepro" RescaleImage "$T1n4" 0 100

    # If no T1native pro exit_status "something is wrong" exit
    if [ ! -f "$T1nativepro" ]; then Error "$T1str_nat was not generated"; Do_cmd exit; fi

    # Brainmask
    Do_cmd mri_synthstrip -i "$T1nativepro" -o "$T1nativepro_brain" -m "$T1nativepro_mask" --no-csf

    # If no T1native pro exit_status "something is wrong" exit
    if [ ! -f "$T1nativepro_brain" ]; then Error "$T1str_nat masked was not generated"; Do_cmd exit; else ((Nsteps++)); fi

    # Create json file for T1native
    json_nativepro_T1w "$T1nativepro" "$Nimgs" "${bids_T1ws[*]}" "${procstruct_json}"
else
    Info "Subject ${id} has a T1w_nativepro and T1w_nativepro_brain"; ((Nsteps++)); ((N++))
fi

# Loop over all requested templates.
# mmTemplates is a fixed value=(0.8 2) of the MNI152 atlas resolution
for mm in 2 0.8; do
  # Only runs if the output doesn't exist
  if [ ! -f "${dir_warp}/${idBIDS}_from-nativepro_brain_to-MNI152_${mm}mm_mode-image_desc-SyN_1Warp.nii.gz" ]; then
      ((N++))
      Info "Registration of T1w nativepro to MNI152 ${mm}mm"
      # ---------------------------------------------------------
      # MNI152 brain templates
      MNI152_brain="${util_MNIvolumes}/MNI152_T1_${mm}mm_brain.nii.gz"
      # T1 to MNI152 transformation matrix and warp field names
      str_MNI152_SyN="${dir_warp}/${idBIDS}_from-nativepro_brain_to-MNI152_${mm}mm_mode-image_desc-SyN_"
      # ANTs - 3 steps registration: rigid, affine and SyN
      Do_cmd antsRegistrationSyN.sh -d 3 -f "$MNI152_brain" -m "$T1nativepro_brain" -o "$str_MNI152_SyN" -t s -n "$threads"
      ((Nsteps++))
  else
      Info "Subject $id has nativepro_T1w on MNI152 space"; ((Nsteps++)); ((N++))
  fi
done

# Variables
mat_MNI152_SyN="${dir_warp}/${idBIDS}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_"    # transformation strings nativepro to MNI152_0.8mm
T1_MNI152_InvWarp="${mat_MNI152_SyN}1InverseWarp.nii.gz"                      # Inversewarp - nativepro to MNI152_0.8mm
T1_MNI152_Warp="${mat_MNI152_SyN}1Warp.nii.gz"                      # Inversewarp - nativepro to MNI152_0.8mm
T1_MNI152_affine="${mat_MNI152_SyN}0GenericAffine.mat"                        # Affine matrix - nativepro to MNI152_0.8mm
MNI152_2mm="${util_MNIvolumes}/MNI152_T1_2mm_brain.nii.gz"
MNI152_08mm="${util_MNIvolumes}/MNI152_T1_0.8mm_brain.nii.gz"
MNI151_to_T1nativepro="-t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp}"
T1nativepro_to_MNI152="-t ${T1_MNI152_Warp} -t ${T1_MNI152_affine}"

# Save transformations in json file
proc_struct_transformations "$T1nativepro_brain" "${MNI152_08mm}" "${MNI152_2mm}" "${dir_warp}/${idBIDS}_transformations-proc_structural.json" '${T1nativepro_to_MNI152}' '${MNI151_to_T1nativepro}'

# Update the T1native mask and T1native_brain
T1nativepro_maskjson="${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_mask.json"
if [ ! -f "$T1nativepro_maskjson" ]; then ((N++))
    Do_cmd antsApplyTransforms -d 3 -n GenericLabel -i "$MNI152_mask" -r "$T1nativepro_brain" \
            "${MNI151_to_T1nativepro}" -o "$T1nativepro_mask" -v
    Do_cmd ImageMath 3 "$T1nativepro_brain" m "$T1nativepro" "$T1nativepro_mask"
    json_nativepro_mask "$T1nativepro_brain" "$MNI152_mask" "$T1nativepro_maskjson"\
        "antsApplyTransforms -d 3 -n GenericLabel -i ${MNI152_mask} -r ${T1nativepro_brain} -t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp} -o ${T1nativepro_mask}"
    ((Nsteps++))
else
        Info "Subject $id has a masked T1w_${mm}mm_nativepro"; ((Nsteps++)); ((N++))
fi

# FSL first on the T1w_nativepro
unset SGE_ROOT
export FSLPARALLEL=0
T1str_atlas="${idBIDS}_space-nativepro_T1w_atlas"
T1_seg_subcortex="${dir_volum}/${T1str_atlas}-subcortical.nii.gz"
firstout=${T1nativepro_first/.nii.gz/_all_fast_firstseg.nii.gz}

if [ ! -f "$T1_seg_subcortex" ] || [ ! -f "$T1nativepro_5tt" ]; then ((N++))
    Info "FSL first is running"
    Note "output file:" "${firstout}/${T1str_nat}"
    Do_cmd run_first_all -i "$T1nativepro_brain" -o "$T1nativepro_first" -b &
    wait $!
    until [ -f "$firstout" ]; do sleep 120; done
    if [ -f "$firstout" ]; then Do_cmd cp "$firstout" "$T1_seg_subcortex"; ((Nsteps++)); fi
else
    Info "Subject $id has a Subcortical parcellation"; ((Nsteps++)); ((N++))
fi

# FSL FAST on the T1w_nativepro
if [ ! -f "${T1nativepro_brain/.nii.gz/_pve_0.nii.gz}" ]; then ((N++))
    Do_cmd fast -N "$T1nativepro_brain"
    if [ -f "${T1nativepro_brain/.nii.gz/_pve_0.nii.gz}" ]; then ((Nsteps++)); fi
else
    Info "Subject $id has FSL-fast"; ((Nsteps++)); ((N++))
fi

# Bias field correction weighted by white matter (e.g. 7T data or lost of signal in temporal areas
N4wmStatus_check=$(grep "N4wmProcessed" "${procstruct_json}" | awk -F '"' '{print $4}')
if [ "$N4wmStatus_check" == "FALSE" ]; then ((N++))
  Info "N4 bias field corretion weighted by white matter"
  pve2=${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_pve_2.nii.gz
  T1_n4="${tmp}/${idBIDS}_space-nativepro_T1w_N4w.nii.gz"
  Do_cmd N4BiasFieldCorrection -r -d 3 -w "$pve2" -i "$T1nativepro" -o "$T1_n4"
  Do_cmd ImageMath 3 "$T1nativepro" RescaleImage "$T1_n4" 0 100
  export N4wmStatus="TRUE"
  # Update json file for T1native
  json_nativepro_T1w "$T1nativepro" "$Nimgs" "${bids_T1ws[*]}" "${procstruct_json}"; ((Nsteps++))
else
    Info "Subject ${id} has N4 bias field corretion weighted by white matter"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# Cerebellar parcellation from MNI152_0.8mm to T1-Nativepro
Info "Cerebellum parcellation to T1-nativepro Volume"
atlas_cerebellum="${util_MNIvolumes}/MNI152_T1_0.8mm_cerebellum.nii.gz"       # cerebellar lobules atlas
T1_seg_cerebellum="${dir_volum}/${T1str_atlas}-cerebellum.nii.gz"             # Cerebellar output

# Apply inverse transfrom from MNI152-cerebellum to T1-nativepro
if [[ ! -f "$T1_seg_cerebellum" ]]; then ((N++))
    Do_cmd antsApplyTransforms -d 3 -i "$atlas_cerebellum" \
                -r "$T1nativepro" \
                -n GenericLabel "${MNI151_to_T1nativepro}" \
                -o "$T1_seg_cerebellum" -v -u int
    if [[ -f "$T1_seg_cerebellum" ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a Cerebellum parcellation on T1-nativepro"; ((Nsteps++)); ((N++))
fi

# Generate a five-tissue-type image for anatomically constrained tractography
if [[ ! -f "$T1nativepro_5tt" ]]; then ((N++))
    # --------------------------------------------------------------
    # Step by step 5tt
    # Process data in tmp directory
    cd "$tmp"
    # Convert to mif format
    Do_cmd mrconvert "$T1nativepro_brain" "${tmp}/input.mif"
    # Change strides
    first_input='T1.nii'
    Do_cmd mrconvert input.mif "$first_input" -strides -1,+2,+3
    Info "Convert FIRST meshes to partial volume images:"
    firstr="${proc_struct}/first/${T1str_nat}-"
    sgm_structures=(L_Accu R_Accu L_Caud R_Caud L_Pall R_Pall L_Puta R_Puta L_Thal R_Thal)
    echo "         	${sgm_structures[@]}"
    for struct in "${sgm_structures[@]}"; do
      Info "Mesh to volume: $struct"
      pve_image_path="mesh2voxel_${struct}.mif"
      vtk_in_path="${firstr}${struct}_first.vtk"
      vtk_temp_path="${struct}.vtk"
      Do_cmd meshconvert "$vtk_in_path" "$vtk_temp_path" -transform first2real "$first_input"
      Do_cmd mesh2voxel "$vtk_temp_path" "$T1nativepro_brain" "$pve_image_path"
    done
    mrmath mesh2voxel_* sum - | mrcalc - 1.0 -min all_sgms.mif

    Info "Generating partial volume images for SGM structures"
    T1_pve_1="${T1nativepro/.nii.gz/_brain_pve_1.nii.gz}"
    T1_pve_2="${T1nativepro/.nii.gz/_brain_pve_2.nii.gz}"
    T1_pve_0="${T1nativepro/.nii.gz/_brain_pve_0.nii.gz}"

    mrthreshold "$T1_pve_2" - -abs 0.001 | maskfilter - connect - -connectivity | mrcalc 1 - 1 -gt -sub remove_unconnected_wm_mask.mif -datatype bit
    Do_cmd mrcalc "$T1_pve_0" remove_unconnected_wm_mask.mif -mult csf.mif
    Do_cmd mrcalc 1.0 csf.mif -sub all_sgms.mif -min sgm.mif
    Do_cmd mrcalc 1.0 csf.mif sgm.mif -add -sub "$T1_pve_1" "$T1_pve_2" -add -div multiplier.mif
    Do_cmd mrcalc multiplier.mif -finite multiplier.mif 0.0 -if multiplier_noNAN.mif
    Do_cmd mrcalc "$T1_pve_1" multiplier_noNAN.mif -mult remove_unconnected_wm_mask.mif -mult cgm.mif
    Do_cmd mrcalc "$T1_pve_2" multiplier_noNAN.mif -mult remove_unconnected_wm_mask.mif -mult wm.mif
    Do_cmd mrcalc 0 wm.mif -min path.mif
    mrcat cgm.mif sgm.mif wm.mif csf.mif path.mif - -axis 3 | mrconvert - combined_precrop.mif -strides +2,+3,+4,+1
    Do_cmd mv combined_precrop.mif result.mif
    Do_cmd mrconvert result.mif "$T1nativepro_5tt"
    Do_cmd cd "$here"
    if [[ -f "$T1nativepro_5tt" ]]; then ((Nsteps++)); slim_proc_struct; fi
else
    Info "Subject ${id} has 5TT nifti"; ((Nsteps++)); ((N++))
fi

# -----------------------------------------------------------------------------------------------
# Notification of completition
micapipe_completition_status proc_structural
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_structural" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus_json "${id}" "${SES/ses-/}" "proc_structural" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"