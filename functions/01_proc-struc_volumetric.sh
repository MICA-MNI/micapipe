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
umask 003
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
t1wStr=$8
N4wm=$9
PROC=${10}
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

# Manage manual inputs: T1w images
if [[ "$t1wStr" != "DEFAULT" ]]; then
  IFS=',' read -ra bids_t1wStr <<< "$t1wStr"
  for i in "${!bids_t1wStr[@]}"; do bids_t1wStr[i]=$(ls "${subject_bids}/anat/${idBIDS}_${bids_t1wStr[$i]}.nii"* 2>/dev/null); done
  bids_T1ws=("${bids_t1wStr[@]}")
  N="${#bids_t1wStr[*]}"
  Info "Manually selected T1w string(s): $t1wStr, N=${N}"
fi

#------------------------------------------------------------------------------#
Title "Structural processing\n\t\tmicapipe $Version, $PROC"
micapipe_software
# print the names on the terminal
bids_print.variables
Info "Saving temporal dir: $nocleanup"

# GLOBAL variables for this script
Info "ANTs will use $threads threads"

#	Timer
aloita=$(date +%s)
Nsteps=0
umask 011

# Create script specific temp directory
tmp=${tmpDir}/${RANDOM}_micapipe_proc_struc-vol_${id}
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# BIDS T1w processing
N=${#bids_T1ws[@]} # total number of T1w
n=$((N - 1))
T1str_nat="${idBIDS}_space-nativepro_t1w"
T1n4="${tmp}/${T1str_nat}_n4.nii.gz"
T1nativepro="${proc_struct}/${T1str_nat}.nii.gz"
T1nativepro_brain="${T1nativepro/.nii.gz/_brain.nii.gz}"
T1nativepro_first="${proc_struct}/first/${T1str_nat}.nii.gz"
T1nativepro_5tt="${T1nativepro/.nii.gz/_5TT.nii.gz}"
T1nativepro_mask="${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_mask.nii.gz"

# Creates the t1w_nativepro for structural processing
if [ ! -f "${proc_struct}/${T1str_nat}".nii.gz ] || [ ! -f "${proc_struct}/${T1str_nat}"_brain.nii.gz ]; then
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
    if [ "$N" -gt 1 ]; then
      reo_T1ws=(${tmp}/reo_*)
      ref=${reo_T1ws[0]} # reference to registration
      ref_run=$(echo "${reo_T1ws[0]}" | awk -F 'run-' '{print $2}'| sed 's:_T1w.nii.gz::g')
      t1ref="run-${ref_run}"
      # Loop over each T1
      for ((i=1; i<=n; i++)); do
          run=$(echo "${reo_T1ws[i]}" | awk -F 'run-' '{print $2}'| sed 's:_T1w.nii.gz::g')
          T1mat_str="${dir_warp}/${idBIDS}_t1w_from-run-${run}_to_${t1ref}_"
          T1mat="${T1mat_str}0GenericAffine.mat"
          T1run_2_T1="${tmp}/${id}_t1w_from-run-${run}_to_${t1ref}.nii.gz"
          Info "Registering T1w_run-${run} to ${t1ref}"
          Do_cmd antsRegistrationSyN.sh -d 3 -m "${reo_T1ws[i]}" -f "$ref"  -o "$T1mat_str" -t a -n "$threads" -p d
          Do_cmd antsApplyTransforms -d 3 -i "${reo_T1ws[i]}" -r "$ref" -t "$T1mat" -o "$T1run_2_T1" -v -u int
      done
      # Calculate the mean over all T1w registered to the 1st t1w run
      t1s_reg=$(find "${tmp}"/*_to_"${t1ref}".nii.gz)
      t1s_add=$(echo "-add $(echo ${t1s_reg} | sed 's: : -add :g')")
      Do_cmd fslmaths "$ref" "$t1s_add" -div "$N" "$T1n4" -odt float

    # If only one T1w is provided
    elif [ "$N" -eq 1 ]; then
      Do_cmd fslchfiletype NIFTI_GZ "$t1_reo"
      Do_cmd cp -v "$t1_reo"* "$T1n4"
    fi

    Info "T1w_nativepro biasfield correction and intensity rescaling"

    # Intensity Non-uniform correction - N4
    Do_cmd N4BiasFieldCorrection  -d 3 -i "$T1n4" -r -o "$T1n4" -v

    # Rescale intensity [100,0]
    Do_cmd ImageMath 3 "$T1nativepro" RescaleImage "$T1n4" 0 100

    # If no T1native pro exit_status "something is wrong" exit
    if [ ! -f "$T1nativepro" ]; then Error "$T1str_nat was not generated"; Do_cmd exit; fi

    # Brainmask
    Do_cmd bet "$T1nativepro" "$T1nativepro_brain" -B -f 0.25 -v

    # If no T1native pro exit_status "something is wrong" exit
    if [ ! -f "$T1nativepro_brain" ]; then Error "$T1str_nat masked was not generated"; Do_cmd exit; else ((Nsteps++)); fi

    # Create json file for T1native
    json_nativepro_t1w "$T1nativepro" "$N" "${bids_T1ws[*]}" "${proc_struct}/${T1str_nat}.json"
else
    Info "Subject ${id} has a t1w_nativepro and t1w_nativepro_brain"; ((Nsteps++))
fi

# FSL first on the t1w_nativepro
unset SGE_ROOT
export FSLPARALLEL=0
firstout=${T1nativepro_first/.nii.gz/_all_fast_firstseg.nii.gz}
if [ ! -f "$firstout" ]; then
    Info "FSL first is running, output file: ${T1str_nat}\n\t\t\t ${firstout}"
    Do_cmd run_first_all -i "$T1nativepro_brain" -o "$T1nativepro_first" -b &
    wait $!
    until [ -f "$firstout" ]; do sleep 120; done
    if [ -f "$firstout" ]; then ((Nsteps++)); fi
  else
    Info "Subject $id has FSL-first"; ((Nsteps++))
fi

# FSL FAST on the t1w_nativepro
if [ ! -f "${T1nativepro_brain/.nii.gz/_pve_0.nii.gz}" ]; then
    Do_cmd fast -N "$T1nativepro_brain"
    if [ -f "${T1nativepro_brain/.nii.gz/_pve_0.nii.gz}" ]; then ((Nsteps++)); fi
else
    Info "Subject $id has FSL-fast"; ((Nsteps++))
fi

if [[ "$N4wm" == "TRUE" ]]; then
  Info "N4 bias field corretion weighted by white matter"
  pve2=${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_pve_2.nii.gz
  T1_n4="${tmp}/${idBIDS}_space-nativepro_t1w_N4w.nii.gz"
  Do_cmd N4BiasFieldCorrection -r -d 3 -w ${pve2} -i "$T1nativepro" -o "$T1_n4"
  Do_cmd ImageMath 3 "$T1nativepro" RescaleImage "$T1_n4" 0 100
fi

# Loop over all requested templates.
# mmTemplates is a fixed value=(0.8 2) of the MNI152 atlas resolution
for mm in 2 0.8; do
  # Only runs if the output doesn't exist
  if [ ! -f "${dir_warp}/${idBIDS}_from-nativepro_brain_to-MNI152_${mm}mm_mode-image_desc-SyN_1Warp.nii.gz" ]; then
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
      Info "Subject $id has nativepro_t1w on MNI152 space"; ((Nsteps++))
  fi
done
Do_cmd rm -rf ${dir_warp}/*Warped.nii.gz 2>/dev/null

# Update the T1native mask and T1native_brain
T1nativepro_maskjson="${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_mask.json"
if [ ! -f "$T1nativepro_maskjson" ]; then
    Do_cmd antsApplyTransforms -d 3 -n GenericLabel -i "$MNI152_mask" -r "$T1nativepro_brain" \
            -t ["$T1_MNI152_affine",1] -t "$T1_MNI152_InvWarp" -o "$T1nativepro_mask" -v
    Do_cmd ImageMath 3 "$T1nativepro_brain" m "$T1nativepro" "$T1nativepro_mask"
    json_nativepro_mask "$T1nativepro_brain" "$MNI152_mask" "$T1nativepro_maskjson"\
        "antsApplyTransforms -d 3 -n GenericLabel -i ${MNI152_mask} -r ${T1nativepro_brain} -t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp} -o ${T1nativepro_mask}"
    ((Nsteps++))
else
        Info "Subject $id has a masked t1w_${mm}mm_nativepro"; ((Nsteps++))
fi

# Generate a five-tissue-type image for anatomically constrained tractography
if [[ ! -f "$T1nativepro_5tt" ]]; then
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
    if [[ -f "$T1nativepro_5tt" ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has 5TT nifti"; ((Nsteps++))
fi

# -----------------------------------------------------------------------------------------------
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
N=7
if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "Volumetric structural processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:\n\tlogs:
\tSteps completed : $(printf "%02d" "$Nsteps")/07
\tStatus          : ${status}
\tCheck logs      : $(ls "$dir_logs"/proc_structural_*.txt)"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_structural" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_structural" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
cleanup "$tmp" "$nocleanup" "$here"
