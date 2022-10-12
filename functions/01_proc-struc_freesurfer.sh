#!/bin/bash
#
# T1w Structural processing with bash:
#
# Preprocessing workflow for structural T1w.
#
# This workflow makes use of FREESURFER, FSL (fslchfiletype)
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
# FastSurfer
# https://doi.org/10.1016/j.neuroimage.2020.117012
# For future implementation: https://github.com/Deep-MI/FastSurfer
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
export threads=$6
tmpDir=$7
FSdir=$8
hires=$9
t1wStr=${10}
UNI=${11}
MF=${12}
PROC=${13}
here=$(pwd)
export OMP_NUM_THREADS=$threads

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Manage manual inputs: T1w images
if [[ "$t1wStr" != "DEFAULT" ]]; then
  IFS=',' read -ra bids_t1wStr <<< "$t1wStr"
  for i in "${!bids_t1wStr[@]}"; do bids_t1wStr[i]=$(ls "${subject_bids}/anat/${idBIDS}_${bids_t1wStr[$i]}.nii"* 2>/dev/null); done
  bids_T1ws=("${bids_t1wStr[@]}")
fi

# End script if no T1 are found
Nimgs="${#bids_T1ws[*]}"  # total number of T1w
if [ "$Nimgs" -lt 1 ]; then Error "Subject $id doesn't have T1 on: \n\t\t\t${subject_bids}/anat"; exit; fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_surf.json"
micapipe_check_json_status "${module_json}" "proc_surf"

# If UNi is selected and multiple t1Str (3) are included the script will assing possitional values:
# 1:UNI, 2:INV1, 3:INV2
if [[ "${UNI}" == "TRUE" ]]; then
  if [ "$Nimgs" -gt 1 ]; then
    N4wm="TRUE"
    bids_inv1=${bids_T1ws[1]}
    bids_inv2=${bids_T1ws[2]}
    bids_T1ws=(${bids_T1ws[0]})
  fi
  # Look for Inversion-1 and Inversion-2 if UNI is selected
  if [ ! -f "${bids_inv1}" ]; then Error "Subject $id doesn't have INV1 on: \n\t${bids_inv1}"; exit; fi
  if [ ! -f "${bids_inv2}" ]; then Error "Subject $id doesn't have INV2 on: \n\t${bids_inv2}"; exit; fi
fi

# For hires check that proc_structural nativepro exists
if [[ "$hires" = "TRUE" ]] && [[ ! -f "${proc_struct}/${idBIDS}"_space-nativepro_t1w.nii.gz ]]; then
  Error "Submilimetric (hires) processing of fresurfer requires the T1_nativepro: RUN -proc_structural first"; exit
fi

#------------------------------------------------------------------------------#
Title "Structural processing: Freesurfer\n\t\tmicapipe $Version, $PROC "
micapipe_software
# print the names on the terminal
bids_print.variables
bids_print.variables-structural
Note "Preprocessed freesurfer directory: $FSdir"

# # Create script specific temp directory
tmp="${tmpDir}/${id}_micapipe_proc-freesurfer_${RANDOM}"
Do_cmd mkdir -p "${tmp}/nii"
Note "Saving temporal dir: $nocleanup"
Note "\t\ttmp:" "${tmp}"

# GLOBAL variables for this script
Info "Freesurfer will use $threads threads"

#	Timer and steps progress
aloita=$(date +%s)
N=0
Nsteps=0
recon_log="${dir_subjsurf}/scripts/recon-all.log"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# IF FREESURFER directory is provided create a symbolic link
if [[ "$FSdir" != "FALSE" ]]; then ((N++))
    if [[ -d "$FSdir" ]]; then
        Info "Copying from freesurfer_dir"
        Do_cmd mkdir "$dir_subjsurf"
        Do_cmd ln -s "$FSdir"/* "$dir_subjsurf"
    elif [[ ! -d "$FSdir" ]]; then
        Error "The provided freesurfer directory does not exist: $FSdir"
        exit
    fi

# If not, prepare the BIDS T1w for freesurfer
elif [[ "$FSdir" == "FALSE" ]]; then ((N++))
    Info "Running Freesurfer"
    # Define SUBJECTS_DIR for freesurfer processing as a global variable
    # Will work on a temporal directory
    export SUBJECTS_DIR=${tmp}
    T1fs="${tmp}/${idBIDS}_t1w_mean_n4.nii.gz"

    # Differential workflow if is mp2rage or mprage
    if  [[ "${UNI}" == "TRUE" ]]; then
      Info "Removing background noise from mp2rage UNI-T1map"
      UNIdns="${tmp}/${idBIDS}_mp2rage-uni_denoised.nii.gz"
      # UNI mp2rage denoising
      Do_cmd ${MICAPIPE}/functions/mp2rage_denoise.py "${bids_T1ws[0]}" "${bids_inv1[0]}" "${bids_inv2[0]}" "$UNIdns" --mf "${MF}"
      bids_T1ws=("$UNIdns")
    fi

    # transform to NIFTI_GZ (if == NIFTI)
    for t1 in ${bids_T1ws[@]}; do
        t1_name=$(echo "$t1" | awk -F 'anat/' '{print $2}')
        Do_cmd fslchfiletype NIFTI_GZ "$t1" "$tmp/nii/${t1_name/.gz/}.gz"
    done

    # List of files for freesurfer
    fs_T1ws=(${tmp}/nii/*)

    #  If multiple T1w were provided, Register and average to the first T1w
    if [ "${#fs_T1ws[@]}" -gt 1 ]; then
        ref=${fs_T1ws[0]} # reference to registration
        t1ref="T1wRef"
        # Loop over each T1
        Nimgs=${#fs_T1ws[@]} # total number of T1w
        n=$((Nimgs - 1))
        for ((i=1; i<=n; i++)); do
            run=$((i+1))
            T1mat_str="${tmp}/t1w_from-run-${run}_to_T1wRef_"
            T1mat="${T1mat_str}0GenericAffine.mat"
            T1run_2_T1ref="${tmp}/t1w_from-run-${run}_to_T1wRef.nii.gz"
            Do_cmd antsRegistrationSyN.sh -d 3 -m "${bids_T1ws[i]}" -f "$ref" -o "$T1mat_str" -t a -n "$threads" -p d
            Do_cmd antsApplyTransforms -d 3 -i "${bids_T1ws[i]}" -r "$ref" -t "$T1mat" -o "$T1run_2_T1ref" -u int
        done
        # Calculate the mean over all T1w registered to the 1st t1w run
        t1s_reg=$(find "${tmp}/t1w_from-run-"*"_to_T1wRef.nii.gz")
        t1s_add=$(echo "-add $(echo ${t1s_reg} | sed 's: : -add :g')")
        Do_cmd fslmaths "$ref" "$t1s_add" -div "$N" "$T1fs" -odt float
    # If only one T1w is provided
    elif [ "$N" -eq 1 ]; then
      mv ${fs_T1ws[0]} $T1fs
    fi

    # Perform freesurfer-all surface registration
    T1_n4="${tmp}/${idBIDS}_space-native_t1w_N4wm.nii.gz"
    if [[ "$hires" == "TRUE" ]]; then
        # Affine transformation between T1nativepro and native T1w for freesurfer
        mov="${proc_struct}/${idBIDS}"_space-nativepro_t1w.nii.gz
        T1mat_str="${tmp}/${idBIDS}_from-nativepro_to_native_"
        T1mat="${T1mat_str}0GenericAffine.mat"
        pve2="${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_pve_2.nii.gz"
        pve2_nat="${tmp}/space-native_t1w_brain_pve_2.nii.gz"
        Info "Registering T1w_run-${run} to ${t1ref}"
        Do_cmd antsRegistrationSyN.sh -d 3 -m "${mov}" -f "${T1fs}" -o "${T1mat_str}" -t a -n "${threads}" -p d
        Do_cmd antsApplyTransforms -d 3 -i "${pve2}" -r "${T1fs}" -t "${T1mat}" -o "${pve2_nat}" -u int -n GenericLabel

        # 7T - Optimize the contrast of the T1w with the White matter probabilistic mask
        N4wm="TRUE";  N4bfc="FALSE"
        Do_cmd N4BiasFieldCorrection -r 1 -d 3 -w ${pve2_nat} -i "${T1fs}" -o "${T1_n4}"

        Info "Running recon-all for 7T (hires)"
        export EXPERT_FILE=${tmp}/expert.opts
        echo "mris_inflate -n 100" > "$EXPERT_FILE"

        # Run freesurfer
        Do_cmd recon-all -cm -all -parallel -openmp ${threads} -expert "$EXPERT_FILE" -i "${T1_n4}" -s "$idBIDS"
        # Fix the inflation
        Do_cmd mris_inflate -n 100 "${tmp}/${idBIDS}"/surf/?h.smoothwm "${tmp}/${idBIDS}"/surf/?h.inflated
    else
        # Normalize intensities
        N4wm="FALSE"; N4bfc="TRUE"
        Do_cmd N4BiasFieldCorrection -r 1 -d 3 -i "${T1fs}" -o "${T1_n4}"

        # Run FREESURFER recon-all
        Do_cmd recon-all -cm -all -parallel -openmp ${threads} -i "${T1_n4}" -s "$idBIDS"
    fi

    # Copy the freesurfer log to our MICA-log Directory
    Do_cmd cp "${tmp}/${idBIDS}/scripts/recon-all.log" "${dir_logs}/recon-all.log"

    # Copy results to  freesurfer's SUBJECTS_DIR directory
    Do_cmd cp -r "${tmp}/${idBIDS}" "$dir_surf"

    Info "Check log file:\n\t\t\t ${dir_logs}/recon-all.log"
fi

# -----------------------------------------------------------------------------------------------
# Check proc_freesurfer status
if [[ -f "${dir_logs}/recon-all.log" ]] && grep -q "finished without error" "${dir_logs}/recon-all.log"; then ((Nsteps++)); fi

# Create json file for T1native
if [[ "$FSdir" == "FALSE" ]]; then
    freesurfer_json="${proc_struct}/${idBIDS}_proc_freesurfer.json"
    json_freesurfer "$T1_n4" "$Nimgs" "${bids_T1ws[*]}" "${freesurfer_json}"
fi

# Notification of completition
micapipe_completition_status proc_freesurfer
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_freesurfer" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "proc_freesurfer" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
