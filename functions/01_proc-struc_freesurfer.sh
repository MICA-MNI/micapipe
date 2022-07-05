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
umask 003
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
PROC=${11}
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
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
  N="${#bids_t1wStr[*]}"
  Info "Manually selected T1w string(s): $t1wStr, N=${N}"
fi

# BIDS T1w processing
N=${#bids_T1ws[@]} # total number of T1w
# End script if no T1 are found
if [ "$N" -lt 1 ]; then Error "Subject $id doesn't have T1 on: \n\t\t\t${subject_bids}/anat"; exit; fi

# # Create script specific temp directory
tmp="${tmpDir}/${id}_micapipe_proc-freesurfer_${RANDOM}"
Do_cmd mkdir -p "${tmp}/nii"

# Stop if freesurfer has finished without errors
if grep -q "finished without error" "${dir_freesurfer}/scripts/recon-all.log"; then
status="COMPLETED"; Nsteps=01
grep -v "${id}, ${SES/ses-/}, proc_freesurfer" "${out}/micapipe_processed_sub.csv" > "${tmp}/tmpfile" && mv "${tmp}/tmpfile" "${out}/micapipe_processed_sub.csv"
echo "${id}, ${SES/ses-/}, proc_freesurfer, ${status}, ${N}/01, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipe_processed_sub.csv"
Warning "Subject ${id} has Freesurfer
                    > If you want to re-run for QC purposes try it manually
                    > If you want to run again this step first erase all the outputs with:
                      mica_cleanup -sub <subject_id> -out <derivatives> -bids <BIDS_dir> -proc_fresurfer";
exit
fi
if [[ "$hires" = "TRUE" ]] && [[ ! -f "${proc_struct}/${idBIDS}"_space-nativepro_t1w.nii.gz ]]; then
  Error "Submilimetric (hires) processing of fresurfer requires the T1_nativepro: RUN -proc_structural first"; rm -rf "${tmp}"; exit
fi

#------------------------------------------------------------------------------#
Title "Structural processing: Freesurfer\n\t\tmicapipe $Version, $PROC "
micapipe_software
# print the names on the terminal
bids_print.variables
Info "Saving temporal dir: $nocleanup"
Note "\t\ttmp:" "${tmp}"

#	Timer
aloita=$(date +%s)

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM
Info "Preprocessed freesurfer directory: $FSdir"

# IF FREESURFER directory is provided create a symbolic link
if [[ "$FSdir" != "FALSE" ]]; then
    if [[ -d "$FSdir" ]]; then
        Info "Copying from freesurfer_dir"
        Do_cmd mkdir "$dir_freesurfer"
        Do_cmd ln -s "$FSdir"/* "$dir_freesurfer"
    elif [[ ! -d "$FSdir" ]]; then
        Error "The provided freesurfer directory does not exist: $FSdir"
        exit
    fi

# If not, prepare the BIDS T1w for freesurfer
elif [[ "$FSdir" == "FALSE" ]]; then
    Info "Running Freesurfer"
    # Define SUBJECTS_DIR for freesurfer processing as a global variable
    # Will work on a temporal directory
    export SUBJECTS_DIR=${tmp}
    T1fs="${tmp}/${idBIDS}_t1w_mean_n4.nii.gz"

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
        N=${#fs_T1ws[@]} # total number of T1w
        n=$((N - 1))
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

    # Perform recon-all surface registration
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
        T1_n4="${tmp}/${idBIDS}_space-nativepro_t1w_N4wm.nii.gz"
        Do_cmd N4BiasFieldCorrection -r 1 -d 3 -w ${pve2_nat} -i "${T1fs}" -o "${T1_n4}"

        Info "Running recon with native submillimeter resolution"
        export EXPERT_FILE=${tmp}/expert.opts
        echo "mris_inflate -n 100" > "$EXPERT_FILE"
        # Run freesurfer
        Do_cmd recon-all -cm -all -i "${T1_n4}" -s "$idBIDS" -expert "$EXPERT_FILE"
        # Fix the inflation
        Do_cmd mris_inflate -n 100 "${tmp}/${idBIDS}"/surf/?h.smoothwm "${tmp}/${idBIDS}"/surf/?h.inflated
    else
        # Run FREESURFER recon-all
        Do_cmd recon-all -cm -all -i "${T1fs}" -s "$idBIDS"
    fi

    # Copy the recon-all log to our MICA-log Directory
    Do_cmd cp -v "${tmp}/${idBIDS}/scripts/recon-all.log" "${dir_logs}/recon-all.log"

    # Copy results to  freesurfer's SUBJECTS_DIR directory
    Do_cmd cp -rv "${tmp}/${idBIDS}" "$dir_surf"

    Info "Check log file:\n\t\t\t ${dir_logs}/recon-all.log"
fi

# -----------------------------------------------------------------------------------------------
# Notification of completition
N=1 # total number of steps
if grep -q "finished without error" "${dir_freesurfer}/scripts/recon-all.log"; then status="COMPLETED"; Nsteps=01; else status="INCOMPLETE"; Nsteps=00; fi

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

Title "Freesurfer recon-all processing ended:
\tStatus          : ${status}
\tCheck logs      : $(ls "$dir_logs"/proc_freesurfer*.txt)"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_freesurfer" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_freesurfer" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
cleanup "$tmp" "$nocleanup" "$here"
