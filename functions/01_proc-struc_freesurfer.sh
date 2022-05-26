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
tmp="${tmpDir}/${RANDOM}_micapipe_proc-freesurfer_${id}"
Do_cmd mkdir -p "$tmp"

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
Note "\t\ttmp:" "${tmpDir}"

#	Timer
aloita=$(date +%s)

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM
Info "Preprocessed freesurfer directory: $FSdir"
if [[ "$FSdir" != "FALSE" ]]; then
    if [[ -d "$FSdir" ]]; then
        Info "Copying from freesurfer_dir"
        Do_cmd mkdir "$dir_freesurfer"
        Do_cmd cp -Rf "$FSdir"/* "$dir_freesurfer"
    elif [[ ! -d "$FSdir" ]]; then
        Error "The provided freesurfer directory does not exist: $FSdir"
        exit
    fi
elif [[ "$FSdir" == "FALSE" ]]; then
    Info "Running Freesurfer"
    # Define SUBJECTS_DIR for freesurfer processing as a global variable
    # Will work on a temporal directory
    export SUBJECTS_DIR=${tmp}
    if [ ! -d "$tmp/nii" ]; then Do_cmd mkdir "$tmp/nii"; fi

    # Perform recon-all surface registration
    if [[ "$hires" == "TRUE" ]]; then
        Info "Running recon with native submillimeter resolution"
        # export EXPERT_FILE=${tmp}/expert.opts
        # echo "mris_inflate -n 100" > "$EXPERT_FILE"
        # Optimize the contrast of the T1nativepro
        pve2=${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_pve_2.nii.gz
        T1_n4="${tmp}/${idBIDS}_space-nativepro_t1w_N4w.nii.gz"
        Do_cmd N4BiasFieldCorrection -r -d 3 -w ${pve2} -i "$T1nativepro" -o "$T1_n4"
        # Run freesurfer
        # Do_cmd recon-all -all -s "$idBIDS" -cm -i "$T1_n4" -expert "$EXPERT_FILE" -openmp ${threads}
        Do_cmd recon-all -cm -all -i "$T1_n4" -s "$idBIDS" -openmp ${threads}
        # Fix the inflation
        #Do_cmd mris_inflate -n 15 "${tmp}/${idBIDS}"/surf/?h.smoothwm "${tmp}/${idBIDS}"/surf/?h.inflated
    else
        # Copy all the T1 from the BIDS directory to the TMP
        # transform to NIFTI (if == NIFTI_GZ)
        for t1 in ${bids_T1ws[@]}; do
            t1_name=$(echo "$t1" | awk -F 'anat/' '{print $2}')
            if [[ "$t1" == *'.gz'* ]]; then
                  Do_cmd fslchfiletype NIFTI "$t1" "$tmp/nii/${t1_name/.gz/}"
            else  Do_cmd cp "$t1" "$tmp/nii/"
            fi
        done

        # List of Files for processing
        fs_cmd=$(echo "-i $(echo "$tmp"/nii/*nii | sed 's: : -i :g')")
        Do_cmd recon-all -cm -all "$fs_cmd" -s "$idBIDS" -openmp ${threads}
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
