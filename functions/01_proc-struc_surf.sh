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
# Github: https://github.com/Deep-MI/FastSurfer
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
export threads=$6
tmpDir=$7
surfdir=$8
FastSurfer=$9
PROC=${10}
here=$(pwd)
export OMP_NUM_THREADS=$threads

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

# Check inputs: Nativepro T1
if [ ! -f "${T1nativepro}" ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_surf.json"
if [ -f "${module_json}" ] && [ $(grep "Status" "${module_json}" | awk -F '"' '{print $4}')=="COMPLETED" ]; then
Warning "Subject ${idBIDS} has been processed with -proc_surf
                If you want to re-run this step again, first erase all the outputs with:
                micapipe_cleanup -sub <subject_id> -out <derivatives> -bids <BIDS_dir> -proc_surf"; exit; fi

#------------------------------------------------------------------------------#
Title "Surface processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
# print the names on the terminal
bids_print.variables
bids_print.variables-structural
Note "Preprocessed surface directory: $surfdir"

# # Create script specific temp directory
tmp="${tmpDir}/${id}_micapipe_proc-surf_${RANDOM}"
Do_cmd mkdir -p "${tmp}/nii"
Note "Saving temporal dir: $nocleanup"
Note "\t\ttmp:" "${tmp}"

# GLOBAL variables for this script
Info "Proc_surf will use $threads threads"

#	Timer and steps progress
aloita=$(date +%s)
N=0
Nsteps=0
recon_log="${dir_subjsurf}/scripts/recon-all.log"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# IF SURFACE directory is provided create a symbolic link
if [[ "$surfdir" != "FALSE" ]]; then ((N++))
    if [[ -d "$surfdir" ]]; then
        Info "Copying from surface directory"
        Do_cmd mkdir "$dir_subjsurf"
        Do_cmd ln -s "$surfdir"/* "$dir_subjsurf"
    elif [[ ! -d "$surfdir" ]]; then
        Error "The provided surface directory does not exist: $surfdir"
        exit
    fi

# If not, prepare the BIDS T1w for freesurfer
elif [[ "$surfdir" == "FALSE" ]]; then ((N++))
    Info "Running Freesurfer"
    # Define SUBJECTS_DIR for freesurfer processing as a global variable
    # Will work on a temporal directory
    export SUBJECTS_DIR=${tmp}

    # Run freesurfer
    Do_cmd recon-all -cm -all -parallel -openmp ${threads} -expert "$EXPERT_FILE" -i "${T1_n4}" -s "$idBIDS"

    # Copy the freesurfer log to our MICA-log Directory
    Do_cmd cp "${tmp}/${idBIDS}/scripts/recon-all.log" "${dir_logs}/recon-all.log"

    # Copy results from TMP to deviratives/SUBJECTS_DIR directory
    Do_cmd cp -r "${tmp}/${idBIDS}" "$dir_surf"

    Info "Check log file:\n\t\t\t ${dir_logs}/recon-all.log"
fi

# -----------------------------------------------------------------------------------------------
# Check proc_surf status
if [[ -f "${dir_logs}/recon-all.log" ]] && grep -q "finished without error" "${dir_logs}/recon-all.log"; then ((Nsteps++)); fi

# Create json file for T1native
if [[ "$surfdir" == "FALSE" ]]; then
    proc_surf_json="${proc_struct}/${idBIDS}_proc_surf-${}.json"
    json_surf "$T1_n4" "$Nimgs" "${bids_T1ws[*]}" "${proc_surf_json}"
fi

# Notification of completition
micapipe_completition_status proc_surf
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_surf" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "proc_surf" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
