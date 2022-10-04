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
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v1.0.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check inputs: Nativepro T1
if [ ! -f "${T1nativepro}" ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi
if [[ "$FastSurfer" == "TRUE" ]]; then recon="fastsurfer"; else recon="freesurfer"; fi

# Freesurfer Directory
set_surface_directory "${recon}"
# Surface Directories
if [ ! -d "${dir_surf}" ]; then mkdir "${dir_surf}" && chmod -R 770 "${dir_surf}"; fi
if [ ! -L "${dir_surf}/fsaverage5" ]; then Do_cmd ln -s "$FREESURFER_HOME/subjects/fsaverage5/" "${dir_surf}"; fi
if [ ! -d "${dir_surf}/conte69" ]; then Do_cmd mkdir -p "${dir_surf}"/conte69/surf; cp ${MICAPIPE}/surfaces/*conte69.gii "${dir_surf}"/conte69/surf; fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_surf-${recon}.json"
if [ -f "${module_json}" ]; then
  status=$(grep "Status" "${module_json}" | awk -F '"' '{print $4}')
  if [ "$status" == "COMPLETED" ]; then
  Note "Proc_surf json" "${module_json}"
  Warning "Subject ${idBIDS} has been processed with -proc_surf
            If you want to re-run this step again, first erase all the outputs with:
            micapipe_cleanup -sub <subject_id> -out <derivatives> -bids <BIDS_dir> -proc_surf"; exit
  else
      Info "proc_surf is INCOMPLETE, processing will continute"
  fi
fi

#------------------------------------------------------------------------------#
Title "Surface processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
# print the names on the terminal
bids_print.variables
Note "Preprocessed surface directory: $surfdir"

# # Create script specific temp directory
tmp="${tmpDir}/${id}_micapipe_proc-surf_${RANDOM}"
Do_cmd mkdir -p "${tmp}/nii"
Note "Saving temporal dir:" "$nocleanup"
Note "Temporal dir:" "${tmp}"

# GLOBAL variables for this script
Info "Proc_surf will use $threads threads"

#	Timer and steps progress
aloita=$(date +%s)
N=0
Nsteps=0

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# IF SURFACE directory is provided create a symbolic link
if [[ "$surfdir" != "FALSE" ]]; then ((N++))
    if [[ -d "$surfdir" ]]; then
        Info "Copying from surface directory"
        Do_cmd mkdir -p "$dir_subjsurf"
        Do_cmd ln -s "$surfdir"/* "$dir_subjsurf"
    elif [[ ! -d "$surfdir" ]]; then
        Error "The provided surface directory does not exist: $surfdir"
        exit
    fi

# If not, get ready to run the surface reconstruccion
elif [[ "$surfdir" == "FALSE" ]]; then ((N++))

    # FIX FOV greater than 256
    dim=($(mrinfo ${T1nativepro} -size))
    res=($(mrinfo ${T1nativepro} -spacing))
    fov=$(printf "%.0f" $(bc -l <<< "scale=2; ${dim[2]}*${res[2]}"))

    if [ ${fov} -gt 256 ]; then
      Info "Cropping structural image to 256 for surface reconstruction compatibility"
      crop=$(bc -l <<< "scale=0; (${fov}-256)/${res[2]}")
      Do_cmd mrgrid ${T1nativepro} crop -axis 2 0,${crop} ${tmp}/space-nativepro_t1w_croped.nii.gz
      T1nativepro=${tmp}/space-nativepro_t1w_croped.nii.gz
    else
      crop="FALSE"
    fi

    # Define SUBJECTS_DIR for surface processing as a global variable
    export SUBJECTS_DIR=${tmp} # Will work on a temporal directory

    # Recontruction method
    if [[ "$FastSurfer" == "TRUE" ]]; then
        Info "Running FastSurfer"
        Do_cmd run_fastsurfer.sh --t1 "${T1nativepro}" \
                  --sid "${idBIDS}" --sd "${SUBJECTS_DIR}" \
                  --fs_license "${fs_licence}" \
                  --parallel --threads "${threads}" --py python3.7 --no_cuda \
		              --surfreg
    else
        Info "Running Freesurfer"
        Do_cmd recon-all -cm -all -parallel -openmp ${threads} -i "${T1nativepro}" -s "${idBIDS}"
    fi

    # Copy the freesurfer log to our MICA-log Directory
    Do_cmd cp "${SUBJECTS_DIR}/${idBIDS}/scripts/recon-all.log" "${dir_logs}/recon-all.log"

    # Copy results from TMP to deviratives/SUBJECTS_DIR directory
    Do_cmd cp -r "${SUBJECTS_DIR}/${idBIDS}" "${dir_surf}"

    Info "Check log file:\n\t\t\t ${dir_logs}/recon-all.log"
fi

# -----------------------------------------------------------------------------------------------
# Check proc_surf status
if [[ -f "${dir_logs}/recon-all.log" ]] && grep -q "finished without error" "${dir_logs}/recon-all.log"; then ((Nsteps++)); fi

# Create json file for T1native
proc_surf_json="${proc_struct}/${idBIDS}_proc_surf-${recon}.json"
json_surf "${T1nativepro}" "${dir_surf}" "${recon}" "${proc_surf_json}"

# Notification of completition
micapipe_completition_status proc_surf
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_surf" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus_json "${id}" "${SES/ses-/}" "proc_surf" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
