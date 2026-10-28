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
FreeSurfer=$9
fs_licence=${10}
t1=${11}
AltFS=${12}
AltFS_TTA=${13}
PROC=${14}
here=$(pwd)
export OMP_NUM_THREADS=$threads

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check inputs: Nativepro T1 or custom T1
if [[ "$surfdir" == "FALSE" ]]; then
  if [[ "$t1" != "DEFAULT" ]]; then
      if [ ! -f "${t1}" ]; then Error "The provided T1 volume does not exist:\n\t${t1}"; exit; fi
      t1_2proc=${t1}
    else
      if [ ! -f "${T1nativepro}" ]; then Error "Subject $id doesn't have T1_nativepro"; exit; fi
      t1_2proc=${T1nativepro}
  fi
fi

# Surface Directory
if [[ "$FreeSurfer" == "FALSE" ]]; then recon="fastsurfer"; else recon="freesurfer"; fi
set_surface_directory "${recon}"
Note "Surface software" "${recon}"

# -alt_fs: FastSurfer -> FastScrubber brain mask -> FastSurfer surfaces (experimental)
if [[ "$AltFS" == "TRUE" ]]; then
  Note "Alternative FastSurfer:" "FastScrubber brain mask (experimental), TTA: ${AltFS_TTA}"
  # MICAPIPE_FASTSCRUBBER overrides the command, e.g. for an installation outside the container
  read -r -a fastscrubber_cmd <<< "${MICAPIPE_FASTSCRUBBER:-fastscrubber}"
  if ! command -v "${fastscrubber_cmd[0]}" > /dev/null; then
    Error "-alt_fs needs FastScrubber: '${fastscrubber_cmd[0]}' was not found"; exit 1
  fi
  if [[ -L "${dir_subjsurf}" ]]; then
    Error "-alt_fs rebuilds the surfaces in place, but ${dir_subjsurf} links to another directory.\n\t\tRemove the link and pass that directory with -surf_dir: it will be copied, not modified"; exit 1
  fi
fi

# FastSurfer runs in its own environment: CUDA in the micapipe CUDA container, CPU otherwise
fastsurfer_env=fastsurfer_cpu
if [[ -d "/opt/environments/fastsurfer_cuda" ]]; then fastsurfer_env=fastsurfer_cuda; fi
function run_fastsurfer() {
  if [[ -d "/opt/environments/${fastsurfer_env}" ]]; then
      fastsurfer_runner=(micapipe_env "${fastsurfer_env}")
  else
      source activate "${fastsurfer_env}"
      fastsurfer_runner=()
  fi
  "${fastsurfer_runner[@]}" run_fastsurfer.sh --allow_root \
   --fs_license "${fs_licence}" \
   --sid "${idBIDS}" --sd "${dir_surf}" \
   --parallel --threads "${threads}" "$@"
  local status=$?
  if [[ ! -d "/opt/environments/${fastsurfer_env}" ]]; then source activate micapipe; fi
  return $status
}

# Surface Directories
if [ ! -d "${dir_surf}" ]; then mkdir -m 777 "${dir_surf}"; fi
if [ ! -d "${dir_surf}/fsaverage5" ]; then Do_cmd mkdir -m 777 -p "${dir_surf}"/fsaverage5/surf; cp  "$FREESURFER_HOME/subjects/fsaverage5/surf/"* "${dir_surf}"/fsaverage5/surf; fi
if [ ! -d "${dir_surf}/fsLR-32k" ]; then Do_cmd mkdir -m 777 -p "${dir_surf}"/fsLR-32k/surf; cp "${MICAPIPE}"/surfaces/fsLR-32k*.gii "${dir_surf}"/fsLR-32k/surf; fi
if [ ! -d "${dir_surf}/fsLR-5k" ]; then Do_cmd mkdir -m 777 -p "${dir_surf}"/fsLR-5k/surf; cp "${MICAPIPE}"/surfaces/fsLR-5k*.gii "${dir_surf}"/fsLR-5k/surf; fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_surf-${recon}.json"
micapipe_check_json_status "${module_json}" "proc_surf"

#------------------------------------------------------------------------------#
Title "Surface processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
# print the names on the terminal
bids_print.variables
Note "Preprocessed surface directory:" "$surfdir"

# # Create script specific temp directory
tmp="${tmpDir}/${id}_micapipe_proc-surf_${RANDOM}"
Info "Processing surfaces"
Note "Saving temporal dir:" "$nocleanup"
Note "Temporal dir:" "${tmp}"
Note "Parallel processing:" "$threads threads"
Note "fs licence:" "${fs_licence}"
Note "t1 to process:" "${t1_2proc}"
Do_cmd mkdir -p "${tmp}/nii"
export FS_LICENSE="${fs_licence}"

#	Timer and steps progress
aloita=$(date +%s)
N=0
Nsteps=0

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# IF SURFACE directory is provided create a symbolic link
if [[ "$surfdir" != "FALSE" ]]; then ((N++))
    if [[ -d "$surfdir" ]]; then
        if [[ -d "$dir_subjsurf" ]]; then
            # dir_subjsurf IS under ${out}/fastsurfer
            Info "Current surface directory has a compatible naming and structure with micapipe"
        elif [[ "$AltFS" == "TRUE" ]]; then
            # -alt_fs rebuilds the surfaces, so work on a copy and leave the provided directory untouched
            if [[ ! -f "${surfdir}/mri/aparc.DKTatlas+aseg.deep.mgz" ]]; then
                Error "-alt_fs needs a FastSurfer directory, but ${surfdir} has no mri/aparc.DKTatlas+aseg.deep.mgz"
                cleanup "$tmp" "$nocleanup" "$here"; exit
            fi
            Info "Copying the provided FastSurfer directory for -alt_fs"
            Do_cmd cp -a "$surfdir" "$dir_subjsurf"
        else
            Info "Creating links for micapipe compatibility"
            # dir_subjsurf IS NOT at ${out}/fastsurfer
            Do_cmd ln -s "$surfdir" "$dir_subjsurf"
        fi
    elif [[ ! -d "$surfdir" ]]; then
        Error "The provided surface directory does not exist: $surfdir"
        cleanup "$tmp" "$nocleanup" "$here"; exit
    fi
    Do_cmd cp "${dir_subjsurf}/scripts/recon-all.log" "${dir_logs}/recon-all.log"
    t1_2proc="${dir_subjsurf}/mri/orig.mgz"
    # If mri/T1.mgz is not generated by fastsurfer cp mri/norm.mgz
    T1surf="${dir_subjsurf}"/mri/T1.mgz
    if [[ ! -f "${T1surf}" ]]; then Do_cmd cp "${dir_subjsurf}"/mri/norm.mgz "${dir_subjsurf}"/mri/T1.mgz; fi

# If not, get ready to run the surface reconstruccion
elif [[ "$surfdir" == "FALSE" ]]; then ((N++))
    # Define SUBJECTS_DIR for surface processing as a global variable
    export SUBJECTS_DIR="${tmp}" # Will work on a temporal directory

    # Recontruction method
    Note "PROC" "$PROC"
    if [[ "$recon" == "fastsurfer" ]]; then
        Do_cmd mkdir -m 777 -p "${dir_surf}/${idBIDS}"
          Info "FastSurfer: running ${fastsurfer_env} environment"
          run_fastsurfer --t1 "${t1_2proc}"
          chmod aug+wr -R "${dir_surf}/${idBIDS}"
          Do_cmd cp "${dir_surf}/${idBIDS}/scripts/recon-all.log" "${dir_logs}/recon-all.log"
    else
        Info "Running FreeSurfer: conform volume to minimum"

        # FIX FOV greater than 256
        t1nii="${t1_2proc}"
        dim=($(mrinfo "${t1nii}" -size))
        res=($(mrinfo "${t1nii}" -spacing))
        fov=$(printf "%.0f" $(bc -l <<< "scale=2; ${dim[2]}*${res[2]}"))
        if [ ${fov} -gt 256 ]; then
          Info "Cropping structural image to 256 for surface reconstruction compatibility"
          crop=$(bc -l <<< "scale=0; (${fov}-256)/${res[2]}")
          Do_cmd mrgrid "${t1nii}" crop -axis 2 0,"${crop}" "${tmp}"/_T1w_croped.nii.gz
          t1nii=${tmp}/_T1w_croped.nii.gz
        else
          crop="FALSE"
        fi
        if (( $(echo "${res[1]} > 1" |bc -l) )); then cm=""; else cm="-cm"; fi
        # Run recon-all -autorecon1
        Do_cmd recon-all -all "${cm}" -parallel -openmp "${threads}" -i "${t1nii}" -s "${idBIDS}"
        # Copy the freesurfer log to our MICA-log Directory
        Do_cmd cp "${SUBJECTS_DIR}/${idBIDS}/scripts/recon-all.log" "${dir_logs}/recon-all.log"
        # Copy results from TMP to deviratives/SUBJECTS_DIR directory
        Do_cmd cp -r "${SUBJECTS_DIR}/${idBIDS}" "${dir_surf}"
    fi
fi

# -alt_fs: replace FastSurfer's brain mask with FastScrubber's and rebuild the surfaces from it
if [[ "$AltFS" == "TRUE" ]]; then ((N++))
    mri_dir="${dir_subjsurf}/mri"
    scripts_dir="${dir_subjsurf}/scripts"
    Info "FastScrubber: predicting the brain mask from FastSurfer's orig.mgz"
    scrub_opts=(--device auto --num_cpus "${threads}")
    if [[ "$AltFS_TTA" != "TRUE" ]]; then scrub_opts+=(--disable_tta); fi
    Do_cmd "${fastscrubber_cmd[@]}" "${mri_dir}/orig.mgz" "${tmp}/fastscrubber_mask.mgz" "${scrub_opts[@]}"
    if [[ ! -f "${tmp}/fastscrubber_mask.mgz" ]]; then
        Error "FastScrubber did not write a mask"; cleanup "$tmp" "$nocleanup" "$here"; exit
    fi
    # Keep FastSurfer's own mask and logs beside the rebuilt ones
    if [[ ! -f "${mri_dir}/mask.fastsurfer.mgz" ]]; then Do_cmd cp "${mri_dir}/mask.mgz" "${mri_dir}/mask.fastsurfer.mgz"; fi
    Do_cmd cp "${tmp}/fastscrubber_mask.mgz" "${mri_dir}/mask.mgz" && ((Nsteps++))
    # Only the rebuild's log may mark this module as finished
    for log in recon-all recon-surf; do
        if [[ -f "${scripts_dir}/${log}.log" ]]; then Do_cmd mv "${scripts_dir}/${log}.log" "${scripts_dir}/${log}.fastsurfer.log"; fi
    done
    Do_cmd rm -f "${dir_logs}/recon-all.log"
    # recon-surf.sh refuses to run over a previous surface run while these exist
    Do_cmd rm -f "${mri_dir}/wm.mgz" "${mri_dir}/aparc.DKTatlas+aseg.orig.mgz"
    Info "FastSurfer: rebuilding the surfaces from the FastScrubber mask"
    run_fastsurfer --surf_only
    chmod aug+wr -R "${dir_subjsurf}"
    Do_cmd cp "${scripts_dir}/recon-all.log" "${dir_logs}/recon-all.log"
    echo -e "{
  \"micapipeVersion\": \"${Version}\",
  \"LastRun\": \"$(date)\",
  \"Workflow\": \"FastSurfer, FastScrubber brain mask, FastSurfer surfaces (-alt_fs, experimental)\",
  \"FastScrubberOptions\": \"${scrub_opts[*]}\",
  \"TestTimeAugmentation\": \"${AltFS_TTA}\",
  \"FastSurferMask\": \"mri/mask.fastsurfer.mgz\",
  \"FastScrubberMask\": \"mri/mask.mgz\"
}" > "${scripts_dir}/micapipe_alt_fs.json"
fi
Note "Check log file:" "${dir_logs}/recon-all.log"

# -----------------------------------------------------------------------------------------------
# Check proc_surf status
if [[ -f "${dir_logs}/recon-all.log" ]] && grep -q "finished without error" "${dir_logs}/recon-all.log"; then ((Nsteps++)); fi

# Create json file for T1native
proc_surf_json="${subject_dir}/surf/${idBIDS}_proc_surf-${recon}.json"
json_surf "${t1_2proc}" "${dir_surf}" "${recon}" "${proc_surf_json}"

# Notification of completition
micapipe_completition_status proc_surf
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_surf-${recon}" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus_json "${id}" "${SES/ses-/}" "proc_surf-${recon}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
