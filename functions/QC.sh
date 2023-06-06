#!/bin/bash
#
# MICA pipe Quality Check script
#
# This script will create a basic html file for QC of the processing
#

#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#
# ONLY for scripting and debugging:
# TEST=ON
Version="(v0.2.0 'Northern flicker')"
version() {
  echo -e "\nMICAPIPE March 2023 ${Version}\n"
}
Error() {
echo -e "\033[38;5;9m\n-------------------------------------------------------------\n\n[ ERROR ]..... $1\n
-------------------------------------------------------------\033[0m\n"
}
#---------------- FUNCTION: HELP ----------------#
help() {
  echo -e "
  \033[38;5;141mCOMMAND:\033[0m
  $(basename $0)

  \033[38;5;141mARGUMENTS:\033[0m
  \t\033[38;5;197m-sub\033[0m 	          : Subject identification
  \t\033[38;5;197m-out\033[0m 	          : Output directory for the processed files <derivatives>.
  \t\033[38;5;197m-bids\033[0m 	          : Path to BIDS directory
  \t\033[38;5;120m-ses <str>\033[0m 	  : OPTIONAL flag that indicates the session name (if omitted will manage as SINGLE session)

  \033[38;5;141mOPTIONS:\033[0m
  \t\033[38;5;197m-h|-help\033[0m          : Print help
  \t\033[38;5;197m-tmpDir\033[0m           : Specify location of temporary directory <path> (Default is /tmp)
  \t\033[38;5;197m-version\033[0m 	  : Print software version

  \033[38;5;141mUSAGE:\033[0m
      \033[38;5;141m$(basename $0)\033[0m \033[38;5;197m-sub\033[0m <subject_id> \033[38;5;197m-out\033[0m <outputDirectory> \033[38;5;197m-bids\033[0m <BIDS-directory>\n

  McGill University, MNI, MICA-lab, June, 2023
  https://github.com/MICA-MNI/micapipe
  http://mica-mni.github.io/
  "
}


# -----------------------------------------------------------------------------------------------#
#			ARGUMENTS
# Create VARIABLES
for arg in "$@"
do
  case "$arg" in
  -h|-help)
    help
    exit 1
  ;;
  -version)
    version
    exit 1
  ;;
  -sub)
    id=$2
    shift;shift
  ;;
  -out)
    out=$2
    shift;shift
  ;;
  -bids)
    BIDS=$2
    shift;shift
  ;;
  -ses)
    SES=$2
    shift;shift
  ;;
  -tmpDir)
    tmpDir=$2
    shift;shift;
  ;;
  -PROC)
    PROC=$2
    shift;shift;
  ;;
  -*)
    Error "Unknown option ${2}"
    help
    exit 1
  ;;
    esac
done

# argument check out & WARNINGS
arg=($id $out $BIDS)
if [ "${#arg[@]}" -lt 3 ]; then
Error "One or more mandatory arguments are missing:
               -sub  : $id
               -out  : $out
               -bids : $BIDS"
help; exit 1; fi

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
else
  # Source utilities functions from MICAPIPE
  MICAPIPE=$(dirname $(dirname $(realpath "$0")))
fi
source "${MICAPIPE}/functions/utilities.sh"

# Get the real path of the Inputs
out=$(realpath $out)/micapipe_v0.2.0
BIDS=$(realpath $BIDS)
id=${id/sub-/}
here=$(pwd)

# Number of session (Default is "ses-pre")
if [ -z ${SES} ]; then SES="SINGLE"; else SES="ses-${SES/ses-/}"; fi

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Exit if subject is not found
if [ ! -d "${subject_dir}" ]; then Error "$id was not found on the OUTPUT directory\n\t Check ls ${subject_dir}"; exit 1; fi

# Optional arguments number of tracts
if [ -z "${tracts}" ]; then tracts=40M; else tracts="$tracts"; fi

# Temporal directory
if [ -z "${tmpDir}" ]; then export tmpDir="/tmp/${RANDOM}_micapipe_QC_${id}"; else tmpDir=$(realpath "$tmpDir"); fi

procDirs=$(ls -d "${subject_dir}/anat" "${subject_dir}/dwi" "${subject_dir}/func" | wc -l)
if [ "${procDirs}" -lt 3 ]; then
Error "Wrong path to subject_dir, Did you forget to set the '-ses' flag (SINGLE by default)?:
               -ses  : $SES"
help; exit 1; fi

# Variables
parcellations=$(find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*" | sort)
workflow="${dir_QC}/${idBIDS}_desc-qc_micapipe_workflow.html"

qc_jsons=$(ls ${subject_dir}/QC/${idBIDS}_module-*.json 2>/dev/null | wc -l)
Note "Modules processed:" $qc_jsons
if [[ "$qc_jsons" -lt 1 ]]; then exit; fi

#------------------------------------------------------------------------------#
Title "MICAPIPE: Creating a QC rport for $idBIDS"
#	Timer
aloita=$(date +%s)

#------------------------------------------------------------------------------#
# Create files and png for QC
# Create tmp dir
if [ ! -d ${tmpDir} ]; then Do_cmd mkdir -p $tmpDir; fi

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Calculate everythin on a tmpDir dir
cd $tmpDir

Title "Generating necessary files for QC report"

# -----------------------------------------------------------------------------------------------
# Structural processing
# -----------------------------------------------------------------------------------------------

# PROC_STRUC ------------------------------------------------------------------------------------
if [ -f ${subject_dir}/QC/${idBIDS}_module-proc_structural.json ]; then
  # T1w nativepro 5 tissue segmentation (5tt)
  Do_cmd mrconvert "$T15ttgen" -coord 3 0 -axes 0,1,2  "${tmpDir}/nativepro_T1w_brain_5tt.nii.gz" -force

  # Registration: T1wnatipro to MNI152 (2mm and 0.8mm)
  xfm_proc_struc_json=${dir_warp}/${idBIDS}_transformations-proc_structural.json
  for mm in 2 0.8; do
      T1w_in_MNI=${tmpDir}/${idBIDS}_space-MNI152_${mm}_T1w_brain.nii.gz

      if [[ ${mm} == 2 ]] ; then
        transformation=$(grep transformation $xfm_proc_struc_json | awk -F '"' 'NR==3{print $4}')
      else
        transformation=$(grep transformation $xfm_proc_struc_json | awk -F '"' 'NR==1{print $4}')
      fi
      MNI152_brain="${util_MNIvolumes}/MNI152_T1_${mm}mm_brain.nii.gz"
      Do_cmd antsApplyTransforms -d 3 -v -u int -o "${T1w_in_MNI}" \
              -i "${T1nativepro_brain}" \
              -r "${MNI152_brain}" \
              "${transformation}"
  done
fi

# -----------------------------------------------------------------------------------------------
# Functional processing
# -----------------------------------------------------------------------------------------------

# PROC_FUNC -------------------------------------------------------------------------------------
func_acq=($(ls -d ${subject_bids}/func/*.nii.gz | awk -F 'func/' '{print $2}' | awk -F 'task-' '{print $2}' | awk -F '_' '{print "task-"$1}' | sort -u))
if [ -f ${subject_dir}/QC/${idBIDS}_module-proc_func-*${func_acq}*.json ]; then
  for func_scan in $(ls -d ${subject_bids}/func/${idBIDS}_${func_acq}*_bold.nii.gz); do
    func_scan_mean=$(basename $func_scan | sed "s/.nii.gz/_mean.nii.gz/")
    Do_cmd fslmaths "${func_scan}" -Tmean "${tmpDir}/${func_scan_mean}"
  done

  if [ -d ${subject_bids}/fmap/ ]; then
    mainPhase_scan=($(ls "${subject_bids}/fmap/${idBIDS}"_*AP*.nii* 2>/dev/null))
    reversePhase_scan=($(ls "${subject_bids}/fmap/${idBIDS}"_*PA*.nii* 2>/dev/null))
  else
    mainPhase_scan=${bids_mainPhase[0]}
    reversePhase_scan=${bids_reversePhase[0]}
  fi

  mainPhase_scan_mean=$(basename $mainPhase_scan | sed "s/.nii.gz/_mean.nii.gz/")
  Do_cmd fslmaths "${mainPhase_scan}" -Tmean "${tmpDir}/${mainPhase_scan_mean}"

  reversePhase_scan_mean=$(basename $reversePhase_scan | sed "s/.nii.gz/_mean.nii.gz/")
  Do_cmd fslmaths "${reversePhase_scan}" -Tmean "${tmpDir}/${reversePhase_scan_mean}"

  export default_mainPhase=${bids_mainPhase[0]}
  export default_reversePhase=${bids_reversePhase[0]}

  Do_cmd antsApplyTransforms -d 3 -v -u int -o "${T1w_in_MNI}" \
          -i "${T1nativepro_brain}" \
          -r "${MNI152_brain}" \
          "${transformation}"
fi

# -----------------------------------------------------------------------------------------------
# Diffusion processing
# -----------------------------------------------------------------------------------------------
# PROC_DWI --------------------------------------------------------------------------------------
if [ -f ${subject_dir}/QC/${idBIDS}_module-proc_dwi.json ]; then
  for dwi_scan in $(ls ${subject_bids}/dwi/${idBIDS}*.nii.gz); do
      dwi_scan_mean=$(basename $dwi_scan | sed "s/.nii.gz/_mean.nii.gz/")
      Do_cmd fslmaths "${dwi_scan}" -Tmean "${tmpDir}/${dwi_scan_mean}"
  done

  Do_cmd mrmath ${proc_dwi}/${idBIDS}_space-dwi_desc-preproc_dwi.mif mean ${tmpDir}/${idBIDS}_space-dwi_desc-preproc_dwi_mean.nii.gz -axis 3

  dwi_fod="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
  Do_cmd mrconvert "$dwi_fod" -coord 3 0 -axes 0,1,2  "${tmpDir}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz" -force

  dwi_5tt="${proc_dwi}/${idBIDS}_space-dwi_desc-5tt.nii.gz"
  Do_cmd mrconvert "$dwi_5tt" -coord 3 0 -axes 0,1,2  "${tmpDir}/${idBIDS}_space-dwi_desc-5tt.nii.gz" -force
fi

# SC --------------------------------------------------------------------------------------------
if [ $(ls ${subject_dir}/QC/${idBIDS}_module-SC-*.json 2>/dev/null | wc -l) -gt 0 ]; then
  Do_cmd mrconvert "$fod" -coord 3 0 -axes 0,1,2  "${tmpDir}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"

  for tdi in $(ls ${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-*_tdi.nii.gz); do
      tdi_mean=$(basename $tdi | sed "s/.nii.gz/_mean.nii.gz/")
      Do_cmd mrmath "${tdi}" mean "${tmpDir}/${tdi_mean}" -axis 3
  done
fi

# -----------------------------------------------------------------------------------------------
# Generate QC PDF
# -----------------------------------------------------------------------------------------------
Title "Generating QC report"

Do_cmd python "$MICAPIPE"/functions/QC.py -sub ${subject} -out ${out} -bids ${BIDS} -ses ${SES/ses-/} -tmpDir ${tmpDir}

# -----------------------------------------------------------------------------------------------
# QC notification of completion
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print $eri/60 | perl)

# Cleanup if processing was local
#if [ -d $tmpDir ]; then
#    cleanup $tmpDir $nocleanup $here
#fi

Title "QC html creation ended in \033[38;5;220m $(printf "%0.3f\n" ${eri}) minutes \033[38;5;141m:
\t\tOutput file path: $QC_html"
