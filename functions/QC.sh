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
  echo -e "\nMICAPIPE November 2022 ${Version}\n"
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
  \t\033[38;5;120m-tracts <int>\033[0m     : OPTIONAL Number of streamlines, where 'M' stands for millions (default=40M)

  \033[38;5;141mOPTIONS:\033[0m
  \t\033[38;5;197m-h|-help\033[0m          : Print help
  \t\033[38;5;197m-tmpDir\033[0m           : Specify location of temporary directory <path> (Default is /tmp)
  \t\033[38;5;197m-quiet\033[0m 	          : Do not print comments
  \t\033[38;5;197m-nocleanup\033[0m 	  : Do not delete temporal directory at script completion
  \t\033[38;5;197m-version\033[0m 	  : Print software version

  \033[38;5;141mUSAGE:\033[0m
      \033[38;5;141m$(basename $0)\033[0m \033[38;5;197m-sub\033[0m <subject_id> \033[38;5;197m-out\033[0m <outputDirectory> \033[38;5;197m-bids\033[0m <BIDS-directory>\n

  \033[38;5;141mDEPENDENCIES:\033[0m

  McGill University, MNI, MICA-lab, May-September 2020
  https://github.com/MICA-MNI/micapipe
  http://mica-mni.github.io/
  "
}

# Source utilities functions from MICAPIPE
MICAPIPE=$(dirname $(dirname $(realpath "$0")))
source "${MICAPIPE}/functions/utilities.sh"
umask 003

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
  -tracts)
    tracts=$2
    shift;shift
  ;;
  -mica)
    mica=TRUE
    shift
  ;;
  -nocleanup)
    nocleanup=TRUE
    shift
  ;;
  -tmpDir)
    tmpDir=$2
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

# Erase temporal files by default
if [ -z "${nocleanup}" ]; then nocleanup=FALSE; fi

# Processing
if [[ -z $PROC ]]; then export PROC="LOCAL"; fi
if [ "$mica" = "TRUE" ]; then source "${MICAPIPE}/functions/init.sh"; fi

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

procDirs=$(ls -d "${subject_dir}/anat" "${subject_dir}/dwi" "${subject_dir}/func" | wc -l)
if [ "${procDirs}" -lt 3 ]; then
Error "Wrong path to subject_dir, Did you forget to set the '-ses' flag (SINGLE by default)?:
               -ses  : $SES"
help; exit 1; fi

# Variables
parcellations=$(find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*" | sort)
workflow="${dir_QC}/${idBIDS}_desc-qc_micapipe_workflow.html"

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

if false; then
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
if false; then
for func_scan in `ls ${subject_bids}/func/${idBIDS}_task-rest_echo-*_bold.nii.gz`; do
    func_scan_mean=$(basename $func_scan | sed "s/.nii.gz/_mean.nii.gz/")
    Do_cmd fslmaths "${func_scan}" -Tmean "${tmpDir}/${func_scan_mean}"
done

for fmap_scan in `ls ${subject_bids}/fmap/${idBIDS}_acq-fmri_dir-*_epi.nii.gz`; do
    fmap_scan_mean=$(basename $fmap_scan | sed "s/.nii.gz/_mean.nii.gz/")
    Do_cmd fslmaths "${fmap_scan}" -Tmean "${tmpDir}/${fmap_scan_mean}"
done
fi

# -----------------------------------------------------------------------------------------------
# Diffusion processing
# -----------------------------------------------------------------------------------------------

# STRUCTRAL CONNECTOME --------------------------------------------------------------------------
if false; then
dwi_fod="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.mif"
Do_cmd mrconvert "$dwi_fod" -coord 3 0 -axes 0,1,2  "${tmpDir}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"

for tdi in `ls ${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-*_tdi.nii.gz`; do
    tdi_mean=$(basename $tdi | sed "s/.nii.gz/_mean.nii.gz/")
    Do_cmd mrmath "${tdi}" mean "${tmpDir}/${tdi_mean}" -axis 3
done
fi

# -----------------------------------------------------------------------------------------------
# Generate QC PDF
# -----------------------------------------------------------------------------------------------
Title "Generating QC report"

Do_cmd python "$MICAPIPE"/functions/QC.py -sub ${subject} -out ${out} -bids ${BIDS} -ses ${SES/ses-/} -tmpDir ${tmpDir}

exit
# -----------------------------------------------------------------------------------------------
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print $eri/60 | perl)

# Cleanup if processing was local
if [ -d $tmpDir ]; then
    cleanup $tmpDir $nocleanup $here
fi

Title "QC html creation ended in \033[38;5;220m $(printf "%0.3f\n" ${eri}) minutes \033[38;5;141m:
\t\tOutput file path: $QC_html"
