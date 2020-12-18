#!/bin/bash
#
# MICA BIDS structural processing
#
# Utilities

bids_variables() {
  # This functions assignes variables names acording to:
  #     BIDS directory = $1
  #     participant ID = $2
  #     Out directory  = $3
  #     session        = $4
  BIDS=$1
  id=$2
  out=$3
  SES=$4

  #   Define UTILITIES directories
  export scriptDir=${MICAPIPE}/functions
  # Directory with the templates for the processing
  export util_MNIvolumes=${MICAPIPE}/MNI152Volumes
  # Directory with all the parcellations
  export util_parcelations=${MICAPIPE}/parcellations
  export util_lut=${MICAPIPE}/parcellations/lut
  # Directory with the resampled freesurfer surfaces
  export util_surface=${MICAPIPE}/surfaces # utilities/resample_fsaverage
  export util_mics=${MICAPIPE}/MICs60_T1-atlas

  export subject=sub-${id}
  export subject_dir=$out/${subject}/${SES}     # Output directory
  export subject_bids=${BIDS}/${subject}/${SES} # Input BIDS directory

  # Structural directories derivatives/
  export proc_struct=$subject_dir/proc_struct # structural processing directory
  	 export dir_first=$proc_struct/first      # FSL first
  	 export dir_volum=$proc_struct/volumetric # Cortical segmentations
  	 export dir_surf=$proc_struct/surfaces    # surfaces
  			     export dir_freesurfer=${dir_surf}/${id}  # freesurfer dir
  			     export dir_conte69=${dir_surf}/conte69   # conte69
  export proc_dwi=$subject_dir/proc_dwi               # DWI processing directory
    export dwi_cnntm=$proc_dwi/connectomes
    export autoTract_dir=$proc_dwi/auto_tract
  export proc_rsfmri=$subject_dir/proc_rsfmri
    export rsfmri_ICA=$proc_rsfmri/ICA_MELODIC
    export rsfmri_volum=$proc_rsfmri/volumetric
    export rsfmri_surf=$proc_rsfmri/surfaces
  export dir_warp=$subject_dir/xfms              # Transformation matrices
  export dir_logs=$subject_dir/logs              # directory with log files
  export dir_QC=$subject_dir/QC                  # directory with QC files
  export dir_QC_png=$subject_dir/QC/png                  # directory with QC files

  # post structural Files (the resolution might vary depending on the dataset)
  if [ -f "${proc_struct}"/"${id}"_t1w_*mm_nativepro.nii.gz ]; then
      export T1nativepro=${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz
      export T1nativepro_brain=${proc_struct}/${id}_t1w_*mm_nativepro_brain.nii.gz
      export T1freesurfr=${dir_freesurfer}/mri/T1.mgz
      export T15ttgen=${proc_struct}/${id}_t1w_*mm_nativepro_5TT.nii.gz
      export T1fast_seg=$proc_struct/first/${id}_t1w_*mm_nativepro_all_fast_firstseg.nii.gz
      export res=$(mrinfo ${T1nativepro} -spacing | awk '{printf "%.1f\n", $2}')
  fi

  # Native midsurface in gifti format
  export lh_midsurf=${dir_freesurfer}/surf/lh.midthickness.surf.gii
  export rh_midsurf=${dir_freesurfer}/surf/rh.midthickness.surf.gii

  # Registration from MNI152 to Native pro
  export T1str_nat=${id}_t1w_${res}mm_nativepro
  export mat_MNI152_SyN=${dir_warp}/${T1str_nat}_brain_to_0.8mm_MNI152_SyN_brain_    # transformation strings nativepro to MNI152_0.8mm
  export T1_MNI152_InvWarp=${mat_MNI152_SyN}1InverseWarp.nii.gz                      # Inversewarp - nativepro to MNI152_0.8mm
  export T1_MNI152_affine=${mat_MNI152_SyN}0GenericAffine.mat
  export MNI152_mask=${util_MNIvolumes}/MNI152_T1_0.8mm_brain_mask.nii.gz

  # BIDS Files: resting state
  export bids_mainScan=${subject_bids}/func/${subject}_${SES}_task-rest_acq-AP_*.nii*       # main rsfMRI scan
  export bids_mainScanJson=${subject_bids}/func/${subject}_${SES}_task-rest_acq-AP_*.json   # main rsfMRI scan json
  export bids_mainPhase=${subject_bids}/func/${subject}_${SES}_task-rest_acq-APse*.nii*     # main phase scan
  export bids_reversePhase=${subject_bids}/func/${subject}_${SES}_task-rest_acq-PAse*.nii*  # reverse phase scan

  # Resting state proc files
  export topupConfigFile=${FSLDIR}/etc/flirtsch/b02b0_1.cnf                                    # TOPUP config file default
  export icafixTraining=${MICAPIPE}/functions/MICAMTL_training_15HC_15PX.RData                 # ICA-FIX training file default

  # BIDS Files
  bids_T1ws=($(ls ${subject_bids}/anat/*T1w.nii* 2>/dev/null))
  bids_dwis=($(ls ${subject_bids}/dwi/*acq-b*_dir-*_dwi.nii* 2>/dev/null))
  export bids_T1map=${subject_bids}/anat/*mp2rage*.nii*
  export bids_inv1=${subject_bids}/anat/*inv1*.nii*
  export dwi_reverse=${subject_bids}/dwi/*_${SES}_acq-PA_dir-*_dwi.nii*
}

bids_print.variables() {
  # This functions prints BIDS variables names
  # IF they exist
  Info "mica-pipe inputs:"
  Note "id   =" "$id"
  Note "BIDS =" "$BIDS"
  Note "out  =" "$out"
  Note "ses  =" "$SES"

  Info "BIDS naming:"
  Note "subject_bids =" "$subject_bids"
  Note "bids_T1ws    =" "N-${#bids_T1ws[@]}, e.g. ${bids_T1ws[0]}"
  Note "bids_dwis    =" "N-${#bids_dwis[@]}, e.g. ${bids_dwis[0]}"
  Note "subject      =" "$subject"
  Note "subject_dir  =" "$subject_dir"
  Note "proc_struct  =" "$proc_struct"
  Note "dir_warp     =" "$dir_warp"
  Note "logs         =" "$dir_logs"

  Info "Processing directories:"
  Note "subject_dir     =" "$subject_dir"
  Note "proc_struct     =" "$proc_struct"
  Note "dir_freesurfer  =" "$dir_freesurfer"
  Note "dir_conte69     =" "$dir_conte69"
  Note "dir_volum       =" "$dir_volum"
  Note "dir_warp        =" "$dir_warp"
  Note "dir_logs        =" "$dir_logs"
  Note "dir_QC          =" "$dir_QC"

  Info "Utilities directories:"
  Note "scriptDir         =" "$scriptDir"
  Note "util_MNIvolumes   =" "$util_MNIvolumes"
  Note "util_lut          =" "$util_lut"
  Note "util_parcelations =" "$util_parcelations"
  Note "util_surface      =" "$util_surface"
  Note "util_mics         =" "$util_mics"
}

bids_print.variables-post() {
  # This functions prints BIDS variables names and files if found
  Info "Structural processing output variables:"
  Note "T1 nativepro    =" "$(find $T1nativepro 2>/dev/null)"
  Note "T1 5tt          =" "$(find $T15ttgen 2>/dev/null)"
  Note "T1 fast_all     =" "$(find $T1fast_seg 2>/dev/null)"
  Note "T1 resolution   =" "$res"
}

bids_print.variables-dwi() {
  # This functions prints BIDS variables names and files if found
  Info "mica-pipe variables for DWI processing:"
  Note "bids_dwis       =" "N-${#bids_dwis[@]}, $bids_dwis"
  Note "dwi_reverse     =" "$(find $dwi_reverse 2>/dev/null)"
  Note "proc_dwi        =" "$proc_dwi"

  Note "T1 nativepro    =" "$(find $T1nativepro 2>/dev/null)"
  Note "T1 5tt          =" "$(find $T15ttgen 2>/dev/null)"
  Note "MNI152_mask     =" "$MNI152_mask"
}

bids_print.variables-rsfmri() {
  # This functions prints BIDS variables names and files if found
  Info "mica-pipe variables for rs-fMRI processing:"
  Note "T1 nativepro       =" "$(find $T1nativepro 2>/dev/null)"
  Note "Main rsfMRI        =" "$(find ${bids_mainScan} 2>/dev/null)"
  Note "Main rsfMRI json   =" "$(find ${bids_mainScanJson} 2>/dev/null)"
  Note "Main phase scan    =" "$(find ${bids_mainPhase} 2>/dev/null)"
  Note "Main reverse phase =" "$(find ${bids_reversePhase} 2>/dev/null)"
  Note "TOPUP config file  =" $(find "$topupConfigFile" 2>/dev/null)
  Note "ICA-FIX training   =" $(find "$icafixTraining" 2>/dev/null)
}

bids_variables_unset() {
  # This function unsets all the enviromentalk variables defined by
  # bids_variables
  unset scriptDir
  unset util_MNIvolumes
  unset util_parcelations
  unset util_lut
  unset util_surface
  unset util_mics
  unset subject
  unset subject_dir
  unset subject_bids
  unset proc_struct
  unset dir_first
  unset dir_volum
  unset dir_surf
  unset dir_freesurfer
  unset dir_conte69
  unset proc_dwi
  unset dwi_cnntm
  unset proc_rsfmri
  unset rsfmri_ICA
  unset rsfmri_volum
  unset rsfmri_surf
  unset dir_warp
  unset dir_logs
  unset dir_QC
  unset dir_QC_png
  unset T1nativepro
  unset T1nativepro_brain
  unset T1freesurfr
  unset T15ttgen
  unset T1fast_seg
  unset res
  unset lh_midsurf
  unset rh_midsurf
  unset T1str_nat
  unset mat_MNI152_SyN
  unset T1_MNI152_InvWarp
  unset T1_MNI152_affine
  unset MNI152_mask
  unset bids_mainScan
  unset bids_mainScanJson
  unset bids_mainPhase
  unset bids_reversePhase
  unset topupConfigFile
  unset icafixTraining
  unset bids_T1ws
  unset bids_dwis
  unset bids_T1map
  unset bids_inv1
  unset dwi_reverse
}

t1w_str() {
  # This function aims to create a NAME strig with homogeneous format
  # NEW name format:
  #      <id>_<source>_<resolution>_<brain-space>_<class>_<extras>.<extension>
  #
  id=$1
  t1w_full=$2
  space=$3
  res=$(mrinfo "${t1w_full}" -spacing | awk '{printf "%.1f\n", $2}')
  echo "${id}_t1w_${res}mm_${space}${run}"
}

micapipe_software() {
  Info "MICA pipe - Software versions"
  Note "MRtrix3....." "$(mrinfo -version | awk 'NR==1 {print $3}')"
  Note "            " "$(which mrinfo)"
  Note "FSL........." "$(flirt -version | awk '{print $3}')"
  Note "            " "$FSLDIR"
  Note "ANFI........" "$(afni -version | awk -F ':' '{print $2}')"
  Note "            " "$(which 3dresample)"
  Note "ANTS........" "$(antsRegistration --version | awk -F ':' 'NR==1{print $2}')"
  Note "            " "$ANTSPATH"
  Note "WorkBench..." "$(wb_command -version | awk 'NR==3{print $2}')"
  Note "            " "$(which wb_command)"
  Note "FreeSurfer.." "$(recon-all -version)"
  Note "            " "$FREESURFER_HOME"
  Note "fix........." "$(which fix)"
  Note "            " "$FIXPATH"
  Note "python......" "$(python --version)"
  Note "            " "$(which python)"
  Note "R..........." "$(R --version | awk 'NR==1{print $3}')"
  Note "            " "$(which R)"
}

#---------------- FUNCTION: PRINT ERROR & Note ----------------#
# The following functions are only to print on the terminal colorful messages:
# This is optional on the pipelines
#     Error messages
#     Warning messages
#     Note messages
#     Warn messages
#     Title messages
Error() {
echo -e "\033[38;5;9m\n-------------------------------------------------------------\n\n[ ERROR ]..... $1\n
-------------------------------------------------------------\033[0m\n"
}
Note(){
# I replaced color \033[38;5;197m to \033[38;5;122m
if [[ ${quiet} != TRUE ]]; then echo -e "\t\t$1\t\033[38;5;122m$2\033[0m"; fi
}
Info() {
Col="38;5;75m" # Color code
if [[ ${quiet} != TRUE ]]; then echo  -e "\033[$Col\n[ INFO ]..... $1 \033[0m"; fi
}
Warning() {
Col="38;5;184m" # Color code
if [[ ${quiet} != TRUE ]]; then echo  -e "\033[$Col\n[ WARNING ]..... $1 \033[0m"; fi
}
Warn() {
Col="38;5;184m" # Color code
if [[ ${quiet} != TRUE ]]; then echo  -e "\033[$Col
-------------------------------------------------------------\n
[ WARNING ]..... $1
\n-------------------------------------------------------------\033[0m"; fi
}
Title() {
if [[ ${quiet} != TRUE ]]; then echo -e "\n\033[38;5;141m
-------------------------------------------------------------
\t$1
-------------------------------------------------------------\033[0m"; fi
}

#---------------- FUNCTION: PRINT COLOR COMMAND ----------------#
function Do_cmd() {
# do_cmd sends command to stdout before executing it.
str="$(whoami) @ $(uname -n) $(date)"
local l_command=""
local l_sep=" "
local l_index=1
while [ ${l_index} -le $# ]; do
    eval arg=\${$l_index}
    if [ "$arg" = "-fake" ]; then
      isFake=1
      arg=""
    fi
    if [ "$arg" = "-no_stderr" ]; then
      stderr=0
      arg=""
    fi
    if [ "$arg" == "-log" ]; then
      nextarg=$(expr "${l_index}" + 1)
      eval logfile=\${"${nextarg}"}
      arg=""
      l_index=$[${l_index}+1]
    fi
    l_command="${l_command}${l_sep}${arg}"
    l_sep=" "
    l_index=$[${l_index}+1]
   done
if [[ ${quiet} != TRUE ]]; then echo -e "\033[38;5;118m\n${str}:\nCOMMAND -->  \033[38;5;122m${l_command}  \033[0m"; fi
if [ -z "$TEST" ]; then $l_command; fi
}

cmd() {
text=$1
if [[ ${quiet} != TRUE ]]; then echo -e "\033[38;5;118mCOMMAND -->  \033[38;5;122m${text}  \033[0m"; fi
eval "$text"
}

function cleanup() {
  # ITS A TRAP! Specifically a trap that will clean the temporal directory
  # and reset the old user path upon,
  # interrupts, and termination.
  Error "something went wrong, check the logs"
  export PATH=$OLD_PATH
  unset OLD_PATH
  rm -rf "$tmp"
  if [ ! -z "$here" ]; then cd "$here"; fi
  bids_variables_unset
}
