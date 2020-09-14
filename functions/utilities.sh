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
  BIDS=$1
  id=$2
  out=$3

  id=${id/sub-/}
  subject=sub-${id}
  subject_dir=$out/${subject}/ses-pre     # Output directory
  subject_bids=${BIDS}/${subject}/ses-pre # Input BIDS directory

  # Structural directories derivatives/
  proc_struct=$subject_dir/proc_struct # structural processing directory
  	 dir_first=$proc_struct/first      # FSL first
  	 dir_volum=$proc_struct/volumetric # Cortical segmentations
  	 dir_patc=$proc_struct/surfpatch   # Surfpatch
  	 dir_surf=$proc_struct/surfaces    # surfaces
  			     dir_freesurfer=${dir_surf}/${id}     # freesurfer dir
  			     dir_conte69=${dir_surf}/conte69    # conte69
  proc_dwi=$subject_dir/proc_dwi      # DWI processing directory
  proc_rsfmri=$subject_dir/proc_rsfmri
    rsfmri_ICA=$proc_rsfmri/ICA_MELODIC
    rsfmri_volum=$proc_rsfmri/volumetric
    rsfmri_surf=$proc_rsfmri/surfaces
  dir_warp=$subject_dir/xfms              # Transformation matrices
  dir_logs=$subject_dir/logs          # directory with log files

  # post structural Files (the resolution might vary depending on the dataset)
  if [ -f ${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz ]; then
    T1nativepro=${proc_struct}/${id}_t1w_*mm_nativepro.nii.gz
    T1nativepro_brain=${proc_struct}/${id}_t1w_*mm_nativepro_brain.nii.gz
    T1freesurfr=${dir_freesurfer}/mri/T1.mgz
    T15ttgen=${proc_struct}/${id}_t1w_*mm_nativepro_5TT.nii.gz
    T1fast_seg=$proc_struct/first/${id}_t1w_*mm_nativepro_all_fast_firstseg.nii.gz
    res=`mrinfo ${T1nativepro} -spacing | awk '{printf "%.1f\n", $2}'`
  fi

  # rsfMRI processing
  mainScan=${subject_bids}/func/${subject}_ses-pre_task-rest_acq-AP_bold.nii.gz           # Main rsfMRI scan
  mainPhaseScan=${subject_bids}/func/${subject}_ses-pre_task-rest_acq-APse_bold.nii.gz    # main phase scan
  reversePhaseScan=${subject_bids}/func/${subject}_ses-pre_task-rest_acq-PAse_bold.nii.gz # Reverse phase scan

  # BIDS Files
  bids_T1ws=(`ls ${subject_bids}/anat/*T1w.nii*`)
  bids_T1map=(`ls ${subject_bids}/anat/*mp2rage*.nii*`)
  bids_inv1=(`ls ${subject_bids}/anat/*inv1*.nii*`)
  bids_dwis=(`ls ${subject_bids}/dwi/*dwi.nii*`)
  dwi_PAreverse=${subject_bids}/dwi/${id}_ses-pre_acq-PA_dir-*_dwi.nii.gz

  # Utilities
  # -----------------------------------------------------------------------------------------------#
  #   Define UTILITIES directories
  scriptDir=${MICAPIPE}/functions
  # Directory with the templates for the processing
  export util_MNIvolumes=${MICAPIPE}/MNI152Volumes
  # Directory with all the parcellations
  export util_parcelations=${MICAPIPE}/parcellations
  # Directory with the resampled freesurfer surfaces
  export util_surface=${MICAPIPE}/surfaces # utilities/resample_fsaverage

}

bids_print.variables() {
  # This functions prints BIDS variables names
  # IF they exist
  Info "Inputs:"
  Note "id   =" $id
  Note "BIDS =" $BIDS
  Note "out  =" $out

  Info "BIDS naming:"
  Note "subject_bids =" $subject_bids
  Note "bids_T1ws    =" "N-${#bids_T1ws[@]}, $bids_T1ws"
  Note "bids_dwis    =" "N-${#bids_dwis[@]}, $bids_dwis"
  Note "subject      =" $subject
  Note "subject_dir  =" $subject_dir
  Note "proc_struct  =" $proc_struct
  Note "dir_warp     =" $dir_warp
  Note "logs         =" $dir_logs

  Info "Utilities:"
  Note "util_MNIvolumes   =" $util_MNIvolumes
  Note "util_parcelations =" $util_parcelations
  Note "util_surface      =" $util_surface
}

bids_print.variables-post() {
  # This functions prints BIDS variables names and files if found
  Info "Inputs:"
  Note "id   =" $id
  Note "BIDS =" $BIDS
  Note "out  =" $out

  Info "mica-pipe variables:"
  Note "T1 nativepro    =" `ls $T1nativepro`
  Note "T1 5tt          =" `ls $T15ttgen`
  Note "T1 fast_all     =" `ls $T1fast_seg`
  Note "T1 resolution   =" $res

  Info "mica-pipe directories:"
  Note "subject_dir     =" $subject_dir
  Note "proc_struct     =" $proc_struct
  Note "dir_freesurfer  =" $dir_freesurfer
  Note "dir_conte69     =" $dir_conte69
  Note "dir_volum       =" $dir_volum
  Note "dir_warp        =" $dir_warp
  Note "logs            =" $dir_logs

  Info "Utilities:"
  Note "util_MNIvolumes   =" $util_MNIvolumes
  Note "util_parcelations =" $util_parcelations
  Note "util_surface      =" $util_surface
}

bids_print.variables-dwi() {
  # This functions prints BIDS variables names and files if found
  Info "mica-pipe variables for rsfMRI processing:"
  Note "T1 nativepro    =" `find $T1nativepro`
  Note "T1 5tt          =" `find $T15ttgen`
  Note "T1 fast_all     =" `find $T1fast_seg`
  Note "T1 resolution   =" $res
}

bids_print.variables-rsfmri() {
  # This functions prints BIDS variables names and files if found
  Info "mica-pipe variables for DWI processing:"
  Note "T1 nativepro       =" `find $T1nativepro`
  Note "Main rsfMRI        =" `find $mainScan`
  Note "Main phase scan    =" `find $mainPhaseScan`
  Note "Main reverse phase =" `find $reversePhaseScan`
}

t1w_str() {
  # This function aims to create a NAME strig with homogeneous format
  # NEW name format:
  #      <id>_<source>_<resolution>_<brain-space>_<class>_<extras>.<extension>
  #
  id=$1
  t1w_full=$2
  space=$3
  res=`mrinfo ${bids_T1ws[i]} -spacing | awk '{printf "%.1f\n", $2}'`
  echo ${id}_t1w_${res}mm_${space}${run}
}

register_QC() {
  f=$1 # back
  m=$2 # red border
  f_str=${f/"sub-${id}_ses-pre_"/}
  m_str=${m/"sub-${id}_ses-pre_"/}
  nom=${m_str/.nii.gz/}_in_${f_str/.nii.gz/}
  QC=tmp_QCreg-
  QCpng=${QC}${nom}.png
  QCjpg=${id}_${nom}.jpg
  slicer $m $f -a ${QCpng}
  montage -quality 100 -fill white -label '' ${QCpng} -tile 1x1 -background '#000000' -geometry '640' -title ${QCjpg/.jpg/} ${QCjpg}
  rm ${QC}*.png
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
echo -e "\t\t$1\t\033[38;5;122m$2\033[0m"
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
str="`whoami` @ `uname -n` `date`"
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
      nextarg=`expr ${l_index} + 1`
      eval logfile=\${${nextarg}}
      arg=""
      l_index=$[${l_index}+1]
    fi
    l_command="${l_command}${l_sep}${arg}"
    l_sep=" "
    l_index=$[${l_index}+1]
   done
if [[ ${quiet} != TRUE ]]; then echo -e "\033[38;5;118m\n${str}:\nCOMMAND -->  \033[38;5;122m${l_command}  \033[0m"; fi
if [ -z $TEST ]; then $l_command; fi
}

cmd() {
text=$1
if [[ ${quiet} != TRUE ]]; then echo -e "\033[38;5;118mCOMMAND -->  \033[38;5;122m${text}  \033[0m"; fi
eval $text
}
