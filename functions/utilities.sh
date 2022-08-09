#!/bin/bash
#
# MICA BIDS structural processing
#
# Utilities
export Version="v0.2.0 'dev'"

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

  # Handle Single Session
  if [ "$SES" == "SINGLE" ]; then
      export subject_dir=$out/${subject}     # Output directory
      export subject_bids=${BIDS}/${subject} # Input BIDS directory
      ses=""
  else
      export subject_dir=$out/${subject}/${SES}     # Output directory
      export subject_bids=${BIDS}/${subject}/${SES} # Input BIDS directory
      ses="_${SES}"
  fi
export idBIDS="${subject}${ses}"

  # Structural directories derivatives/
  export dir_surf=${out/\/micapipe/}/freesurfer    # surfaces
  	 export dir_freesurfer=${dir_surf}/${idBIDS}  # freesurfer dir
  export proc_struct=$subject_dir/anat # structural processing directory
  	 export dir_volum=$proc_struct/volumetric # Cortical segmentations
  	 export dir_conte69=${proc_struct}/surfaces/conte69   # conte69
  export proc_dwi=$subject_dir/dwi               # DWI processing directory
    export dwi_cnntm=$proc_dwi/connectomes
    export autoTract_dir=$proc_dwi/auto_tract
  export proc_func=$subject_dir/func
    export func_ICA=$proc_func/ICA_MELODIC
    export func_volum=$proc_func/volumetric
    export func_surf=$proc_func/surfaces
  export proc_asl=$subject_dir/perf # ASL processing directory
    export asl_volum=$proc_asl/volumetric
    export asl_surf=$proc_asl/surfaces
  export dir_warp=$subject_dir/xfm              # Transformation matrices
  export dir_logs=$subject_dir/logs              # directory with log files
  export dir_QC=$subject_dir/QC                  # directory with QC files
  export dir_QC_png=$subject_dir/QC/png                  # directory with QC files

  # post structural Files (the resolution might vary depending on the dataset)
  if [ -f "${proc_struct}"/"${idBIDS}"_space-nativepro_t1w.nii.gz ]; then
      export res=$(mrinfo "${proc_struct}"/"${idBIDS}"_space-nativepro_t1w.nii.gz -spacing | awk '{printf "%.1f\n", $2}')
      export T1nativepro=${proc_struct}/${idBIDS}_space-nativepro_t1w.nii.gz
      export T1nativepro_brain=${proc_struct}/${idBIDS}_space-nativepro_t1w_brain.nii.gz
      export T1nativepro_mask=${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_mask.nii.gz
      export T1freesurfr=${dir_freesurfer}/mri/T1.mgz
      export T15ttgen=${proc_struct}/${idBIDS}_space-nativepro_t1w_5TT.nii.gz
      export T1fast_seg=$proc_struct/first/${idBIDS}_space-nativepro_t1w_all_fast_firstseg.nii.gz
  fi

  # Native midsurface in gifti format
  export lh_midsurf=${dir_freesurfer}/surf/lh.midthickness.surf.gii
  export rh_midsurf=${dir_freesurfer}/surf/rh.midthickness.surf.gii

  # Registration from MNI152 to Native pro
  export T1str_nat=${idBIDS}_space-nativepro_t1w
  export mat_MNI152_SyN=${dir_warp}/${idBIDS}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_    # transformation strings nativepro to MNI152_0.8mm
  export T1_MNI152_InvWarp=${mat_MNI152_SyN}1InverseWarp.nii.gz                      # Inversewarp - nativepro to MNI152_0.8mm
  export T1_MNI152_affine=${mat_MNI152_SyN}0GenericAffine.mat
  export MNI152_mask=${util_MNIvolumes}/MNI152_T1_0.8mm_brain_mask.nii.gz

  # BIDS Files: resting state
  bids_mainScan=($(ls "${subject_bids}/func/${subject}${ses}"_task-rest_acq-AP_*bold.nii* 2>/dev/null))       # main func scan
  bids_mainScanJson=($(ls "${subject_bids}/func/${subject}${ses}"_task-rest_acq-AP_*bold.json 2>/dev/null))   # main func scan json
  bids_mainPhase=($(ls "${subject_bids}/func/${subject}${ses}"_task-rest_acq-APse_*bold.nii* 2>/dev/null))     # main phase scan
  bids_reversePhase=($(ls "${subject_bids}/func/${subject}${ses}"_task-rest_acq-PAse_*bold.nii* 2>/dev/null))  # reverse phase scan

  # Resting state proc files
  export topupConfigFile=${FSLDIR}/etc/flirtsch/b02b0_1.cnf                                    # TOPUP config file default
  export icafixTraining=${MICAPIPE}/functions/MICAMTL_training_15HC_15PX.RData                 # ICA-FIX training file default

  # BIDS files: arterial spin labelling (ASL)
  bids_asl=($(ls "$subject_bids/perf/${subject}${ses}"_*asl.nii* 2>/dev/null))
  bids_m0=($(ls "$subject_bids/perf/${subject}${ses}"_*m0scan.nii* 2>/dev/null))

  # BIDS Files
  bids_T1ws=($(ls "$subject_bids"/anat/*T1w.nii* 2>/dev/null))
  bids_dwis=($(ls "${subject_bids}/dwi/${subject}${ses}"*_dir-AP_*dwi.nii* 2>/dev/null))
  bids_T1map=$(ls "$subject_bids"/anat/*mp2rage*.nii* 2>/dev/null)
  bids_inv1=$(ls "$subject_bids"/anat/*inv1*.nii* 2>/dev/null)
  bids_inv2=$(ls "$subject_bids"/anat/*inv2*.nii* 2>/dev/null)
  bids_flair=$(ls "$subject_bids"/anat/*FLAIR*.nii* 2>/dev/null)
  dwi_reverse=($(ls "${subject_bids}/dwi/${subject}${ses}"_dir-PA_*dwi.nii* 2>/dev/null))
}

bids_print.variables() {
  # This functions prints BIDS variables names
  # IF they exist
  Info "mica-pipe inputs"
  Note "id   :" "$id"
  Note "BIDS :" "$BIDS"
  Note "out  :" "$out"
  Note "ses  :" "$SES"

  Info "BIDS naming"
  Note "subject_bids :" "$subject_bids"
  Note "bids_T1ws    :" "N-${#bids_T1ws[@]}, e.g. ${bids_T1ws[0]}"
  Note "bids_dwis    :" "N-${#bids_dwis[@]}, e.g. ${bids_dwis[0]}"
  Note "subject      :" "$subject"
  Note "subject_dir  :" "$subject_dir"
  Note "proc_struct  :" "$proc_struct"
  Note "dir_warp     :" "$dir_warp"
  Note "logs         :" "$dir_logs"

  Info "Processing directories"
  Note "subject_dir     :" "$subject_dir"
  Note "proc_struct     :" "$proc_struct"
  Note "dir_conte69     :" "$dir_conte69"
  Note "dir_volum       :" "$dir_volum"
  Note "dir_warp        :" "$dir_warp"
  Note "dir_logs        :" "$dir_logs"
  Note "dir_QC          :" "$dir_QC"
  Note "dir_freesurfer  :" "$dir_freesurfer"

  Info "Utilities directories"
  Note "scriptDir         :" "$scriptDir"
  Note "util_MNIvolumes   :" "$util_MNIvolumes"
  Note "util_lut          :" "$util_lut"
  Note "util_parcelations :" "$util_parcelations"
  Note "util_surface      :" "$util_surface"
  Note "util_mics         :" "$util_mics"
}

file.exist(){
  if [[ ! -z "${2}" ]] && [[ -f "${2}" ]]; then
    Note "$1" "$(find $2 2>/dev/null)"
  else
    Note "$1" "file not found"
  fi
}

bids_print.variables-post() {
  # This functions prints BIDS variables names and files if found
  Info "Structural processing output variables"
  Note "T1 nativepro    :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 5tt          :" "$(find "$T15ttgen" 2>/dev/null)"
  Note "T1 fast_all     :" "$(find "$T1fast_seg" 2>/dev/null)"
  Note "T1 resolution   :" "$res"
}

bids_print.variables-dwi() {
  # This functions prints BIDS variables names and files if found
  Info "mica-pipe variables for DWI processing"
  Note "proc_dwi dir    :" "$proc_dwi"
  Note "bids_dwis       :" "N-${#bids_dwis[@]}, $bids_dwis"
  Note "dwi_reverse     :" "N-${#dwi_reverse[@]}, $dwi_reverse"

  Note "T1 nativepro    :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 5tt          :" "$(find "$T15ttgen" 2>/dev/null)"
  Note "MNI152_mask     :" "$MNI152_mask"
}

bids_print.variables-func() {
  # This functions prints BIDS variables names and files if found
  Info "mica-pipe variables for rs-fMRI processing"
  Note "T1 nativepro       :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 freesurfer      :" "$(find "$T1freesurfr" 2>/dev/null)"
  for i in "${!mainScan[@]}"; do
  file.exist "mainScan${i/0/}         :" ${mainScan[i]}
  file.exist "mainScan${i/0/} json    :" ${mainScanJson[i]}
  done
  file.exist "Main phase scan    :" $mainPhaseScan
  file.exist "Main reverse phase :" $reversePhaseScan
  Note "TOPUP config file  :" $(find "$topupConfigFile" 2>/dev/null)
  Note "ICA-FIX training   :" $(find "$icafixTraining" 2>/dev/null)
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
  unset dir_volum
  unset dir_surf
  unset dir_freesurfer
  unset dir_conte69
  unset proc_dwi
  unset dwi_cnntm
  unset proc_func
  unset func_ICA
  unset func_volum
  unset func_surf
  unset proc_asl
  unset asl_volum
  unset asl_surf
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
  unset bids_asl
  unset bids_m0
  unset topupConfigFile
  unset icafixTraining
  unset bids_T1ws
  unset bids_dwis
  unset bids_T1map
  unset bids_inv1
  unset dwi_reverse
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

function micapipe_procStatus() {
  id=$1
  session=$2
  mod=$3
  outfile=$4
  grep -v "${id}, ${session}, ${mod}" "${outfile}" > "${tmp}/tmpfile" && mv "${tmp}/tmpfile" "${outfile}"
  echo "${id}, ${session}, ${mod}, ${status}, $(printf "%02d" "$Nsteps")/$(printf "%02d" "$N"), $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${outfile}"
}

function micapipe_json() {
  # Name is the name of the raw-BIDS directory
  if [ -f "${BIDS}/dataset_description.json" ]; then
    Name=$(grep Name "${BIDS}"/dataset_description.json | awk -F '"' '{print $4}')
    BIDSVersion=$(grep BIDSVersion "${BIDS}"/dataset_description.json | awk -F '"' '{print $4}')
  else
    Name="BIDS dataset_description NOT found"
    BIDSVersion="BIDS dataset_description NOT found"
  fi

  echo -e "{
    \"Name\": \"${Name}\",
    \"BIDSVersion\": \"${BIDSVersion}\",
    \"DatasetType\": \"derivative\",
    \"GeneratedBy\": [
      {
        \"Name\": \"micapipe\",
        \"Version\": \"${Version}\",
        \"Container\": {
          \"Type\": \"github\",
          \"Tag\": \"MICA-MNI/micapipe:${Version}\"
          }
      },
      {
        \"Name\": \"$(whoami)\",
        \"Workstation\": \"$(uname -n)\"
        \"LastRun\": \"$(date)\"
        \"Processing\": \"${PROC}\"
      }
    ],
    \"SourceDatasets\": [
      {
        \"DOI\": \"doi:\",
        \"URL\": \"https://micapipe.readthedocs.io/en/latest/\",
        \"Version\": \"${Version}\"
      }
    ]
  }" > "${out}/pipeline-description.json"
}

function tck_json() {
  qform=$(fslhd "$dwi_b0" | grep qto_ | awk -F "\t" '{print $2}')
  sform=$(fslhd "$dwi_b0" | grep sto_ | awk -F "\t" '{print $2}')
  Info "Creating tractography json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${8}\",
    \"inputNIFTI\": [
      {
      \"Name\": \"${fod_wmN}\",
      \"sform\": [
\"${qform}\"
      ],
      \"qform\": [
\"${sform}\"
      ],
      }
    ],
    \"Tractography\": [
      {
        \"TractographyClass\": \"local\",
        \"TractographyMethod\": \"probabilistic\",
        \"TractographyAlgorithm\": \"$1\",
        \"StepSizeUnits\": [\"mm\"],
        \"StepSize\": \"$2\",
        \"AngleCurvature\": \"$3\",
        \"cutoff\": \"$4\",
        \"maxlength\": \"$5\",
        \"minlength\": \"$6\",
        \"SeedingMethod\": \"$7\",
        \"SeedingNumberMethod\": \"${tracts}\",
        \"TerminationCriterion\": [\"reachingTissueType”],
        \"TerminationCriterionTest\": [\"ACT\"],
        \"TractographySaved\": \"${nocleanup}\"
      }
    ]
  }" > "${tckjson}"
}

function json_nativepro_t1w() {
  qform=$(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}')
  sform=$(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}')
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=$(mrinfo "$1" -transform)
  Info "Creating T1w_nativepro json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${1}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"Transform\": \"${Transform}\",
    \"sform\": [
\"${qform}\"
    ],
    \"qform\": [
\"${sform}\"
    ],
    \"inputsRawdata\": \"${3}\",
    \"anatPreproc\": [
      {
        \"Resample\": \"LPI\",
        \"Reorient\": \"fslreorient2std\",
        \"NumberofT1w\": \"$2\",
        \"BiasFieldCorrection\": \"ANTS N4BiasFieldCorrection\",
        \"WMweightedN4BFC\": \"${N4wm}\",
        \"RescaleRange\": \"0:100\"
      }
    ]
  }" > "$4"
}

function json_nativepro_mask() {
  qform=$(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}')
  sform=$(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}')
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=$(mrinfo "$1" -transform)
  Info "Creating T1natipro_brain json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${1}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"Transform\": \"${Transform}\",
    \"inputNIFTI\": [
      {
      \"Name\": \"${T1nativepro}\",
      \"sform\": [
\"${qform}\"
      ],
      \"qform\": [
\"${sform}\"
      ],
      }
    ],
    \"BinaryMask\": [
      {
        \"BinaryMask_original\": \"$2\",
        \"BinaryMask_from\": \"MNI152_0.8mm\",
        \"BinaryMask_antsApplyTransforms\": \"$4\"
      }
    ]
  }" > "$3"
}

function json_func() {
  qform=$(fslhd "$fmri_processed" | grep qto_ | awk -F "\t" '{print $2}')
  sform=$(fslhd "$fmri_processed" | grep sto_ | awk -F "\t" '{print $2}')
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"Tag\": \"${tagMRI}\",
    \"Acquisition\": \"${acq}\",
    \"Name\": \"${fmri_processed}\",
    \"sform\": [
\t\t\"${sform}\"
      ],
    \"qform\": [
  \t\t\"${sform}\"
      ],
    \"Preprocess\": [
      {
        \"MainScan\": \"${mainScan[*]}\",
        \"Resample\": \"LPI\",
        \"Reorient\": \"fslreorient2std\",
        \"MotionCorrection\": \"3dvolreg AFNI $(afni -version | awk -F ':' '{print $2}')\",
        \"MotionCorrection\": [\"${func_volum}/${idBIDS}_space-func_spikeRegressors_FD.1D\"],
        \"MainPhaseScan\": \"${mainPhaseScan}\",
        \"ReversePhaseScan\": \"${reversePhaseScan}\",
        \"TOPUP\": \"${statusTopUp}\",
        \"HighPassFilter\": \"${fmri_HP}\",
        \"Passband\": \"0.01 666\",
        \"RepetitionTime\": \"${RepetitionTime[*]}\",
        \"TotalReadoutTime\": \"${readoutTime[*]}\",
        \"EchoTime\": \"${EchoTime[*]}\",
        \"Melodic\": \"${statusMel}\",
        \"FIX\": \"${statusFIX}\",
        \"Registration\": \"${reg}\",
        \"GlobalSignalRegression\": \"${performGSR}\",
        \"CSFWMSignalRegression\": \"${performNSR}\",
        \"dropTR\": \"${dropTR}\",
        \"procStatus\": \"${status}\",
      }
    ]
  }" > "$1"
}

function json_mpc() {
  qform=$(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}')
  sform=$(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}')
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=$(mrinfo "$1" -transform)
  Info "Creating MPC json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"Class\": \"Microstructural profile covariance\",
    \"acquisition\": \"${mpc_acq}\",
    \"input\": \"${1}\",
    \"freesurferTransformation\": \"${2}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"Transform\": \"${Transform}\",
    \"sform\": [
\"${qform}\"
      ],
    \"qform\": [
\"${sform}\"
      ]
  }" > "$3"
}

function json_dwipreproc() {
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=$(mrinfo "$1" -transform)

  res_rpe=$(mrinfo "$4" -spacing)
  Size_rpe=$(mrinfo "$4" -size)
  Strides_rpe=$(mrinfo "$4" -strides)
  Offset_rpe=$(mrinfo "$4" -offset)
  Multiplier_rpe=$(mrinfo "$4" -multiplier)
  Transform_rpe=$(mrinfo "$4" -transform)

  Info "Creating DWI preproc json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"Class\": \"DWI preprocessing\",
    \"DWIpe\": [
        \"fileName\": \"${bids_dwis[*]}\",
        \"NumberOfInputs\": \"${#bids_dwis[*]}\",
        \"VoxelSizepe\": \"${res}\",
        \"Dimensionspe\": \"${Size}\",
        \"Strides\": \"${Strides}\",
        \"Offset\": \"${Offset}\",
        \"Multiplier\": \"${Multiplier}\",
        \"Transform\": \"${Transform}\"
    ],
    \"DWIrpe\": [
        \"fileName\": \"${dwi_reverse[*]}\",
        \"NumberOfInputs\": \"${#dwi_reverse[*]}\",
        \"VoxelSizepe\": \"${res_rpe}\",
        \"Dimensionspe\": \"${Size_rpe}\",
        \"Strides\": \"${Strides_rpe}\",
        \"Offset\": \"${Offset_rpe}\",
        \"Multiplier\": \"${Multiplier_rpe}\",
        \"Transform\": \"${Transform_rpe}\",
    ],
    \"Denoising\": \"Marchenko-Pastur PCA denoising, dwidenoise\",
    \"GibbsRingCorrection\": \"mrdegibbs\",
    \"dwiflspreproc\": [
        \"input\": \"${dwi_4proc}\",
        \"output\": \"${dwi_corr}\",
        \"Shells\": \"${shells[*]}\",
        \"pe_dir\": \"${pe_dir}\",
        \"ReadoutTime\": \"${ReadoutTime}\",
        \"Options\": \"${opt}\",
        \"slm\": \"linear\"
    ],
    \"B1fieldCorrection\": \"ANTS N4BiasFieldCorrection\",
    \"DWIprocessed\": \"$2\",
  }" > "$3"
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
      arg=""
    fi
    if [ "$arg" = "-no_stderr" ]; then
      arg=""
    fi
    if [ "$arg" == "-log" ]; then
      nextarg=$(("${l_index}" + 1))
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

function cleanup() {
  # This script will clean the temporal directory
  # and reset the old user path upon,
  # interrupts, and termination.
  tmp=$1
  nocleanup=$2
  here=$3
  # Clean temporal directory and temporal fsaverage5
  if [[ $nocleanup == "FALSE" ]]; then
      echo -e "Erasing temporal directory: $tmp"
      rm -Rf "$tmp" 2>/dev/null
  else
      echo -e "Mica-pipe tmp directory was not erased: \n\t\t${tmp}";
  fi
  cd "$here"
  bids_variables_unset
  if [[ ! -z "$OLD_PATH" ]]; then  export PATH=$OLD_PATH; unset OLD_PATH; else echo "OLD_PATH is unset or empty"; fi
}

function QC_proc-func() {
  outname=$1
  html="$dir_QC/${outname}"
  if [ -f "$html" ]; then rm "$html"; fi
  echo -e "            <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">mainScan</span></td>
                <td class=\"tg-8pnm\">BIDS func<br><br></td>
                <td class=\"tg-8pnm\">$(find "$mainScan" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">bids_mainScanJson</span></td>
                <td class=\"tg-8pnm\">BIDS func<br><br></td>
                <td class=\"tg-8pnm\">$(find "$mainScanJson" 2>/dev/null)</td>
              </tr>" >> "$html"
            if [ -f "$mainPhaseScan" ]; then
  echo -e "          <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">mainPhaseScan</span></td>
                <td class=\"tg-8pnm\">BIDS func<br><br></td>
                <td class=\"tg-8pnm\">$(find "$mainPhaseScan" 2>/dev/null)</td>
              </tr>"  >> "$html"
            fi
            if [ -f "$reversePhaseScan" ]; then
  echo -e "          <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">reversePhaseScan</span></td>
                <td class=\"tg-8pnm\">BIDS func<br><br></td>
                <td class=\"tg-8pnm\">$(find "$reversePhaseScan" 2>/dev/null)</td>
              </tr>"  >> "$html"
            fi
  echo -e "          <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">topupConfigFile</span></td>
                <td class=\"tg-8pnm\">Default/Defined<br><br></td>
                <td class=\"tg-8pnm\">$(find "$topupConfigFile" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">icafixTraining</span></td>
                <td class=\"tg-8pnm\">Default/Defined<br><br></td>
                <td class=\"tg-8pnm\">$(find "$icafixTraining" 2>/dev/null)</td>
              </tr>"   > "$html"
}

function QC_SC() {
  html=$1
  if [ -f "$html" ]; then rm "$html"; fi
  echo -e "          <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">fod</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$fod_wmN" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">dwi_b0</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$dwi_b0" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">mat_dwi_affine</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$mat_dwi_affine" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">dwi_5tt</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$dwi_5tt" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">T1_seg_cerebellum</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$T1_seg_cerebellum" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">T1_seg_subcortex</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$T1_seg_subcortex" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">dwi_mask</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$dwi_mask" 2>/dev/null)</td>
              </tr>
              <tr>
                <td class=\"tg-8pnm\"><span style=\"font-weight:bold\">fa</span></td>
                <td class=\"tg-8pnm\">proc-dwi<br><br></td>
                <td class=\"tg-8pnm\">$(find "$dti_FA" 2>/dev/null)</td>
              </tr>"   >> "$html"
}

function micapipe_group_QC() {
  #------------------------------------------------------------------------------#
  # Group QC html file
  here=$(pwd)
  QC_html=${out}/micapipe_progress.html
  if [ ! -d "${out}" ]; then Error "Output path does not contain a /micapipe directory:\n \t${out}/micapipe "; exit; fi
  DataName=$(grep "Name" "${out}/pipeline-description.json" | awk -F '"' 'NR==1{print $4}')
  Title "MICAPIPE: group-level Quality Control"
  table_style=" <style type=\"text/css\">\n
      .tg  {border-collapse:collapse;border-spacing:0;border-top: none;border-bottom: none;}\n
      .tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
        overflow:hidden;padding:10px 5px;word-break:normal;}\n
      .tg th{position:sticky; top:119px; border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:14px;
        font-weight:normal;overflow:hidden;padding:10px 5px;word-break:normal;}\n
      .tg .tg-lp92{background-color:#656565;border-color:#c0c0c0;color:#efefef;
        font-family:\"Lucida Console\", Monaco, monospace !important;;font-size:14px;font-weight:bold;text-align:center;
        vertical-align:top}\n
      .tg .tg-oq6h{background-color:#343434;border-color:#c0c0c0;color:#818181;font-family:\"Courier New\", Courier, monospace !important;;
        text-align:center;vertical-align:top}\n
      .tg .tg-e8zy{background-color:#343434;border-color:#c0c0c0;color:#efefef;font-family:\"Courier New\", Courier, monospace !important;;
        font-weight:bold;text-align:center;vertical-align:top}\n
      .tg .tg-sl9e{background-color:#ee7942;border-color:#c0c0c0;color:#efefef;font-family:\"Courier New\", Courier, monospace !important;;
        text-align:center;vertical-align:top}\n
      .tg .tg-8779{background-color:#039fd3;border-color:#c0c0c0;color:#efefef;
        font-family:\"Lucida Console\", Monaco, monospace !important;;font-size:12px;text-align:center;vertical-align:top}\n
    </style>"

  echo -e "<!doctype html>
  <html>
  <head>
    <meta charset=\"utf-8\">
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
    <title>micapipe progress</title>

    <!-- Headings format -->
    <style>
      h1 {color:#e9e9f7;font-family:\"Lucida Console\", Monaco, monospace !important;;text-align:center;}
      h2 {color:#e9e9f7;font-family:\"Lucida Console\", Monaco, monospace !important;;}
      a {color:#efefef;font-family:\"Lucida Console\";}
    </style>

    <!-- Sticky image -->
    <style>
    img.sticky {
      position: -webkit-sticky;
      position: sticky;
      top: 0;
      width: 250px;
      height: 120px;
    }
    </style>

    <style>
    body {
      background-color: rgba(5, 5, 5, 0.95);
    }
    </style>

  </head>
  <body>
    <img class=\"sticky\", id=\"top\" src=\"${MICAPIPE}/docs/figures/micapipe_long.png\" style=\"width:100%\"  alt=\"micapipe\"> " > "$QC_html"

  #------------------------------------------------------------------------------#
  # MICAPIPE progress table
  echo -e "
  <h1>Pipeline progress</h1>
  <h2>Database: ${DataName}</h2>
  <h2>Last run: $(date)</h2>"  >> "$QC_html"

  echo -e "$table_style" >> "$QC_html"
  echo -e "
  <div class=\"container\">
      <table class=\"tg\">
        <thead>
          <tr>
            <th class=\"tg-lp92\">Subject</th>
            <th class=\"tg-lp92\">Session</th>
            <th class=\"tg-lp92\">Subject QC</th>
            <th class=\"tg-lp92\">proc_structural</th>
            <th class=\"tg-lp92\">proc_freesurfer</th>
            <th class=\"tg-lp92\">post_structural</th>
            <th class=\"tg-lp92\">Morphology</th>
            <th class=\"tg-lp92\">GD</th>
            <th class=\"tg-lp92\">proc_dwi</th>
            <th class=\"tg-lp92\">SC</th>
            <th class=\"tg-lp92\">proc_func</th>
            <th class=\"tg-lp92\">MPC</th>
          </tr>
        </thead>
        <tbody>" >> "$QC_html"

  function fill.table() {
    sub_ses=$1
    sub_html=$2
    echo -e "\n        <tr>" >> "$QC_html"
    echo -e "              <td class=\"tg-e8zy\"><span style=\"font-weight:bold\">sub-${Nsub}</span></td>
              <td class=\"tg-e8zy\"><span style=\"font-weight:bold\">${sub_ses}</span></td>" >> "$QC_html"
    if [ -f "$sub_html" ]; then
        echo -e "              <td class=\"tg-8779\"><a href=\"${sub_html}\">Click here</a></td>" >> "$QC_html"
    else
        echo -e "              <td class=\"tg-oq6h\">Not processed<br><br></td>" >> "$QC_html"
    fi
    for module in proc_structural proc_freesurfer post_structural Morphology GD proc_dwi SC proc_func MPC; do
        Status=$(grep "${Nsub}, ${sub_ses/ses-/}, ${module}" "${pipecsv}" | awk -F ", " '{print $4}')
        Steps=$(grep "${Nsub}, ${sub_ses/ses-/}, ${module}" "${pipecsv}" | awk -F ", " '{print $5}')
        Date=$(grep "${Nsub}, ${sub_ses/ses-/}, ${module}" "${pipecsv}" | awk -F ", " '{print $8}')
        RunV=$(grep "${Nsub}, ${sub_ses/ses-/}, ${module}" "${pipecsv}" | awk -F ", " '{print $11}')
        if [[ "$Status" == "COMPLETED" ]]; then
            echo -e "              <td class=\"tg-8779\">${Steps}<br>${Date}<br>${RunV}</td>" >> "$QC_html"
        elif [[ "$Status" == "INCOMPLETE" ]]; then
            echo -e "              <td class=\"tg-sl9e\">${Steps}<br>${Date}<br>${RunV}</td>" >> "$QC_html"
        else
            echo -e "              <td class=\"tg-oq6h\">Not processed<br><br></td>" >> "$QC_html"
        fi
    done
    echo -e "        </tr>" >> "$QC_html"
  }

  pipecsv=${out}/micapipe_processed_sub.csv
  cd "$out"
  for Subj in sub*; do
      Info "Processing $Subj"
      Nsub=$(echo "${Subj/sub-/}" | awk -F '/' '{print $1}')
      NumSes=($(ls "$Subj"/ses-* 2>/dev/null | wc -l))
      if [[ "$NumSes" -eq 0 ]]; then
          fill.table "SINGLE" "${out}/${Subj}/QC/sub-${Nsub}_micapipe_qc.html"
      elif [[ "$NumSes"  -gt 0 ]]; then
        for SubSes in "$Subj"/ses-*; do
            Nses=$(echo "${SubSes/ses-/}" | awk -F '/' '{print $2}')
            fill.table "ses-${Nses}" "${out}/${SubSes}/QC/sub-${Nsub}_ses-${Nses}_micapipe_qc.html"
        done
      fi
  done

  echo -e "      </tbody>
      </table>
</div>" >> "$QC_html"
  #------------------------------------------------------------------------------#
  # closes html document
  echo "
  </body>
  </html>" >> "$QC_html"
  cd "$here"
  Title "micapipe group-level QC ended succesfully\033[38;5;141m:
  \t\tOutput to file: $QC_html"
}
