#!/bin/bash
#
# MICA BIDS structural processing
#
# Utilities
export Version="v0.2.3 'Northern flicker'"
export MRTRIX_QUIET=TRUE
export CUDA_VISIBLE_DEVICES=''
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
  umask 002

  #   Define UTILITIES directories
  export scriptDir=${MICAPIPE}/functions
  # Directory with the templates for the processing
  export util_MNIvolumes=${MICAPIPE}/MNI152Volumes
  # Directory with all the parcellations
  export util_parcelations=${MICAPIPE}/parcellations
  export util_lut=${MICAPIPE}/parcellations/lut
  # Directory with the resampled surfaces
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
  export proc_struct=$subject_dir/anat # structural processing directory
  	 export dir_volum=$subject_dir/parc # Cortical segmentations
  	 export dir_conte69=$subject_dir/surf   # conte69
     export dir_maps=$subject_dir/maps
  export proc_dwi=$subject_dir/dwi               # DWI processing directory
    export dwi_cnntm=$proc_dwi/connectomes
    export autoTract_dir=$proc_dwi/auto_tract
  export proc_func=$subject_dir/func
    export func_ICA=$proc_func/ICA_MELODIC
    export func_volum=$proc_func/volumetric
    export func_surf=$proc_func/surf
  export proc_asl=$subject_dir/perf # ASL processing directory
    export asl_volum=$proc_asl/volumetric
    export asl_surf=$proc_asl/surf
  export dir_warp=$subject_dir/xfm              # Transformation matrices
  export dir_logs=$subject_dir/logs              # directory with log files
  export dir_QC=$subject_dir/QC                  # directory with QC files

  # post structural Files (the resolution might vary depending on the dataset)
  if [ -f "${proc_struct}"/"${idBIDS}"_space-nativepro_T1w.nii.gz ]; then
      export res=$(mrinfo "${proc_struct}"/"${idBIDS}"_space-nativepro_T1w.nii.gz -spacing | awk '{printf "%.1f\n", $2}')
      export T1nativepro=${proc_struct}/${idBIDS}_space-nativepro_T1w.nii.gz
      export T1nativepro_brain=${proc_struct}/${idBIDS}_space-nativepro_T1w_brain.nii.gz
      export T1nativepro_mask=${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_mask.nii.gz
      export T15ttgen=${proc_struct}/${idBIDS}_space-nativepro_T1w_5tt.nii.gz
      export T1fast_seg=$subject_dir/parc/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz
  fi

  # Registration from MNI152 to Native pro
  export T1str_nat=${idBIDS}_space-nativepro_T1w
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
  bids_T1map=$(ls "$subject_bids"/anat/*mp2rage*T1map.nii* 2>/dev/null)
  bids_inv1=$(ls "$subject_bids"/anat/*inv1*T1map.nii* 2>/dev/null)
  bids_inv2=$(ls "$subject_bids"/anat/*inv2*T1map.nii* 2>/dev/null)
  bids_flair=$(ls "$subject_bids"/anat/*FLAIR*.nii* 2>/dev/null)
  dwi_reverse=($(ls "${subject_bids}/dwi/${subject}${ses}"_dir-PA_*dwi.nii* 2>/dev/null))
}

set_surface_directory() {
  local recon=${1}
  export dir_surf=${out/\/micapipe_v0.2.0/}/${recon}    # surf
  export dir_subjsurf=${dir_surf}/${idBIDS}  # Subject surface dir
  export T1surf=${dir_subjsurf}/mri/orig.mgz

  # Native midsurface in gifti format
  export lh_midsurf="${dir_conte69}/${idBIDS}_hemi-L_surf-fsnative_label-midthickness.surf.gii"
  export rh_midsurf="${dir_conte69}/${idBIDS}_hemi-R_surf-fsnative_label-midthickness.surf.gii"
}

bids_print.variables() {
  # This functions prints BIDS variables names
  # IF they exist
  Info "micapipe inputs"
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
  Note "dir_subjsurf  :" "$dir_subjsurf"

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

bids_print.variables-structural() {
  Info "Main structural image(s) for processing"
  Note "T1w string(s):" "$T1wStr, N=${Nimgs}"
  if [[ "${UNI}" == "TRUE" ]]; then
    Note "UNI" "${bids_T1ws[0]}"
    Note "INV1" "${bids_inv1}"
    Note "INV2" "${bids_inv2}"
    Note "Multiplying Factor" "${MF}"
  fi
  Note "mp2rage-UNI" "${UNI}"
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
  Info "Variables for DWI processing"
  Note "proc_dwi dir    :" "$proc_dwi"
  Note "bids_dwis       :" "N-${#bids_dwis[@]}, $bids_dwis"
  Note "dwi_reverse     :" "N-${#dwi_reverse[@]}, $dwi_reverse"

  Note "T1 nativepro    :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 5tt          :" "$(find "$T15ttgen" 2>/dev/null)"
  Note "MNI152_mask     :" "$MNI152_mask"
}

bids_print.variables-func() {
  # This functions prints BIDS variables names and files if found
  Info "Variables for functional processing"
  Note "T1 nativepro       :" "$(find "$T1nativepro" 2>/dev/null)"
  Note "T1 surface      :" "$(find "$T1surf" 2>/dev/null)"
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
  unset dir_subjsurf
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
  unset T1surf
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
  Note "python......" "$(python --version | awk 'NR==1{print $2}')"
  Note "            " "$(which python)"
  Note "R..........." "$(R --version | awk 'NR==1{print $3}')"
  Note "            " "$(which R)"
  Note "conda......." "$(conda --version)"
  Note "            " "$(which conda)"
}

function micapipe_procStatus() {
  id=$1
  session=$2
  mod=$3
  outfile=$4
  grep -v "${id}, ${session}, ${mod}" "${outfile}" > "${tmp}/tmpfile" && mv "${tmp}/tmpfile" "${outfile}"
  echo "${id}, ${session}, ${mod}, ${status}, $(printf "%02d" "$Nsteps")/$(printf "%02d" "$N"), $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${outfile}"
}

function micapipe_completition_status() {
  if [ -z "${2}" ]; then logaqc=""; else logaqc="${2}"; fi
    # Processing time
    lopuu=$(date +%s)
    eri=$(echo "$lopuu - $aloita" | bc)
    eri=$(echo print "$eri"/60 | perl)

    # Print logs
    if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
    Title "${1} processing ended in $(printf "%0.3f\n" "$eri") minutes:\n\tlogs:
    \tSteps completed : $(printf "%02d" "$Nsteps")/$(printf "%02d" "$N")
    \tStatus          : ${status}
    \tCheck logs      : $(ls "$dir_logs"/${1}_*"${logaqc}.txt")"
}

function micapipe_procStatus_json() {
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"Module\": \"${3}\",
    \"Subject\": \"${1}\",
    \"Session\": \"${2}\",
    \"Status\": \"${status}\",
    \"Progress\": \"$(printf "%02d" "$Nsteps")/$(printf "%02d" "$N")\",
    \"User\": \"$(whoami)\",
    \"Workstation\": \"$(uname -n)\",
    \"Date\": \"$(date)\",
    \"Processing.time\": \"$(printf "%0.3f\n" "$eri")\",
    \"Processing\": \"${PROC}\",
    \"Threads\": \"${threads}\"
  }" > "${4}"
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
    \"GeneratedBy\": [{
        \"Name\": \"micapipe\",
        \"Version\": \"${Version}\",
        \"Reference\": \"Raúl R. Cruces, Jessica Royer, Peer Herholz, Sara Larivière, Reinder Vos de Wael, Casey Paquola, Oualid Benkarim, Bo-yong Park, Janie Degré-Pelletier, Mark Nelson, Jordan DeKraker, Ilana Leppert, Christine Tardif, Jean-Baptiste Poline, Luis Concha, Boris C. Bernhardt. (2022) Micapipe: a pipeline for multimodal neuroimaging and connectome analysis. NeuroImage, 2022, 119612, ISSN 1053-8119.\",
        \"DOI\": \"https://doi.org/10.1016/j.neuroimage.2022.119612\",
        \"URL\": \"https://micapipe.readthedocs.io/en/latest\",
        \"GitHub\": \"https://github.com/MICA-MNI/micapipe\",
        \"Container\": {
          \"Type\": \"docker\",
          \"Tag\": \"micalab/micapipe:$(echo ${Version} | awk '{print $1}')\"
        },
        \"RunBy\": \"$(whoami)\",
        \"Workstation\": \"$(uname -n)\",
        \"LastRun\": \"$(date)\",
        \"Processing\": \"${PROC}\"
      }]
  }" > "${out}/dataset_description.json"
}

function micapipe_check_json_status() {
  local mod_json="${1}"
  local mod_func="${2}"
  if [ -f "${module_json}" ]; then
    status=$(grep "Status" "${module_json}" | awk -F '"' '{print $4}')
    if [ "$status" == "COMPLETED" ]; then
    Note "${mod_func} json" "${module_json}"
    Warning "Subject ${idBIDS} has been processed with -${mod_func}
              If you want to re-run this step again, first erase all the outputs with:

              micapipe_cleanup -sub <subject_id> -out <derivatives> -bids <BIDS_dir> -${mod_func}"; exit
    else
        Info "${mod_func} is ${status}, processing will continute"
    fi
  fi
}

function micapipe_check_dependency() {
  local mod_str="${1}"
  local module_json="${2}"
  if [[ ! -f "${module_json}" ]]; then Error "${mod_str} outputs were not found: \n\t\t\tRun -${mod_str}"; exit 1; fi
  status=$(grep "Status" "${module_json}" | awk -F '"' '{print $4}')
  if [ "$status" != "COMPLETED" ]; then
    Error "${mod_str} output has an $status status, try re-running -${mod_str}"; exit 1
  fi
}

function tck_json() {
  qform=($(fslhd "$dwi_b0" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$dwi_b0" | grep sto_ | awk -F "\t" '{print $2}'))
  Info "Creating tractography json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${8}\",
    \"fileInfo\": {
      \"Name\": \"${fod_wmN}\",
      \"qform\": [
          \"${qform[@]:0:4} \",
          \"${qform[@]:4:4} \",
          \"${qform[@]:8:4} \",
          \"${qform[@]:12:8}\"
      ],
      \"sform\": [
          \"${sform[@]:0:4} \",
          \"${sform[@]:4:4} \",
          \"${sform[@]:8:4} \",
          \"${sform[@]:12:8}\"
      ]
    },
    \"Tractography\": {
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
      \"TerminationCriteria\": [\"ACT\"],
      \"weighted_SC\": \"${weighted_SC}\",
      \"TractographySaved\": \"${nocleanup}\"
    }
  }" > "${tckjson}"
}

function json_nativepro_T1w() {
  qform=($(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}'))
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "$1" -transform))
  if [[ "${UNI}" == "FALSE" ]]; then MF="NONE"; fi
  if [[ "${maskbet}" == "TRUE" ]]; then BrainMask="bet"; else BrainMask="mri_synthstrip"; fi
  Info "Creating T1w_nativepro json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${1}\",
    \"fileInfo\": {
        \"VoxelSize\": \"${res}\",
        \"Dimensions\": \"${Size}\",
        \"Strides\": \"${Strides}\",
        \"Offset\": \"${Offset}\",
        \"Multiplier\": \"${Multiplier}\",
        \"Transform\": [
          \"${Transform[@]:0:4} \",
          \"${Transform[@]:4:4} \",
          \"${Transform[@]:8:4} \",
          \"${Transform[@]:12:8}\"
      ],
        \"qform\": [
          \"${qform[@]:0:4} \",
          \"${qform[@]:4:4} \",
          \"${qform[@]:8:4} \",
          \"${qform[@]:12:8}\"
        ],
        \"sform\": [
          \"${sform[@]:0:4} \",
          \"${sform[@]:4:4} \",
          \"${sform[@]:8:4} \",
          \"${sform[@]:12:8}\"
        ]
      },
    \"inputsRawdata\": \"${3}\",
    \"anatPreproc\": [
      {
        \"Resample\": \"LPI\",
        \"Reorient\": \"fslreorient2std\",
        \"NumberOfT1w\": \"$2\",
        \"UNI-T1map\": \"${UNI}\",
        \"UNI-T1map-mf\": \"${MF}\",
        \"BiasFieldCorrection\": \"ANTS N4BiasFieldCorrection\",
        \"WMweightedN4B\": \"${N4wm}\",
        \"N4wmProcessed\": \"${N4wmStatus}\",
        \"RescaleRange\": \"0:100\",
        \"BrainMask\": \"${BrainMask}\"
      }
    ]
  }" > "$4"
}

function json_surf() {
  Info "Creating proc_surf json file"
  qform=($(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}'))
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "${1}" -transform))
  if [[ "$surfdir" == "FALSE" ]]; then
      echo -e "{
        \"micapipeVersion\": \"${Version}\",
        \"LastRun\": \"$(date)\",
        \"fileName\": \"${1}\",
        \"fileInfo\": {
            \"Cropped\": \"${crop}\",
            \"VoxelSize\": \"${res}\",
            \"Dimensions\": \"${Size}\",
            \"Strides\": \"${Strides}\",
            \"Offset\": \"${Offset}\",
            \"Multiplier\": \"${Multiplier}\",
            \"Transform\": [
              \"${Transform[@]:0:4} \",
              \"${Transform[@]:4:4} \",
              \"${Transform[@]:8:4} \",
              \"${Transform[@]:12:8}\"
      ],
            \"qform\": [
              \"${qform[@]:0:4} \",
              \"${qform[@]:4:4} \",
              \"${qform[@]:8:4} \",
              \"${qform[@]:12:8}\"
            ],
            \"sform\": [
              \"${sform[@]:0:4} \",
              \"${sform[@]:4:4} \",
              \"${sform[@]:8:4} \",
              \"${sform[@]:12:8}\"
            ]
          },
        \"SurfaceDir\": \"${2}\",
        \"SurfRecon\": \"${3}\"
      }" > "$4"
  elif [[ "$surfdir" != "FALSE" ]]; then
    echo -e "{
      \"micapipeVersion\": \"${Version}\",
      \"LastRun\": \"$(date)\",
      \"originalDir\": \"${surfdir}\",
      \"SurfaceDir\": \"${2}\",
      \"SurfRecon\": \"${3}\"
    }" > "$4"
  fi
}

function proc_struct_transformations() {
  Info "Creating transformations file: MNI152 >><< T1nativepro"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"Module\": \"proc_structural\",
    \"LastRun\": \"$(date)\",
    \"T1nativepro\": \"${1}\",
    \"MNI152_0p8mm\": \"${2}\",
    \"MNI152_2mm\": \"${3}\",
    \"from-t1nativepro_to-MNI152_0p8mm\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"$1\",
        \"reference\": \"$2\",
        \"transformations\": \"-t ${T1_MNI152_Warp} -t ${T1_MNI152_affine}\",
        \"output\": \"-o from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ],
    \"from-MNI152_0p8mm_to-t1nativepro\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"$2\",
        \"reference\": \"$1\",
        \"transformations\": \"-t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp}\",
        \"output\": \"-o from-MNI152_0.8mm_to-nativepro_mode-image_desc-SyN.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ],
    \"from-t1nativepro_to-MNI152_2mm\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"${1}\",
        \"reference\": \"${3}\",
        \"transformations\": \"-t ${T1_MNI152_Warp/0.8mm/2mm} -t ${T1_MNI152_affine/0.8mm/2mm}\",
        \"output\": \"-o from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ],
    \"from-MNI152_2mm_to-t1nativepro\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"${3}\",
        \"reference\": \"${1}\",
        \"transformations\": \"-t [${T1_MNI152_affine/0.8mm/2mm},1] -t ${T1_MNI152_InvWarp/0.8mm/2mm}\",
        \"output\": \"-o from-MNI152_2mm_to-nativepro_mode-image_desc-SyN.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ]
  }" > "$4"
}

function post_struct_transformations() {
  Info "Creating transformations file: Surface Native >><< T1nativepro"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"Module\": \"post_structural\",
    \"LastRun\": \"$(date)\",
    \"T1nativepro\": \"${1}\",
    \"T1surf\": \"${T1surf}\",
    \"from-t1nativepro_to-fsnative\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"$1\",
        \"reference\": \"$2\",
        \"transformations\": \"[${T1_fsnative_affine},1]\",
        \"output\": \"-o from-nativepro_to-fsnative_mode-image_desc-affine.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ],
    \"from-fsnative_to-t1nativepro\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"$2\",
        \"reference\": \"$1\",
        \"transformations\": \"-t ${T1_fsnative_affine}\",
        \"output\": \"-o from-fsnative_to-nativepro_mode-image_desc-affine.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ]
  }" > "$3"
}

function proc_func_transformations() {
  if [[ ${regAffine}  == "FALSE" ]]; then Mode="SyN"; else Mode="affine"; fi
  Info "Creating transformations file: func space <<>> T1nativepro"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"Module\": \"proc_func\",
    \"LastRun\": \"$(date)\",
    \"Only affine\": \"${regAffine}\",
    \"T1nativepro brain\": \"${T1nativepro_brain}\",
    \"func brain\": \"${fmri_brain}\",
    \"from-t1nativepro_to-func\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"${T1nativepro_brain}\",
        \"reference\": \"${fmri_brain}\",
        \"transformations\": \"$(echo ${2} | sed 's/:/ /g')\",
        \"output\": \"-o from-nativepro_to-${tagMRI}_mode-image_desc-${Mode}.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ],
    \"from-func_to-t1nativepro\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"${fmri_brain}\",
        \"reference\": \"${t1bold}\",
        \"transformations\": \"$(echo ${3} | sed 's/:/ /g')\",
        \"output\": \"-o from-${tagMRI}_to-nativepro_mode-image_desc-${Mode}.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ]
  }" > ${1}
}

function proc_dwi_transformations() {
  if [[ ${regAffine}  == "FALSE" ]]; then Mode="SyN"; else Mode="affine"; fi
  Info "Creating transformations file: DWI space <<>> T1nativepro"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"Module\": \"proc_dwi\",
    \"LastRun\": \"$(date)\",
    \"transform\": \"${Mode}\",
    \"T1nativepro\": \"${T1nativepro}\",
    \"DWI b0\": \"${dwi_b0}\",
    \"from-t1nativepro_to-dwi\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"${T1nativepro}\",
        \"reference\": \"${fod}\",
        \"transformations\": \"$(echo ${2} | sed 's/:/ /g')\",
        \"output\": \"-o from-nativepro_to-dwi${dwi_str_}_mode-image_desc-${Mode}.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ],
    \"from-dwi_to-t1nativepro\": [
      {
        \"Command\": \"antsApplyTransforms\",
        \"input\": \"${dwi_b0}\",
        \"reference\": \"${T1nativepro_brain}\",
        \"transformations\": \"$(echo ${3} | sed 's/:/ /g')\",
        \"output\": \"-o from-dwi${dwi_str_}_to-nativepro_mode-image_desc-${Mode}.nii.gz\",
        \"options\": \"-d 3 -v -u int\"
      }
    ]
  }" > ${1}
}

function slim_proc_struct(){
  Info "Erasing temporary files"
  Do_cmd rm -rf ${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_to_std_sub*
  Do_cmd rm -rf ${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_pveseg.nii.gz
  Do_cmd rm -rf ${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_mixeltype.nii.gz
  Do_cmd rm -rf ${proc_struct}/${idBIDS}_space-nativepro_T1w_brain_seg.nii.gz
  Do_cmd rm -rf ${proc_struct}/first
  Do_cmd rm -rf ${dir_warp}/*Warped.nii.gz 2>/dev/null
}

function json_nativepro_mask() {
  qform=($(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}'))
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "$1" -transform))
  Info "Creating T1nativepro_brain json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${1}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"Transform\": [
        \"${Transform[@]:0:4} \",
        \"${Transform[@]:4:4} \",
        \"${Transform[@]:8:4} \",
        \"${Transform[@]:12:8}\"
      ],
    \"inputNIFTI\": {
      \"Name\": \"${T1nativepro}\",
      \"qform\": [
        \"${qform[@]:0:4} \",
        \"${qform[@]:4:4} \",
        \"${qform[@]:8:4} \",
        \"${qform[@]:12:8}\"
      ],
      \"sform\": [
        \"${sform[@]:0:4} \",
        \"${sform[@]:4:4} \",
        \"${sform[@]:8:4} \",
        \"${sform[@]:12:8}\"
      ]
    },
    \"BinaryMask\": [
      {
        \"BinaryMask_original\": \"$2\",
        \"BinaryMask_from\": \"MNI152_0.8mm\",
        \"BinaryMask_antsApplyTransforms\": \"$4\"
      }
    ]
  }" > "$3"
}

function json_nativepro_flair() {
  qform=($(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}'))
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "$1" -transform))
  Info "Creating Flair json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"flairScanStr\": \"${flairScanStr}\",
    \"fileName\": \"${1}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"regSynth\": \"${synth_reg}\",
    \"mode_wm\": \"${mode_wm}\",
    \"mode_gm\": \"${mode_gm}\",
    \"mode_brain\": \"${mode_brain}\",
    \"TransformCmd\": {
        \"BinaryMask_antsApplyTransforms\": \"$2\"
      },
    \"Transform\": [
        \"${Transform[@]:0:4} \",
        \"${Transform[@]:4:4} \",
        \"${Transform[@]:8:4} \",
        \"${Transform[@]:12:8}\"
      ],
    \"inputNIFTI\": {
      \"Name\": \"$flairScan\",
      \"qform\": [
        \"${qform[@]:0:4} \",
        \"${qform[@]:4:4} \",
        \"${qform[@]:8:4} \",
        \"${qform[@]:12:8}\"
      ],
      \"sform\": [
        \"${sform[@]:0:4} \",
        \"${sform[@]:4:4} \",
        \"${sform[@]:8:4} \",
        \"${sform[@]:12:8}\"
      ]
    }
  }" > "$3"
}

function json_nativepro_qt1() {
  qform=($(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}'))
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "$1" -transform))
  Info "Creating T1nativepro_qt1 json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"fileName\": \"${1}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"TransformCmd\": {
        \"antsApplyTransforms\": \"$2\"
      },
    \"Transform\": [
        \"${Transform[@]:0:4} \",
        \"${Transform[@]:4:4} \",
        \"${Transform[@]:8:4} \",
        \"${Transform[@]:12:8}\"
      ],
    \"inputNIFTI\": {
      \"Name\": \"$bids_T1map\",
      \"qform\": [
        \"${qform[@]:0:4} \",
        \"${qform[@]:4:4} \",
        \"${qform[@]:8:4} \",
        \"${qform[@]:12:8}\"
      ],
      \"sform\": [
        \"${sform[@]:0:4} \",
        \"${sform[@]:4:4} \",
        \"${sform[@]:8:4} \",
        \"${sform[@]:12:8}\"
      ]
    }
  }" > "$3"
}

function json_poststruct() {
  Info "Creating post_structural json file"
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "$1" -transform))
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"SurfRecon\": \"${recon}\",
    \"SurfaceDir\": \"${2}\",
    \"Atlas\": [
        \"${atlas}\"
      ],
    \"NativeSurfSpace\": {
        \"fileName\": \"${1}\",
        \"VoxelSize\": \"${res}\",
        \"Dimensions\": \"${Size}\",
        \"Strides\": \"${Strides}\",
        \"Offset\": \"${Offset}\",
        \"Multiplier\": \"${Multiplier}\",
        \"Transform\": [
          \"${Transform[@]:0:4} \",
          \"${Transform[@]:4:4} \",
          \"${Transform[@]:8:4} \",
          \"${Transform[@]:12:8}\"
          ]
      }
  }" > "$3"
}

function json_func() {
  qform=($(fslhd "$func_processed" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$func_processed" | grep sto_ | awk -F "\t" '{print $2}'))
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"Tag\": \"${tagMRI}\",
    \"Acquisition\": \"${acq}\",
    \"Name\": \"${func_processed}\",
    \"sform\": [
        \"${sform[@]:0:4} \",
        \"${sform[@]:4:4} \",
        \"${sform[@]:8:4} \",
        \"${sform[@]:12:8}\"
      ],
    \"qform\": [
        \"${qform[@]:0:4} \",
        \"${qform[@]:4:4} \",
        \"${qform[@]:8:4} \",
        \"${qform[@]:12:8}\"
      ],
    \"Preprocess\": {
        \"MainScan\": \"${mainScan_orig}\",
        \"Resample\": \"LPI\",
        \"Reorient\": \"fslreorient2std\",
        \"MotionCorrection\": \"3dvolreg AFNI $(afni -version | awk -F ':' '{print $2}')\",
        \"MotionCorrection\": [\"${func_volum}/${idBIDS}_space-func_spikeRegressors_FD.1D\"],
        \"MainPhaseScan\": \"${mainPhaseScan_orig}\",
        \"ReversePhaseScan\": \"${reversePhaseScan_orig}\",
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
        \"dropTR\": \"${dropTR}\"
      }
  }" > "$1"
}

function json_mpc() {
  qform=($(fslhd "$1" | grep qto_ | awk -F "\t" '{print $2}'))
  sform=($(fslhd "$1" | grep sto_ | awk -F "\t" '{print $2}'))
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "${1}" -transform))
  Info "Creating MPC json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"Module\": \"Microstructural profile covariance\",
    \"acquisition\": \"${mpc_str}\",
    \"microstructural_img\": \"${1}\",
    \"microstructural_reg\": \"${regImage}\",
    \"reference_mri\": \"${qMRI_reference}\",
    \"warped_qmri\": \"${qMRI_warped}\",
    \"regSynth\": \"${synth_reg}\",
    \"reg_nonlinear\": \"${reg_nonlinear}\",
    \"registration\": \"${reg}\",
    \"num_surfs\": \"${num_surfs}\",
    \"VoxelSize\": \"${res}\",
    \"Dimensions\": \"${Size}\",
    \"Strides\": \"${Strides}\",
    \"Offset\": \"${Offset}\",
    \"Multiplier\": \"${Multiplier}\",
    \"Transform\": [
        \"${Transform[@]:0:4} \",
        \"${Transform[@]:4:4} \",
        \"${Transform[@]:8:4} \",
        \"${Transform[@]:12:8}\"
      ],
    \"qform\": [
        \"${qform[@]:0:4} \",
        \"${qform[@]:4:4} \",
        \"${qform[@]:8:4} \",
        \"${qform[@]:12:8}\"
      ],
    \"sform\": [
        \"${sform[@]:0:4} \",
        \"${sform[@]:4:4} \",
        \"${sform[@]:8:4} \",
        \"${sform[@]:12:8}\"
      ]
  }" > "$2"
}

function json_dwipreproc() {
  res=$(mrinfo "$1" -spacing)
  Size=$(mrinfo "$1" -size)
  Strides=$(mrinfo "$1" -strides)
  Offset=$(mrinfo "$1" -offset)
  Multiplier=$(mrinfo "$1" -multiplier)
  Transform=($(mrinfo "${1}" -transform))

  res_rpe=$(mrinfo "$4" -spacing)
  Size_rpe=$(mrinfo "$4" -size)
  Strides_rpe=$(mrinfo "$4" -strides)
  Offset_rpe=$(mrinfo "$4" -offset)
  Multiplier_rpe=$(mrinfo "$4" -multiplier)
  Transform_rpe=($(mrinfo "$4" -transform))

  Info "Creating DWI preproc json file"
  echo -e "{
    \"micapipeVersion\": \"${Version}\",
    \"LastRun\": \"$(date)\",
    \"Class\": \"DWI preprocessing\",
    \"rpe_all\": \"${rpe_all}\",
    \"dwi_acq\": \"${dwi_acq}\",
    \"Only Affine\": \"${regAffine}\",
    \"B0 threshold\": \"${b0thr}\",
    \"Bvalue scaling\": \"${bvalscale}\",
    \"regSynth\": \"${synth_reg}\",
    \"dwi_upsample\": \"${dwi_upsample}\",
    \"DWIpe\": {
        \"fileName\": \"${bids_dwis[*]}\",
        \"NumberOfInputs\": \"${#bids_dwis[*]}\",
        \"VoxelSizepe\": \"${res}\",
        \"Dimensionspe\": \"${Size}\",
        \"Strides\": \"${Strides}\",
        \"Offset\": \"${Offset}\",
        \"Multiplier\": \"${Multiplier}\",
        \"Transform\": [
          \"${Transform[@]:0:4} \",
          \"${Transform[@]:4:4} \",
          \"${Transform[@]:8:4} \",
          \"${Transform[@]:12:8}\"
      ]
    },
    \"DWIrpe\": {
        \"fileName\": \"${dwi_reverse[*]}\",
        \"NumberOfInputs\": \"${#dwi_reverse[*]}\",
        \"VoxelSizepe\": \"${res_rpe}\",
        \"Dimensionspe\": \"${Size_rpe}\",
        \"Strides\": \"${Strides_rpe}\",
        \"Offset\": \"${Offset_rpe}\",
        \"Multiplier\": \"${Multiplier_rpe}\",
        \"Transform\": [
          \"${Transform_rpe[@]:0:4} \",
          \"${Transform_rpe[@]:4:4} \",
          \"${Transform_rpe[@]:8:4} \",
          \"${Transform_rpe[@]:12:8}\"
      ]
    },
    \"Denoising\": \"Marchenko-Pastur PCA denoising, dwidenoise\",
    \"GibbsRingCorrection\": \"mrdegibbs\",
    \"dwiflspreproc\": {
        \"input\": \"${dwi_4proc}\",
        \"output\": \"${dwi_corr}\",
        \"Shells\": \"${shells[*]}\",
        \"pe_dir\": \"${pe_dir}\",
        \"ReadoutTime\": \"${ReadoutTime}\",
        \"Options\": \"${opt}\",
        \"slm\": \"linear\"
    },
    \"B1fieldCorrection\": \"ANTS N4BiasFieldCorrection\",
    \"DWIprocessed\": \"$2\"
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
      rm -Rf "$tmp" 2>/dev/null
  else
      echo -e "micapipe tmp directory was not erased: \n\t\t${tmp}";
  fi
  cd "$here"
  bids_variables_unset
  if [[ ! -z "$OLD_PATH" ]]; then  export PATH=$OLD_PATH; unset OLD_PATH; fi
}

function missing_arg() {
  arg=($id $out $BIDS)
  if [ ${#arg[@]} -lt 3 ]; then
  Error "One or more mandatory arguments are missing:
                 -sub  : $id
                 -out  : $out
                 -bids : $BIDS
          -h | -help (print help)"
  exit 1; fi
}

function inputs_realpath() {
  # Get the real path of the Inputs
  out=$(realpath $out)/micapipe_v0.2.0
  BIDS=$(realpath $BIDS)
  id=${id/sub-/}
  here=$(pwd)
}

function map_to-surfaces(){
  # -----------------------------------------------
  # Volume to surface mapping
  # Function that maps a MRI volume to a surfaces
  # Using work bench commands and multiple surfaces:
  # fsnative, fsaverage5, fsLR-5k and fsLR-32k
  # -----------------------------------------------
  # Input variables
  mri_map=$1                      # MRI map from where data will be mapped
  surf_fs=$2                      # Surface to map the MRI (MUST be surf-fsnative, same space ast mri_map)
  map_on_surf=${3}                # Outname of the data mapped on the surface
  H=$4                            # Hemisphere {L, R}
  label_data=$5                   # label of the map (e.g. FA, ADC, flair, T2star, MTR)
  out_map=$6
  # Map to highest resolution surface (fsnative: more vertices)
  wb_command -volume-to-surface-mapping "${mri_map}" "${surf_fs}" "${map_on_surf}" -trilinear
  # Map from volume to surface for each surfaces
  for Surf in "fsLR-32k" "fsaverage5" "fsLR-5k"; do
    surf_id=${idBIDS}_hemi-${H}_surf
    wb_command -metric-resample "${map_on_surf}" \
        "${dir_conte69}/${surf_id}-fsnative_label-sphere.surf.gii" \
        "${util_surface}/${Surf}.${H}.sphere.reg.surf.gii" \
        BARYCENTRIC "${out_map}/${surf_id}-${Surf}_label-${label_data}.func.gii"
  done
}

function steps() {
  Note "N     :" "${N}"
  Note "Nsteps:" "${Nsteps}"
}