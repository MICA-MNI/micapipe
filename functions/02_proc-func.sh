#!/bin/bash
# Resting state preprocessing
# Written by Casey Paquola and Reinder Vos De Wael (Oct 2018).
# and a tiny bit from Sara (Feb 2019)...
# and a whole lot from Sara (August 2019)
# and incorporation to mica-pipe by Raul (August-September 2020)
# and addition of a bunch of fancy flags by Jessica (October-November 2020)
# and updated mapping by Jordan (March 2023)
#
# Resting state fMRI processing with bash:
#
# Preprocessing workflow for func.
#
# This workflow makes use of AFNI, FSL, ANTs, FIX
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out Directory
#
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
changeTopupConfig=$8
changeIcaFixTraining=$9
thisMainScan=${10}
thisPhase=${11}
mainScanStr=${12}
func_pe=${13}
func_rpe=${14}
performNSR=${15}
performGSR=${16}
noFIX=${17}
sesAnat=${18}
dropTR=${19}
noFC=${20}
PROC=${21}
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE"/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check dependencies Status: POST_STRUCTURAL
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfRecon "${post_struct_json}" | awk -F '"' '{print $4}')
set_surface_directory "${recon}"

if [[ "$sesAnat" != FALSE  ]]; then
  sesAnat=${sesAnat/ses-/}
  BIDSanat="${subject}_ses-${sesAnat}"
  dir_anat="${out}/${subject}/ses-${sesAnat}/anat"

  dir_conte69="${dir_anat}/surf/conte69"
  T1nativepro="${dir_anat}/${BIDSanat}_space-nativepro_T1w.nii.gz"
  T1nativepro_brain="${dir_anat}/${BIDSanat}_space-nativepro_T1w_brain.nii.gz"
  T1nativepro_mask="${dir_anat}/${BIDSanat}_space-nativepro_T1w_brain_mask.nii.gz"
  dir_subjsurf="${dir_surf}/${subject}_ses-${sesAnat}"
  T1surf="${dir_subjsurf}/mri/T1.mgz"
else
  BIDSanat="${idBIDS}"
  dir_anat="${proc_struct}"
fi
subject_dir
dir_volum="${subject_dir}/parc"
T1_seg_subcortex="${dir_volum}/${BIDSanat}_space-nativepro_T1w_atlas-subcortical.nii.gz"
T1_seg_cerebellum="${dir_volum}/${BIDSanat}_space-nativepro_T1w_atlas-cerebellum.nii.gz"

### CHECK INPUTS: func, phase encoding, structural proc, topup and ICA-FIX files
Info "Inputs:"
Note "Topup Config     :" "$changeTopupConfig"
Note "ICA fix training :" "$changeIcaFixTraining"
if [[ "$mainScanStr" == DEFAULT ]]; then Note "Main scan        :" "$thisMainScan"; else
Note "Main scan        :" "$mainScanStr"; fi
Note "Phase scan       :" "$func_pe"
Note "Reverse Phase    :" "$func_rpe"
Note "Perform NSR      :" "$performNSR"
Note "Perform GSR      :" "$performGSR"
Note "No FIX           :" "$noFIX"
Note "Longitudinal ses :" "$sesAnat"
Note "Drop TR          :" "${dropTR}"
Note "Surface          :" "${recon}"

#------------------------------------------------------------------------------#
if [[ "$mainScanStr" == DEFAULT ]]; then
    # Main scan
    N_mainScan=${#bids_mainScan[@]}
    if [ "$N_mainScan" -gt 1 ]; then
        if [[ "${thisMainScan}" == "DEFAULT" ]]; then
            Error "Multiple func runs found in BIDS rawdata directory! Please specify which run should be processed using flag -mainScanRun"; exit;
        elif [ "$thisMainScan" -gt "$N_mainScan" ]; then
            Warning "Specified run number (${thisMainScan}) is greater than number of func scans scans found ($N_mainScan). Using first filename in list as default";
            mainScan=${bids_mainScan[0]}
        else
            Info "Found $N_mainScan func scans, processing specified scan # $thisMainScan"
            mainScan=${bids_mainScan[$thisMainScan-1]}
        fi
    else
        mainScan=${bids_mainScan[0]}
        if [[ "$thisMainScan" == "DEFAULT" ]]; then
            Info "No run number specified for func scan and did not find more than one run for main scan - all good!"
        else
            if [ "$thisMainScan" -gt "$N_mainScan" ]; then
                Warning "Found one or less func scan, but specified run number = $thisMainScan). Using first filename in list as default";
            fi
        fi
    fi

    # Main scan json
    N_mainScanJson=${#bids_mainScanJson[@]}
    if [ "$N_mainScanJson" -gt 1 ]; then
        if [[ "${thisMainScan}" == "DEFAULT" ]]; then
            Error "Found multiple .json files for main func scan in BIDS rawdata directory! Please specify which run should be processed using flag -mainScanRun"; exit;
        elif [ "$thisMainScan" -gt "$N_mainScanJson" ]; then
            Warning "Specified run number (${thisMainScan}) is greater than number of func json files found for main scan ($N_mainScan). Using first filename in list as default";
            mainScanJson=${bids_mainScan[0]}
        else
            Info "Found ${N_mainScanJson} func scan json files, using specified run # ${thisMainScan}"
            mainScanJson=${bids_mainScanJson[$thisMainScan-1]}
        fi
    else
        Info "Using default json scan: ${bids_mainScanJson[0]}"
        mainScanJson=${bids_mainScanJson[0]}
    fi
else
    Info "Using user provided main scan string(s): ${mainScanStr}"
    IFS=',' read -ra func_main <<< $mainScanStr
    func_json=(${func_main[*]})
    for i in "${!func_main[@]}"; do
        func_main[i]=$(ls "${subject_bids}/func/${idBIDS}_${func_main[$i]}".nii* 2>/dev/null);
        func_json[i]=$(ls "${subject_bids}/func/${idBIDS}_${func_json[$i]}".json 2>/dev/null)
    done # Full path
    mainScan=(${func_main[*]})
    mainScanJson=(${func_json[*]})
fi
# If no json is found search at the top BIDS directory
if [[ -z ${mainScanJson} ]]; then
    IFS=',' read -ra func_json <<< $mainScanStr
    for i in "${!func_json[@]}"; do
        func_json[i]=$(ls "${BIDS}/${func_json[$i]}.json" 2>/dev/null)
    done # Full path
    mainScanJson=(${func_json[*]})
fi

#------------------------------------------------------------------------------#
# Phase encoding
N_mainPhase=${#bids_mainPhase[@]}
N_revPhase=${#bids_reversePhase[@]}
if [ "$N_mainPhase" -gt 1 ] || [ "$N_revPhase" -gt 1 ]; then
    if [[ "$thisPhase" == "DEFAULT" ]]; then
        Error "Found multiple phase reversal runs in BIDS rawdata directory! Please specify which run should be processed using flag -phaseReversalRun:\n ${bids_reversePhase[*]}"; exit;
    elif [ "$thisPhase" -gt "$N_mainPhase" ] || [ "$thisPhase" -gt "$N_revPhase" ]; then
        Warning "Specified run number ($thisPhase) is greater than number of phase reversal scans scans found ($N_mainPhase and $N_revPhase). Using first filename in list as default";
        mainPhaseScan=${bids_mainPhase[$thisPhase-1]}
        reversePhaseScan=${bids_reversePhase[$thisPhase-1]}
    else
        Info "Found $N_mainPhase and $N_revPhase phase reversal scans, processing specified scan # $thisPhase"
        mainPhaseScan=${bids_mainPhase[$thisPhase-1]}
        reversePhaseScan=${bids_reversePhase[$thisPhase-1]}
    fi
else
    mainPhaseScan=${bids_mainPhase[0]}
    reversePhaseScan=${bids_reversePhase[0]}
    if [[ "$thisPhase" == "DEFAULT" ]]; then
        Info "No run number specified for phase reversals and did not find more than one phase reversal scan - all good!"
    else
        if [ "$thisPhase" -gt "$N_mainPhase" ] || [ "$thisPhase" -gt "$N_revPhase" ]; then
            Warning "Specified run number ($thisPhase) is greater than number of phase reversal scans scans found ($N_mainPhase and $N_revPhase). Using first filename in list as default"; fi
    fi
fi

# Manually defined Phase scan and reverse phase scan
if [[ "$func_pe" != DEFAULT ]] && [[ -f "$func_pe" ]]; then mainPhaseScan="$func_pe"; fi
if [[ "$func_rpe" != DEFAULT ]] && [[ -f "$func_rpe" ]]; then reversePhaseScan="$func_rpe"; fi

# Check inputs
for i in ${mainScan[*]}; do
    if [ ! -f "$i" ]; then Error "Couldn't find $id main func scan : \n\t ls ${i}"; exit; fi
done
for i in ${mainScanJson[*]}; do
    if [ ! -f "$i" ]; then Error "Couldn't find $id main func scan json file: \n\t ls ${i}"; exit; fi
done
if [ -z "$mainPhaseScan" ]; then  Warning "Subject $id doesn't have a Main Phase Scan (pe): TOPUP will run only if a rpe is provided"; fi
if [ -z "$reversePhaseScan" ]; then Warning "Subject $id doesn't have Reverse Phase Scan (rpe): TOPUP will be skipped"; fi

# Check topup input
if [[ ${changeTopupConfig} == "DEFAULT" ]]; then
    Info "Will use default config file for TOPUP: ${topupConfigFile}"
else
    topupConfigFile=${changeTopupConfig}
    Info "Will use specified config file for TOPUP: ${topupConfigFile}"
fi

# Check FIX: run or no?
if [[ "$noFIX" -eq 1 ]]; then
    Info "ICA-FIX will be skipped! Consider performing nuisance signal regression with <-NSR> or <-GSR>"

    # Check ICA-FIX Training input
    if [[ ! ${changeIcaFixTraining} == "DEFAULT" ]]; then
        Error "If ICA-FIX is skipped, <-icafixTraining> must remain empty"; exit; fi
else
    Info "ICA-FIX pipeline will be run!"

    # Check ICA-FIX Training input
    if [[ ${changeIcaFixTraining} == "DEFAULT" ]]; then
        Info "Will use default training file for ICA-FIX: ${icafixTraining}"
    else
        icafixTraining=${changeIcaFixTraining}
        Info "Will use specified training file for ICA-FIX: ${icafixTraining}"
    fi
fi

# Check nuisance signal regression
if [[ $performNSR == 1 ]]; then
    Info "White matter and CSF signals will be regressed from processed timeseries"
elif [[ $performGSR == 1 ]]; then
    Info "Global, white matter and CSF signals will be regressed from processed timeseries"
else
    Info "Global, white matter and CSF signal regression will not be performed (default)"
fi
# gettin dat from mainScanJson exit if Not found
unset readoutTime RepetitionTime EchoNumber EchoTime
for json in ${mainScanJson[*]}; do
    readoutTime+=($(grep TotalReadoutTime "${json}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}'))
    RepetitionTime+=($(grep RepetitionTime "${json}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}'))
    EchoNumber+=($(grep EchoNumber "${json}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}'))
    EchoTime+=($(grep EchoTime "${json}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}'))
done
if [[ -z "$readoutTime" ]]; then Warning "readoutTime is missing in $mainScanJson, if TOPUP was selected it will likely FAIL"; fi
if [[ -z "$RepetitionTime" ]]; then Error "RepetitionTime is missing in $mainScanJson $RepetitionTime"; exit; fi

#------------------------------------------------------------------------------#
# If mainScan is an array with more than one file we'll assume it's multiecho
if [ ${#mainScan[@]} -eq 1 ]; then
  acq="se"
elif [ ${#mainScan[@]} -gt 1 ]; then
  acq="me"; dropTR="FALSE"; noFIX=1
fi

# func directories
fmri_tag=$(echo "${mainScan[0]/echo-1_/}" | awk -F "${idBIDS}_" '{print $2}' | cut -d'.' -f1); fmri_tag="desc-${acq}_${fmri_tag}"
tagMRI="${fmri_tag/desc-/}"
proc_func="$subject_dir/func/${fmri_tag}"

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_func-${fmri_tag}.json"
micapipe_check_json_status "${module_json}" "proc_func"

#------------------------------------------------------------------------------#
Title "functional MRI processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
bids_print.variables-func
Note "Saving temporal dir:" "$nocleanup"
Note "Parallel processing      :" "${threads} threads"
Note "proc_fun outputs:" "${proc_func}"
Note "tagMRI:" "${tagMRI}"

# Save original files for the json file
export mainScan_orig=${mainScan[0]}
export mainPhaseScan_orig=${mainPhaseScan}
export reversePhaseScan_orig=${reversePhaseScan}

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0
# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_proc-func_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Define directories
export SUBJECTS_DIR="$dir_surf"

func_volum="${proc_func}/volumetric"   # volumetricOutputDirectory
func_surf="${proc_func}/surf"      # surfaceOutputDirectory
func_ICA="${tmp}/ICA_MELODIC"      # ICAOutputDirectory

# Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script).
for x in "$func_surf" "$func_volum"; do
    [[ ! -d "${x}" ]] && mkdir -p "${x}"
done

#------------------------------------------------------------------------------#
# Scans to process
toProcess=($reversePhaseScan $mainPhaseScan)
tags=("reversePhaseScan" "mainPhaseScan")
func_lab="_space-func_desc-${acq}"
func_nii="${func_volum}/${idBIDS}${func_lab}".nii.gz

# Processing functions
function func_reoMC() {
  # Function that reorients the func file and motion corrects
  # Positional arguments: 1) nifti, 2) tag,
  # Get basic parameters
  local rawNifti=$1
  local tag=$2
  local EchoN=$3
  # IF FILE IS NOT FOUND DON'T RUN
  if [[ ! -z "${rawNifti}" ]] && [[ -f "${rawNifti}" ]]; then
        Info "Processing $tag scan"
        Note "func file:" "$rawNifti"

        # Drop first five TRs and reorient to standard
        if [ "$tag" == "mainScan" ] && [ "$dropTR" == "TRUE" ]; then
            Do_cmd nifti_tool -cbl -prefix "${tmp}/${tag}_trDrop.nii.gz" -infiles "$rawNifti"'[5..$]'
            Do_cmd 3dresample -orient LPI -prefix "${tmp}/${tag}_reo.nii.gz" -inset "${tmp}/${tag}_trDrop.nii.gz"
            Do_cmd fslreorient2std "${tmp}/${tag}_reo.nii.gz" "${tmp}/${tag}_reo.nii.gz"
        else
            Do_cmd 3dresample -orient LPI -prefix "${tmp}/${tag}_reo.nii.gz" -inset "$rawNifti"
            Do_cmd fslreorient2std "${tmp}/${tag}_reo.nii.gz" "${tmp}/${tag}_reo.nii.gz"
        fi

        # Motion correction within scans <<<  Only first echo for tedana
        if [[ ${EchoN} -eq 1 ]]; then
            Do_cmd fslmaths "${tmp}/${tag}_reo.nii.gz" -Tmean "${tmp}/${tag}_reoMean.nii.gz"
            Do_cmd 3dvolreg -Fourier -twopass -base "${tmp}/${tag}_reoMean.nii.gz" \
                            -zpad 4 -prefix "${tmp}/${tag}_mc.nii.gz" \
                            -1Dfile "${func_volum}/${idBIDS}${func_lab}_${tag}.1D" \
                            -1Dmatrix_save "${tmp}/${tag}"_mat_mc.1D \
                            "${tmp}/${tag}_reo.nii.gz"
            Do_cmd fslmaths "${tmp}/${tag}_mc.nii.gz" -Tmean "${tmp}/${tag}_mcMean.nii.gz"
        fi
  fi
}


function func_MCoutliers() {
  # Function that generates the motion outliers file
  # Calculate motion outliers with FSL
  local outfile=$1
  if [[ ! -f "${outfile}" ]]; then ((N++))
      Do_cmd fsl_motion_outliers -i "${tmp}/mainScan_reo.nii.gz" \
                                 -o "${func_volum}/${idBIDS}${func_lab}_spikeRegressors_FD.1D" \
                                 -s "${func_volum}/${idBIDS}${func_lab}_metric_FD.1D" --fd
      Do_cmd mv "${func_volum}/${idBIDS}${func_lab}_mainScan.1D ${outfile}"; ((Nsteps++))
  else
      Info "Subject ${id} has a ${func_lab/_/}.1D with motion outliers"; ((Nsteps++)); ((N++))
  fi
}

function func_topup() {
  # Function that applies the distortion correction to the REO/MC file(s)
  #
  # If ONLY Reverse phase scan is provided mainScan will be the mainPhaseScan
  if [[ -f "$reversePhaseScan" ]] && [[ ! -f "$mainPhaseScan" ]]; then
        main_pe="NO-main_pe"
        Warning "reversePhaseScan was found but NO mainPhaseScan, using mainScan as mainPhaseScan"
        mainPhaseScan="${tmp}/mainPhaseScan_mc.nii.gz"
        Do_cmd cp "${tmp}/mainScan_mc.nii.gz" "$mainPhaseScan"
        Do_cmd cp "${tmp}/mainScan_mcMean.nii.gz" "${tmp}/func_mainPhaseAlignedMean.nii.gz"
  fi

  # Only do distortion correction if reverse phase encoding images were provided,
  # if not then rename the motion corrected mainScan to $func_nii.
  if [ -z "${mainPhaseScan}" ] || [ -z "${reversePhaseScan}" ]; then ((N++))
      Warning "No pe or rpe (AP, PA) acquisition were found, TOPUP will be skip!!!!!!!"
      export statusTopUp="NO"
      Do_cmd mv -v "${tmp}/mainScan_mc.nii.gz" "${func_nii}"; ((Nsteps++))
  else
      if [[ ! -f "${func_volum}/TOPUP.txt" ]] && [[ ! -f "${func_nii}" ]]; then ((N++))
        # NOTE print readout times
          mainPhaseScanMean=$(find "$tmp"    -maxdepth 1 -name "*mainPhaseScan_mcMean.nii.gz")
          mainPhaseScan=$(find "$tmp"        -maxdepth 1 -name "*mainPhaseScan_mc.nii.gz")
          reversePhaseScanMean=$(find "$tmp" -maxdepth 1 -name "*reversePhaseScan_mcMean.nii.gz")
          reversePhaseScan=$(find "$tmp"     -maxdepth 1 -name "*reversePhaseScan_mc.nii.gz")
          mainScan=$(find "$tmp"             -maxdepth 1 -name "*mainScan_mc.nii.gz")

          Do_cmd flirt -in "$reversePhaseScanMean" -ref "${tmp}/mainScan_mcMean.nii.gz" -omat "${tmp}/func_tmpXfmSecondary.omat"
          Do_cmd flirt -in "$reversePhaseScan" -ref "${tmp}/mainScan_mcMean.nii.gz" -applyxfm -init "${tmp}/func_tmpXfmSecondary.omat" -out "${tmp}/func_secondaryPhaseAligned.nii.gz"
          Do_cmd fslmaths "${tmp}/func_secondaryPhaseAligned.nii.gz" -Tmean "${tmp}/func_secondaryPhaseAlignedMean.nii.gz"

          if [[ "$main_pe" != "NO-main_pe" ]]; then
              Do_cmd flirt -in "$mainPhaseScanMean" -ref "${tmp}/mainScan_mcMean.nii.gz" -omat "${tmp}/func_tmpXfmMain.omat"
              Do_cmd flirt -in "$mainPhaseScan" -ref "${tmp}/mainScan_mcMean.nii.gz" -applyxfm -init "${tmp}/func_tmpXfmMain.omat" -out "${tmp}/func_mainPhaseAligned.nii.gz"
              Do_cmd fslmaths "${tmp}/func_mainPhaseAligned.nii.gz" -Tmean "${tmp}/func_mainPhaseAlignedMean.nii.gz"
          fi

          # Distortion correction
          echo -e "0 1 0 ${readoutTime[0]} \n0 -1 0 ${readoutTime[0]}" > "${tmp}/func_topupDataIn.txt"
          Info "topup datain:\n$(cat "${tmp}/func_topupDataIn.txt")"
          Do_cmd fslmerge -t "${tmp}/func_mergeForTopUp.nii.gz" "${tmp}/func_mainPhaseAlignedMean.nii.gz" "${tmp}/func_secondaryPhaseAlignedMean.nii.gz"
          Do_cmd topup --imain="${tmp}/func_mergeForTopUp.nii.gz" --datain="${tmp}/func_topupDataIn.txt" --config="${topupConfigFile}" --out="${tmp}/func_topup"
          Do_cmd applytopup --imain="${mainScan}" --inindex=1 --datain="${tmp}/func_topupDataIn.txt" --topup="${tmp}/func_topup" --method=jac --out="${func_nii}"

          # Check if it worked
          if [[ ! -f "${func_nii}" ]]; then Error "Something went wrong while running TOPUP check ${tmp} and log:\n\t\t${dir_logs}/proc_func_$(date +'%d-%m-%Y').txt"; cleanup "$tmp" "$nocleanup" "$here"; exit; fi
          export statusTopUp="YES"; ((Nsteps++))
      else
          Info "Subject ${id} has a distortion corrected functional MRI (TOPUP)"; export statusTopUp="YES"; ((Nsteps++)); ((N++))
      fi
  fi
}

#------------------------------------------------------------------------------#
# Begining of the REAL processing
status="INCOMPLETE"
GSR=0
# Processing fMRI acquisitions.
if [[ ! -f "${func_volum}/${idBIDS}${func_lab}_preproc".nii.gz ]]; then
    # Reorient and motion correct main(s) fMRI
    for i in "${!mainScan[@]}"; do n=$((i+1))
      func_reoMC "${mainScan[i]}" "mainScan${n/1/}" $n
    done
    # Reorient and motion correct fMRI rpe and pe
    for i in "${!toProcess[@]}"; do
      func_reoMC "${toProcess[i]}" "${tags[i]}" 1
    done

    # Run Tedana
    if [[ ${acq} == "me" ]]; then
        Info "Multiecho fMRI acquisition will be process with tedana"
        GSR=1
        scans4tedana=($(ls "${tmp}"/mainScan?_reo.nii.gz))
        Info "Apply motion correction to echos"
        for me in ${!scans4tedana[@]}; do
          3dAllineate -base "${tmp}"/mainScan_mcMean.nii.gz \
            -input "${scans4tedana[${me}]}" -1Dmatrix_apply "${tmp}"/mainScan_mat_mc.1D -prefix "${scans4tedana[${me}]/_reo/_mc}"
        done
        scans4tedana=("${tmp}"/mainScan_mc.nii.gz $(ls "${tmp}"/mainScan?_mc.nii.gz))
        Note "Files      :" "${scans4tedana[*]/${tmp}/}" # this will print the string full path is in mainScan
        Note "EchoNumber :" "${EchoNumber[*]}"
        Note "EchoTime   :" "${EchoTime[*]}"

        tedana_dir="${tmp}"/tedana

        # prepare files for tedana
        mkdir -p "${tedana_dir}"
        # create a mask for tedana
        tmp_func_mean="${tmp}/mainScan01_mean.nii.gz"
        tmp_t1w_brain_res="${tmp}/nativepro_brain_rescaled.nii.gz"
        tmp_affineStr="${tmp}/from-nativepro_to-mainScan01_"
        tmp_aff_mat="${tmp_affineStr}0GenericAffine.mat"
        tmp_func_mask="${tmp}/mainScan01_mean_mask.nii.gz"

        # Create a loose mask for tedana
        fslmaths "${scans4tedana[0]}" -Tmean "${tmp_func_mean}"
        voxels=$(mrinfo "${tmp_func_mean}" -spacing); voxels="${voxels// /,}"
        Do_cmd flirt -applyisoxfm "${voxels}" -in "${T1nativepro_brain}" -ref "${T1nativepro_brain}" -out "${tmp_t1w_brain_res}"
        # [fixedImage,movingImage,initializationFeature]
        centeralign="[${tmp_func_mean},${tmp_t1w_brain_res},0]"
        Do_cmd antsRegistrationSyNQuick.sh -d 3 -m "$tmp_t1w_brain_res" -f "$tmp_func_mean" -o "$tmp_affineStr" -t a -n "$threads" -p d -i ${centeralign}
        Do_cmd antsApplyTransforms -d 3 -i "$T1nativepro_mask" -r "$tmp_func_mean" -t "$tmp_aff_mat" -o "${tmp_func_mask}" -u int

        # Run tedana
        tedana -d $(printf "%s " "${scans4tedana[@]}") -e $(printf "%s " "${EchoTime[@]}") --out-dir "${tedana_dir}" --mask "${tmp_func_mask}"

        # Overwite the motion corrected to insert this into topup.
        ## TODO: func_topup should take proper input arguments instead of relying on architecture implemented in other functions.
        mainScan=$(find "$tmp" -maxdepth 1 -name "*mainScan_mc.nii.gz")
        if [ ! -f "${tedana_dir}/desc-optcomDenoised_bold.nii.gz" ]; then Error "tedana failed: $(which tedana)"; cleanup "$tmp" "$nocleanup" "$here"; exit 1; fi
        Do_cmd cp -f "${tedana_dir}/desc-optcomDenoised_bold.nii.gz" "$mainScan"
    fi

    # FSL MC outliers
    fmri_mc="${func_volum}/${idBIDS}${func_lab}.1D"
    func_MCoutliers "${fmri_mc}"

    # Distortion correction with TOPUP
    func_topup

else
    Info "Subject ${id} has a functional MRI processed (reoriented/distorion and motion corrected)"; Nsteps=$((Nsteps + 2)); N=$((N+2))
fi

#------------------------------------------------------------------------------#
Info "!!!!!  goin str8 to ICA-FIX yo  !!!!!"

fmri_mean="${func_volum}/${idBIDS}${func_lab}_mean.nii.gz"
fmri_HP="${func_volum}/${idBIDS}${func_lab}_HP.nii.gz"
fmri_brain="${func_volum}/${idBIDS}${func_lab}_brain.nii.gz"
fmri_mask="${func_volum}/${idBIDS}${func_lab}_brain_mask.nii.gz"

if [[ ! -f "$fmri_brain" ]]; then ((N++))
    Info "Generating a func binary mask"
    # Calculates the mean func volume
    Do_cmd fslmaths "$func_nii" -Tmean "$fmri_mean"
    Do_cmd N4BiasFieldCorrection  -d 3 -i "$fmri_mean" -r -o "$fmri_mean"

    # Creates a mask from the motion corrected time series
    Do_cmd mri_synthstrip -i "$fmri_mean" -o "$fmri_brain" -m "$fmri_mask"
    if [[ -f "${fmri_mask}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a brain masked from func"; ((Nsteps++)); ((N++))
fi

# High-pass filter - Remove all frequencies EXCEPT those in the range
if [[ ! -f "${func_volum}/${idBIDS}${func_lab}_tSNR.nii.gz" ]]; then ((N++))
    Info "High pass filter"
    Do_cmd 3dTproject -input "${func_nii}" -prefix "$fmri_HP" -passband 0.01 666
        if [[ -f "${fmri_HP}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has High-pass filter"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# run MELODIC for ICA-FIX
melodic_IC="${func_ICA}/filtered_func_data.ica/melodic_IC.nii.gz"
fmri_filtered="${func_ICA}/filtered_func_data.nii.gz"

# melodic will run ONLY no FIX option is selected
if [[ ! -f "${func_volum}/${idBIDS}${func_lab}_preproc.nii.gz" ]]; then
if [[ "$noFIX" -eq 0 ]] && [[ ! -f "${melodic_IC}" ]]; then
    [[ ! -d "${func_ICA}" ]] && Do_cmd mkdir -p "${func_ICA}"
    Info "Running melodic"
    Do_cmd cp "$fmri_HP" "$fmri_filtered"
    Do_cmd melodic --in="${fmri_filtered}" \
          --tr="${RepetitionTime[0]}" \
          --nobet \
          --mask="${fmri_mask}" \
          --bgthreshold=3 \
          --mmthresh=0.5 \
          --report \
          --Oall \
          --outdir="${func_ICA}/filtered_func_data.ica" \
          --Omean="${func_ICA}/mean_func.nii.gz"
    if [[ -f "${melodic_IC}" ]]; then export statusMel="YES"; else export statusMel="FAILED"; fi
else
    Info "Subject ${id} has MELODIC outputs"; export statusMel="YES"
fi
fi
if [[ "$noFIX" -eq 1 ]]; then export statusMel="NO"; fi
#------------------------------------------------------------------------------#
fmri_in_T1nativepro="${proc_struct}/${idBIDS}_space-nativepro_desc-${tagMRI}_mean.nii.gz"
T1nativepro_in_func="${func_volum}/${idBIDS}_space-func_desc-T1w.nii.gz"
str_func_affine="${dir_warp}/${idBIDS}_from-${tagMRI}_to-nativepro_mode-image_desc-affine_"
mat_func_affine="${str_func_affine}0GenericAffine.mat"
t1bold="${tmp}/${idBIDS}_space-nativepro_desc-T1wbold.nii.gz"

str_func_SyN="${dir_warp}/${idBIDS}_from-nativepro_func_to-${tagMRI}_mode-image_desc-SyN_"
SyN_func_affine="${str_func_SyN}0GenericAffine.mat"
SyN_func_warp="${str_func_SyN}1Warp.nii.gz"
SyN_func_Invwarp="${str_func_SyN}1InverseWarp.nii.gz"
regAffine="FALSE"
if [[ ${regAffine}  == "FALSE" ]]; then
    # SyN from T1_nativepro to t1-nativepro
    export reg="Affine+SyN"
    transformsInv="-t ${SyN_func_warp} -t ${SyN_func_affine} -t [${mat_func_affine},1]" # T1nativepro to func
    transform="-t ${mat_func_affine} -t [${SyN_func_affine},1] -t ${SyN_func_Invwarp}"  # func to T1nativepro
    xfmat="-t ${SyN_func_affine} -t [${mat_func_affine},1]" # T1nativepro to func only lineal for FIX
elif [[ ${regAffine}  == "TRUE" ]]; then
    export reg="Affine"
    transformsInv="-t [${mat_func_affine},1]"  # T1nativepro to func
    transform="-t ${mat_func_affine}"   # func to T1nativepro
    xfmat="-t [${mat_func_affine},1]" # T1nativepro to func only lineal for FIX
fi

# Registration to native pro
Nreg=$(ls "$mat_func_affine" "$fmri_in_T1nativepro" "$T1nativepro_in_func" 2>/dev/null | wc -l )
if [[ "$Nreg" -lt 3 ]]; then ((N++))
    Info "Creating a synthetic T1nativepro image for registration"
    voxels=$(mrinfo "${fmri_mean}" -spacing); voxels="${voxels// /,}"
    Do_cmd flirt -applyisoxfm "${voxels}" -in "${T1nativepro_brain}" -ref "${T1nativepro_brain}" -out "${t1bold}"

    Info "Registering func MRI to nativepro"
    bold_synth="${tmp}/func_brain_synthsegGM.nii.gz"
    t1_synth="${tmp}/T1bold_synthsegGM.nii.gz"
    Do_cmd mri_synthseg --i "${t1bold}" --o "${tmp}/T1bold_synthseg.nii.gz" --robust --threads "$threads" --cpu
    Do_cmd fslmaths "${tmp}/T1bold_synthseg.nii.gz" -uthr 42 -thr 42 -bin -mul -39 -add "${tmp}/T1bold_synthseg.nii.gz" "${t1_synth}"

    Do_cmd mri_synthseg --i "$fmri_brain" --o "${tmp}/func_brain_synthseg.nii.gz" --robust --threads "$threads" --cpu
    Do_cmd fslmaths "${tmp}/func_brain_synthseg.nii.gz" -uthr 42 -thr 42 -bin -mul -39 -add "${tmp}/func_brain_synthseg.nii.gz" "${bold_synth}"

    # Affine from func to t1-nativepro
    centeralign="[${t1_synth},${bold_synth},0]"
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$t1_synth" -m "$bold_synth" -o "$str_func_affine" -t a -n "$threads" -p d -i ${centeralign}
    Do_cmd antsApplyTransforms -d 3 -i "$t1_synth" -r "$bold_synth" -t ["$mat_func_affine",1] -o "${tmp}/T1bold_in_func.nii.gz" -v -u int

    if [[ ${regAffine}  == "FALSE" ]]; then
        # SyN from T1_nativepro to t1-nativepro
        Do_cmd antsRegistrationSyN.sh -d 3 -m "${tmp}/T1bold_in_func.nii.gz" -f "${bold_synth}" -o "$str_func_SyN" -t s -n "$threads" -p d #-i "$mat_func_affine"
    fi
    Do_cmd rm -rf "${dir_warp}"/*Warped.nii.gz 2>/dev/null
    # func to t1-nativepro
    Do_cmd antsApplyTransforms -d 3 -i "$fmri_brain" -r "$t1bold" "${transform}" -o "$fmri_in_T1nativepro" -v -u int
    # t1-nativepro to func
    Do_cmd antsApplyTransforms -d 3 -i "$T1nativepro_brain" -r "$fmri_brain" "${transformsInv}" -o "${T1nativepro_in_func}" -v -u int

    if [[ -d "${func_ICA}/filtered_func_data.ica" ]]; then Do_cmd cp "${T1nativepro_in_func}" "${func_ICA}/filtered_func_data.ica/T1w2fmri_brain.nii.gz"; fi
    if [[ -f "${T1nativepro_in_func}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a func volume and transformation matrix in T1nativepro space"; ((Nsteps++)); ((N++))
fi
proc_func_transformations "${dir_warp}/${idBIDS}_transformations-proc_func-${tagMRI}.json" ${transformsInv// /:} ${transform// /:}

#------------------------------------------------------------------------------#
# run ICA-FIX IF melodic ran succesfully
fix_output="${func_ICA}/filtered_func_data_clean.nii.gz"
func_processed="${func_volum}/${idBIDS}${func_lab}_preproc.nii.gz"

# Run if fmri_clean does not exist
if [[ "$noFIX" -eq 0 ]]; then
    if [[ ! -f "${func_processed}" ]] ; then
          if  [[ -f "${melodic_IC}" ]] && [[ -f $(which fix) ]]; then
              if [[ ! -f "${fix_output}" ]] ; then
                    Info "Getting ICA-FIX requirements"
                    Do_cmd mkdir -p "${func_ICA}"/{reg,mc}
                    # FIX requirements - https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide
                    # $fmri_filtered                                                                                 preprocessed 4D data
                    # $melodic_IC                                                                                    melodic (command-line program) full output directory
                    Do_cmd cp "${fmri_mc}" "${func_ICA}/mc/prefiltered_func_data_mcf.par"   # motion parameters created by mcflirt
                    Do_cmd cp "$fmri_mask" "${func_ICA}/mask.nii.gz"                                                 # valid mask relating to the 4D data
                    Do_cmd cp "${func_ICA}/filtered_func_data.ica/mean.nii.gz" "${func_ICA}/mean_func.nii.gz"      # temporal mean of 4D data
                    middleSlice=$(mrinfo "$fmri_filtered" -size | awk -F ' ' '{printf "%.0f\n", $4/2}')
                    Do_cmd fslroi "$fmri_filtered" "${func_ICA}/reg/example_func.nii.gz" "$middleSlice" 1          # example middle image from 4D data
                    Do_cmd cp "$T1nativepro_brain" "${func_ICA}/reg/highres.nii.gz"                                  # brain-extracted structural

                    # REQUIRED by FIX - reg/highres2example_func.mat                                               # FLIRT transform from structural to functional space
                    if [[ ! -f "${func_ICA}/reg/highres2example_func.mat" ]]; then
                        # Get transformation matrix T1native to func space (ICA-FIX requirement)
                        Do_cmd antsApplyTransforms -v 1 -o Linear["$tmp/highres2example_func.mat",0] "${xfmat}"
                        # Transform matrix: ANTs (itk binary) to text
                        Do_cmd ConvertTransformFile 3 "$tmp/highres2example_func.mat" "$tmp/highres2example_func.txt"

                        # Fixing the transformations incompatibility between ANTS and FSL
                        tmp_ants2fsl_mat="$tmp/itk2fsl_highres2example_func.mat"
                        # Transform matrix: ITK text to matrix (FSL format)
                        Do_cmd lta_convert --initk "$tmp/highres2example_func.txt" --outfsl "$tmp_ants2fsl_mat" --src "$T1nativepro" --trg "$fmri_brain"
                        # apply transformation with FSL
                        Do_cmd flirt -in "$T1nativepro" -out "$tmp/T1w2fmri_brain_ants2fsl.nii.gz" -ref "$fmri_brain" -applyxfm -init "$tmp_ants2fsl_mat"
                        # correct transformation matrix
                        Do_cmd flirt -in "$tmp/T1w2fmri_brain_ants2fsl.nii.gz" -ref "$T1nativepro_in_func" -omat "$tmp/ants2fsl_fixed.omat" -cost mutualinfo -searchcost mutualinfo -dof 6
                        # concatenate the matrices to fix the transformation matrix
                        Do_cmd convert_xfm -concat "${tmp}/ants2fsl_fixed.omat" -omat "${func_ICA}/reg/highres2example_func.mat" "$tmp_ants2fsl_mat"
                    else Info "Subject ${id} has reg/highres2example_func.mat for ICA-FIX"; fi

                    Info "Running ICA-FIX"
                    Do_cmd fix "$func_ICA" "$icafixTraining" 20 -m -h 100

                    # Replace file if melodic ran correctly - Change single-echo files for clean ones
                    if [[ -f "$fix_output" ]]; then
                        yes | Do_cmd 3dresample -orient LPI -prefix "$func_processed" -inset "$fix_output"
                        export statusFIX="YES"
                    else
                        Error "FIX failed, but MELODIC ran log file:\n\t $(ls "${dir_logs}"/proc_func_*.txt)"; cleanup "$tmp" "$nocleanup" "$here"; exit
                    fi
              else
                    Info "Subject ${id} has filtered_func_data_clean from ICA-FIX already"
                    cp -rf "$fix_output" "$func_processed"; export statusFIX="YES"
              fi
          else
              Warning "!!!!  Melodic Failed and/or FIX was not found, check the software installation !!!!
                             If you've installed FIX try to install required R packages and re-run:
                             'kernlab','ROCR','class','party','e1071','randomForest'"
              Do_cmd cp -rf "$fmri_HP" "$func_processed"
              export statusFIX="NO"
          fi
    else
        Info "Subject ${id} has a clean fMRI processed with FIX"; export statusFIX="YES"
    fi
    json_func "${func_volum}/${idBIDS}${func_lab}_preproc.json"
else
    # Skip FIX processing but rename variables anyways for simplicity
    Info "Clean fMRI image has been processed (no FIX)."
    cp -rf "${fmri_HP}" "$func_processed"
    if [[ "$noFIX" -eq 1 ]]; then export statusFIX="NO"; fi
    json_func "${func_volum}/${idBIDS}${func_lab}_preproc.json"
fi

#------------------------------------------------------------------------------#
# Save transformation of the timeseries from T1nativepro_2mm to func space
# func_processed_in_T1nativepro="${func_volum}/${idBIDS}_space-nativepro_2mm_desc-${acq}_clean.nii.gz"
# Do_cmd antsApplyTransforms -d 3 -e 3 -i "$func_processed" -r "$t1bold" "${transform}" -o "$func_processed_in_T1nativepro" -v -u int

#------------------------------------------------------------------------------#
global_signal="${func_volum}/${idBIDS}${func_lab}_global.txt"
if [[ ! -f "${global_signal}" ]]; then ((N++))
    Info "Calculating tissue-specific and global signals changes"
    tissues=(CSF GM WM)
    for idx in "${!tissues[@]}"; do
        tissue=${tissues[$idx]}
        tissuemap="${dir_anat}/${BIDSanat}_space-nativepro_T1w_brain_pve_${idx}.nii.gz"
        tissue_series="${func_volum}/${idBIDS}${func_lab}_pve_${tissue}.txt"
        if [[ ! -f "${tissue_series}" ]] ; then
            Do_cmd antsApplyTransforms -d 3 -i "$tissuemap" -r "$fmri_mean" "${transformsInv}" -o "${tmp}/${idBIDS}${func_lab}_${tissue}.nii.gz" -v -u int
            Do_cmd fslmaths "${tmp}/${idBIDS}${func_lab}_${tissue}.nii.gz" -thr 0.9 "${tmp}/${idBIDS}${func_lab}_${tissue}.nii.gz"
            Do_cmd fslmeants -i "$func_processed" -o "$tissue_series" -m "${tmp}/${idBIDS}${func_lab}_${tissue}.nii.gz" -w
        else
             Info "Subject ${idBIDS} has $tissue time-series"
        fi
    done
    # Global signal from brain mask
    Do_cmd fslmeants -i "$func_processed" -o "$global_signal" -m "${fmri_mask}" -w
    if [[ -f "${global_signal}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has Global time-series"; ((Nsteps++)); ((N++))
fi

# Motion confound (If no output no outliers??)
spikeRegressors="${func_volum}/${idBIDS}${func_lab}_spikeRegressors_REFRMS.1D"
if [[ ! -f "$spikeRegressors" ]] ; then
    Do_cmd fsl_motion_outliers -i "$func_processed" -o "$spikeRegressors" -s "${func_volum}/${idBIDS}${func_lab}_metric_REFRMS.1D" --refmse --nomoco
else
    Info "Subject ${id} has a spike Regressors from fsl_motion_outliers"
fi

#------------------------------------------------------------------------------#
#                                 C O R T E X
surf_dir="${subject_dir}/surf"

# Transform surface to func space
if [[ ! -f "${func_surf}/${idBIDS}_hemi-R_surf-fsnative.func.gii" ]]; then
    # convert affines
    Do_cmd c3d_affine_tool -itk "${mat_func_affine}" -o "$tmp/affine1.mat"
    mat_func_affine=$tmp/affine1.mat
    if [[ ${regAffine}  == "FALSE" ]]; then
        Do_cmd c3d_affine_tool -itk "${SyN_func_affine}" -o "$tmp/affine2.mat" -inv
        SyN_func_affine="$tmp/affine2.mat"
    fi
    # apply registrations to surface
    for HEMICAP in L R; do
        Do_cmd wb_command -surface-apply-affine \
            "${surf_dir}/${idBIDS}_hemi-${HEMICAP}"_space-nativepro_surf-fsnative_label-midthickness.surf.gii \
            "${mat_func_affine}" \
            "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii
        if [[ ${regAffine}  == "FALSE" ]]; then
            Do_cmd wb_command -surface-apply-affine \
                "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii \
                "${SyN_func_affine}" \
                "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii
            Do_cmd wb_command -surface-apply-warpfield \
                "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii \
                "${SyN_func_warp}" \
                "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii
        fi
        Do_cmd wb_command -volume-to-surface-mapping \
            "${func_processed}" \
            "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii \
            "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_surf-fsnative.func.gii \
            -trilinear
        Info "Subject ${id} hemi-${HEMICAP} mapped to fsnative"; ((Nsteps++)); ((N++))
    done
else
    Info "Subject ${id} was already mapped to fsnative"; Nsteps=$((Nsteps+2)); N=$((N+2))
fi

# Propagate to other surfaces
Nsurf=$(ls "${func_surf}/${idBIDS}"_hemi-*_surf-*.func.gii | wc -l)
SURFLIST='fsLR-5k fsLR-32k fsaverage5'
if [ $Nsurf -lt 8 ]; then
    for HEMICAP in L R; do
        for SURF in $SURFLIST; do
            Do_cmd wb_command -metric-resample \
                "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_surf-fsnative.func.gii \
                "${surf_dir}/${idBIDS}_hemi-${HEMICAP}"_surf-fsnative_label-sphere.surf.gii \
                "${MICAPIPE}/surfaces/${SURF}.${HEMICAP}".sphere.reg.surf.gii \
                ADAP_BARY_AREA \
                "${func_surf}/${idBIDS}_hemi-${HEMICAP}_surf-${SURF}".func.gii \
                -area-surfs \
                "${surf_dir}/${idBIDS}_hemi-${HEMICAP}"_space-nativepro_surf-fsnative_label-midthickness.surf.gii \
                "${surf_dir}/${idBIDS}_hemi-${HEMICAP}"_space-nativepro_surf-"${SURF}"_label-midthickness.surf.gii
            Info "Subject ${id} hemi-${HEMICAP} mapped to $SURF"; ((Nsteps++)); ((N++))
        done
    done
else
    Info "Subject ${id} was already mapped to all $SURFLIST"; Nsteps=$((Nsteps+6)); N=$((N+6))
fi


# Caculate tSNR
if [[ ! -f "${func_volum}/${idBIDS}${func_lab}_tSNR.nii.gz" ]]; then
    Do_cmd fslmaths "$func_nii" -Tmean "${tmp}/${idBIDS}${func_lab}_mean"
    Do_cmd fslmaths "$func_nii" -Tstd "${tmp}/${idBIDS}${func_lab}_std"
    Do_cmd fslmaths  "${tmp}/${idBIDS}${func_lab}_mean.nii.gz" \
        -div "${tmp}/${idBIDS}${func_lab}_std.nii.gz" \
        "${func_volum}/${idBIDS}${func_lab}_tSNR.nii.gz"
    Info "Subject ${id} volumetric tSNR calculated"; ((Nsteps++)); ((N++))
else
    Info "Subject ${id} was already tSNR calculated (volumetric)";  ((Nsteps++)); ((N++))
fi
if [[ ! -f "${func_surf}/${idBIDS}_surf-fsnative_hemi-R_tSNR.shape.gii" ]]; then
    for HEMICAP in L R; do
        Do_cmd wb_command -volume-to-surface-mapping \
            "${func_volum}/${idBIDS}${func_lab}_tSNR.nii.gz" \
            "${func_surf}/${idBIDS}_hemi-${HEMICAP}"_space-func_surf-fsnative_label-midthickness.surf.gii \
            "${func_surf}/${idBIDS}_surf-fsnative_hemi-${HEMICAP}_tSNR.shape.gii" \
            -trilinear
        Info "Subject ${id} hemi-${HEMICAP} tSNR mapped to fsnative"; ((Nsteps++)); ((N++))
    done
else
    Info "Subject ${id} tSNR was already mapped to fsnative"; Nsteps=$((Nsteps+2)); N=$((N+2))
fi

#------------------------------------------------------------------------------#
#                           S U B C O R T E X
# Subcortical segmentation (nativepro) to func space
func_subcortex="${func_volum}/${idBIDS}${func_lab}_subcortical.nii.gz"
timese_subcortex="${func_volum}/${idBIDS}${func_lab}_timeseries_subcortical.txt"

if [[ ! -f "$timese_subcortex" ]]; then ((N++))
      Info "Getting subcortical timeseries"
      Do_cmd antsApplyTransforms -d 3 -i "$T1_seg_subcortex" -r "$fmri_mean" -n MultiLabel "${transformsInv}" -o "$func_subcortex" -v -u int
      # Extract subcortical timeseries
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i "$func_processed" --seg "$func_subcortex" --exclude 0 --exclude 16 --avgwf "$timese_subcortex"
      if [[ -f "$timese_subcortex" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has func subcortical time-series"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
#                           C E R E B E L L U M
func_cerebellum="${func_volum}/${idBIDS}${func_lab}_cerebellum.nii.gz"
timese_cerebellum="${func_volum}/${idBIDS}${func_lab}_timeseries_cerebellum.txt"
stats_cerebellum="${func_volum}/${idBIDS}${func_lab}_cerebellum_roi_stats.txt"

if [[ ! -f "$timese_cerebellum" ]]; then ((N++))
      Info "Getting cerebellar timeseries"
      Do_cmd antsApplyTransforms -d 3 -i "$T1_seg_cerebellum" -r "$fmri_mean" -n MultiLabel "${transformsInv}" -o "$func_cerebellum" -v -u int
      # Extract cerebellar timeseries (mean, one ts per segemented structure, exluding nuclei because many are too small for our resolution)
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i "$func_processed" --seg "$func_cerebellum" --exclude 0 --avgwf "$timese_cerebellum"
      3dROIstats -mask "$func_cerebellum" -nzsum "$func_cerebellum" > "$stats_cerebellum"
      if [[ -f "$timese_cerebellum" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has func cerebellar time-series"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# run post-func
cleanTS="${func_surf}/${idBIDS}_surf-fsLR-32k_desc-timeseries_clean.shape.gii"
if [[ ! -f "$cleanTS" ]]; then ((N++))
    Info "Running func post processing"
    labelDirectory="${MICAPIPE}/parcellations/"
    Do_cmd python "$MICAPIPE"/functions/03_FC.py "$idBIDS" "$proc_func" "$labelDirectory" "$util_parcelations" "$dir_volum" "$performNSR" "$performGSR" "$func_lab" "$noFC" "$GSR"
    if [[ -f "$cleanTS" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has post-processed fsLR time-series"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# a bit of extra cleanup
rm "${func_volum}/${idBIDS}"_*brain_mask.nii.gz "${func_volum}/${idBIDS}"_*HP.nii.gz "${func_volum}/${idBIDS}"_*mean.nii.gz "${func_surf}/${idBIDS}"_*surf-fsnative_NoHP.func.gii "${func_volum}/${idBIDS}${func_lab}".nii.gz 2>/dev/null

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status "proc_func"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_func_${tagMRI}" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "proc_func_${tagMRI}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
