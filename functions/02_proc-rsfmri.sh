#!/bin/bash
# Resting state preprocessing
# Written by Casey Paquola and Reinder Vos De Wael (Oct 2018).
# and a tiny bit from Sara (Feb 2019)...
# and a whole lot from Sara (August 2019)
# and incorporation to mica-pipe by Raul (August-September 2020)
# and addition of a bunch of fancy flags by Jessica (October-November 2020)
#
# Resting state fMRI processing with bash:
#
# Preprocessing workflow for rsfmri.
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
umask 003
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
smooth=${12}
mainScanStr=${13}
fmri_pe=${14}
fmri_rpe=${15}
performNSR=${16}
noFIX=${17}
PROC=${18}
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source ${MICAPIPE}/functions/init.sh;
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"
T1_seg_subcortex="${dir_volum}/${idBIDS}_space-nativepro_t1w_atlas-subcortical.nii.gz"
T1_seg_cerebellum="${dir_volum}/${idBIDS}_space-nativepro_t1w_atlas-cerebellum.nii.gz"

### CHECK INPUTS: rsfMRI, phase encoding, structural proc, topup and ICA-FIX files
Info "Inputs:"
Note "Topup Config     :" "$changeTopupConfig"
Note "ICA fix training :" "$changeIcaFixTraining"
if [[ "$mainScanStr" == DEFAULT ]]; then Note "Main scan        :" "$thisMainScan"; else
Note "Main scan        :" $(ls "${subject_bids}/func/${idBIDS}"_"${mainScanStr}".nii* 2>/dev/null); fi
Note "Phase scan       :" "$fmri_pe"
Note "Reverse Phase    :" "$fmri_rpe"

#------------------------------------------------------------------------------#
if [[ "$mainScanStr" == DEFAULT ]]; then
    # Main scan
    N_mainScan=${#bids_mainScan[@]}
    if [ "$N_mainScan" -gt 1 ]; then
        if [[ "${thisMainScan}" == "DEFAULT" ]]; then
            Error "Multiple rsfMRI runs found in BIDS rawdata directory! Please specify which run should be processed using flag -mainScanRun"; exit;
        elif [ "$thisMainScan" -gt "$N_mainScan" ]; then
            Warning "Specified run number (${thisMainScan}) is greater than number of rsfMRI scans scans found ($N_mainScan). Using first filename in list as default";
            mainScan=${bids_mainScan[0]}
        else
            Info "Found $N_mainScan rsfMRI scans, processing specified scan # $thisMainScan"
            mainScan=${bids_mainScan[$thisMainScan-1]}
        fi
    else
        mainScan=${bids_mainScan[0]}
        if [[ "$thisMainScan" == "DEFAULT" ]]; then
            Info "No run number specified for rsfMRI scan and did not find more than one run for main scan - all good!"
        else
            if [ "$thisMainScan" -gt "$N_mainScan" ]; then
                Warning "Found one or less rsfMRI scan, but specified run number = $thisMainScan). Using first filename in list as default";
            fi
        fi
    fi

    # Main scan json
    N_mainScanJson=${#bids_mainScanJson[@]}
    if [ "$N_mainScanJson" -gt 1 ]; then
        if [[ "${thisMainScan}" == "DEFAULT" ]]; then
            Error "Found multiple .json files for main rsfMRI scan in BIDS rawdata directory! Please specify which run should be processed using flag -mainScanRun"; exit;
        elif [ "$thisMainScan" -gt "$N_mainScanJson" ]; then
            Warning "Specified run number (${thisMainScan}) is greater than number of rsfMRI json files found for main scan ($N_mainScan). Using first filename in list as default";
            mainScanJson=${bids_mainScan[0]}
        else
            Info "Found ${N_mainScanJson} rsfMRI scan json files, using specified run # ${thisMainScan}"
            mainScanJson=${bids_mainScanJson[$thisMainScan-1]}
        fi
    else
        Info "Using default json scan: ${bids_mainScanJson[0]}"
        mainScanJson=${bids_mainScanJson[0]}
    fi
else
    mainScan=$(ls "${subject_bids}/func/${idBIDS}_${mainScanStr}".nii* 2>/dev/null)
    mainScanJson=$(ls "${subject_bids}/func/${idBIDS}_${mainScanStr}".json 2>/dev/null)
fi
#------------------------------------------------------------------------------#
# Phase encoding
N_mainPhase=${#bids_mainPhase[@]}
N_revPhase=${#bids_reversePhase[@]}
if [ "$N_mainPhase" -gt 1 ] || [ "$N_revPhase" -gt 1 ]; then
    if [[ "$thisPhase" == "DEFAULT" ]]; then
        Error "Found multiple phase reversal runs in BIDS rawdata directory! Please specify which run should be processed using flag -phaseReversalRun"; exit;
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
if [[ "$fmri_pe" != DEFAULT ]] && [[ -f "$fmri_rpe" ]]; then mainPhaseScan="$fmri_rpe"; fi
if [[ "$fmri_rpe" != DEFAULT ]] && [[ -f "$fmri_rpe" ]]; then reversePhaseScan="$fmri_rpe"; fi

# Check inputs
if [ ! -f "$mainScan" ]; then Error "Couldn't find $id main rsfMRI scan : \n\t ls ${mainScan}"; exit; fi #Last check to make sure file exists
if [ ! -f "$mainScanJson" ]; then Error "Couldn't find $id main rsfMRI scan json file: \n\t ls ${mainScanJson}"; exit; fi #Last check to make sure file exists
if [ -z "$mainPhaseScan" ]; then  Warning "Subject $id doesn't have acq-APse_bold: TOPUP will be skipped"; fi
if [ -z "$reversePhaseScan" ]; then Warning "Subject $id doesn't have acq-PAse_bold: TOPUP will be skipped"; fi

# Check requirements: Structural nativepro scan and freesurfer, and post_structural
if [ ! -f "$T1nativepro" ]; then Error "Subject $id doesn't have T1_nativepro: run -proc_structural"; exit; fi
if [ ! -f "$T1freesurfr" ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${idBIDS}/mri/T1.mgz"; exit; fi
if [ ! -f "$T1_seg_cerebellum" ]; then Error "Subject $id doesn't have cerebellar segmentation:\n\t\t ls ${T1_seg_cerebellum} \n\t\tRUN -post_structural"; exit; fi
if [ ! -f "$T1_seg_subcortex" ]; then Error "Subject $id doesn't have subcortical segmentation:\n\t\t ls ${T1_seg_subcortex} \n\t\t -post_structural"; exit; fi

# Check topup input
if [[ ${changeTopupConfig} == "DEFAULT" ]]; then
    Info "Will use default config file for TOPUP: ${topupConfigFile}"
else
    topupConfigFile=${changeTopupConfig}
    Info "Will use specified config file for TOPUP: ${topupConfigFile}"
fi

# Check FIX: run or no?
if [[ $noFIX == 1 ]]; then
    Info "ICA-FIX will be skipped! Consider performing white matter and CSF signal regression with <-nuisanceRegression>"

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

# Check smoothing
if [[ $smooth == 1 ]]; then
    Info "Smoothing of native surface timeseries will be performed using workbench command"
else
    Info "Smoothing of native surface timeseries will be performed using FreeSurfer tools (default)"
fi

# Check nuisance signal regression
if [[ $performNSR == 1 ]]; then
    Info "White matter and CSF signal will be regressed from processed timeseries"
else
    Info "White matter and CSF signal regression will not be performed (default)"
fi


#------------------------------------------------------------------------------#
Title "Resting state fMRI processing\n\t\tmicapipe $Version, $PROC "
micapipe_software
bids_print.variables-rsfmri
Info "Saving temporal dir: $nocleanup"
Info "ANTs will use $threads threads"
Info "wb_command will use $OMP_NUM_THREADS threads"

#	Timer
aloita=$(date +%s)
Nsteps=0
# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_proc-rsfmri_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Define directories
export SUBJECTS_DIR="$dir_surf"

# rsfMRI directories
rsfmri_volum="${proc_rsfmri}/volumetric"   # volumetricOutputDirectory
rsfmri_surf="${proc_rsfmri}/surfaces"      # surfaceOutputDirectory
rsfmri_ICA="$proc_rsfmri/ICA_MELODIC"      # ICAOutputDirectory

# Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script).
for x in "$rsfmri_surf" "$rsfmri_volum" "$rsfmri_ICA"; do
    [[ ! -d "${x}" ]] && mkdir -p "${x}"
done

#------------------------------------------------------------------------------#
# Begining of the REAL processing
# gettin dat readout time mainScanJson
readoutTime=$(grep TotalReadoutTime "${mainScanJson}" | grep -Eo [0-9].[0-9]+)
RepetitionTime=$(grep RepetitionTime "${mainScanJson}" | grep -Eo [0-9].[0-9]+)

# Scans to process
toProcess=($mainScan $mainPhaseScan $reversePhaseScan)
tags=(mainScan mainPhaseScan reversePhaseScan)
singleecho="${rsfmri_volum}/${idBIDS}"_space-rsfmri_desc-singleecho.nii.gz

# scan for registration to T1
# rsfmri4reg="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_mainPhaseAlignedTopup_mean.nii.gz"

# Processing single.
if [[ ! -f "${singleecho}" ]]; then
    # Loop over all scans for everything before motion correction across scans.
    for i in {0..2}; do
        # Get basic parameters
        rawNifti=${toProcess[$i]}
        tag=${tags[$i]}
        Info "Processing TAG-$tag scan, readout time: $readoutTime ms"
        # IF FILE NOT FOUND DON'T RUN
        if [[ ! -z "${rawNifti}" ]] && [[ -f "${rawNifti}" ]]; then
              Note "RAWNIFTI:" "$rawNifti"

              # Drop first five TRs and reorient to standard
              if [ "$tag" == "mainScan" ]; then
                  Do_cmd nifti_tool -cbl -prefix "${tmp}/${tag}_trDrop.nii.gz" -infiles "$rawNifti"'[5..$]'
                  Do_cmd 3dresample -orient LPI -prefix "${tmp}/${tag}_reorient.nii.gz" -inset "${tmp}/${tag}_trDrop.nii.gz"
                  Do_cmd fslreorient2std "${tmp}/${tag}_reorient.nii.gz" "${tmp}/${tag}_reorient.nii.gz"
              else
                  Do_cmd 3dresample -orient LPI -prefix "${tmp}/${tag}_reorient.nii.gz" -inset "$rawNifti"
                  Do_cmd fslreorient2std "${tmp}/${tag}_reorient.nii.gz" "${tmp}/${tag}_reorient.nii.gz"
              fi

              # Remove slices to make an even number of slices in all directions (requisite for topup).
              #dimensions=`fslhd ${tmp}/${tag}_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
              #newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
              #Do_cmd fslroi ${tmp}/${tag}_reorient.nii.gz ${tmp}/${tag}_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` # I removed -1 -1

              # Skipping fslroi step. Rename files for simplicity
              mv "${tmp}/${tag}_reorient.nii.gz" "${tmp}/${tag}_sliceCut.nii.gz"

              # Motion correction within scans
              Do_cmd fslmaths "${tmp}/${tag}_sliceCut.nii.gz" -Tmean "${tmp}/${tag}_sliceCutMean.nii.gz"
              Do_cmd 3dvolreg -Fourier -twopass -base "${tmp}/${tag}_sliceCutMean.nii.gz" \
                              -zpad 4 -prefix "${tmp}/${tag}_mc.nii.gz" \
                              -1Dfile "${rsfmri_volum}/${idBIDS}_space-rsfmri_${tag}.1D" \
                              "${tmp}/${tag}_sliceCut.nii.gz"
              Do_cmd fslmaths "${tmp}/${tag}_mc.nii.gz" -Tmean "${tmp}/${tag}_mcMean.nii.gz"
        fi
    done

    # Calculate motion outliers with FSL
    if [[ ! -f "${rsfmri_volum}/${idBIDS}_space-rsfmri_singleecho.1D" ]]; then
        Do_cmd fsl_motion_outliers -i "${tmp}/mainScan_sliceCut.nii.gz" \
                                   -o "${rsfmri_volum}/${idBIDS}_space-rsfmri_spikeRegressors_FD.1D" \
                                   -s "${rsfmri_volum}/${idBIDS}_space-rsfmri_metric_FD.1D" --fd
        Do_cmd mv "${rsfmri_volum}/${idBIDS}_space-rsfmri_mainScan.1D ${rsfmri_volum}/${idBIDS}_space-rsfmri_singleecho.1D"
    else
        Info "Subject ${id} has a singleecho.1D with motion outliers"; ((Nsteps++))
    fi

    # If ONLY Reverse phase scan is provided mainScan will be the mainPhaseScan
    if [[ -f "$reversePhaseScan" ]] && [[ ! -f "$mainPhaseScan" ]]; then;
          main_pe="NO-main_pe"
          Warning "reversePhaseScan was found but NO mainPhaseScan, using mainScan as mainPhaseScan"
          mainPhaseScan="${tmp}/mainPhaseScan_mc.nii.gz"
          Do_cmd cp "${tmp}/mainScan_mc.nii.gz" "$mainPhaseScan"
          Do_cmd cp "${tmp}/mainScan_mcMean.nii.gz" "${tmp}/singleecho_mainPhaseAlignedMean.nii.gz"
    fi

    # Only do distortion correction if field maps were provided, if not then rename the scan to distortionCorrected (just to make the next lines of code easy).
    if [ -z "${mainPhaseScan}" ] || [ -z "${reversePhaseScan}" ]; then
        Warning "No AP or PA acquisition was found, TOPUP will be skip!!!!!!!"
        export statusTopUp="NO"
        Do_cmd mv -v "${tmp}/mainScan_mc.nii.gz" "${singleecho}"
    else
        if [[ ! -f "${rsfmri_volum}/TOPUP.txt" ]] && [[ ! -f "${singleecho}" ]]; then
            mainPhaseScanMean=$(find "$tmp"    -maxdepth 1 -name "*mainPhaseScan*_mcMean.nii.gz")
            mainPhaseScan=$(find "$tmp"        -maxdepth 1 -name "*mainPhaseScan*_mc.nii.gz")
            reversePhaseScanMean=$(find "$tmp" -maxdepth 1 -name "*reversePhaseScan*_mcMean.nii.gz")
            reversePhaseScan=$(find "$tmp"     -maxdepth 1 -name "*reversePhaseScan*_mc.nii.gz")
            mainScan=$(find "$tmp"             -maxdepth 1 -name "*mainScan*_mc.nii.gz")

            Do_cmd flirt -in "$reversePhaseScanMean" -ref "$tmp"/mainScan_mcMean.nii.gz -omat "$tmp"/singleecho_tmpXfmSecondary.omat
            Do_cmd flirt -in "$reversePhaseScan" -ref "$tmp"/mainScan_mcMean.nii.gz -applyxfm -init "$tmp"/singleecho_tmpXfmSecondary.omat -out "$tmp"/singleecho_secondaryPhaseAligned.nii.gz
            Do_cmd fslmaths "$tmp"/singleecho_secondaryPhaseAligned.nii.gz -Tmean "$tmp"/singleecho_secondaryPhaseAlignedMean.nii.gz

            if [[ "$main_pe" != "NO-main_pe" ]]; then
                Do_cmd flirt -in "$mainPhaseScanMean" -ref "$tmp"/mainScan_mcMean.nii.gz -omat "$tmp"/singleecho_tmpXfmMain.omat
                Do_cmd flirt -in "$mainPhaseScan" -ref "$tmp"/mainScan_mcMean.nii.gz -applyxfm -init "$tmp"/singleecho_tmpXfmMain.omat -out "$tmp"/singleecho_mainPhaseAligned.nii.gz
                Do_cmd fslmaths "$tmp"/singleecho_mainPhaseAligned.nii.gz -Tmean "$tmp"/singleecho_mainPhaseAlignedMean.nii.gz
            fi

            # Distortion correction
            echo -e "0 1 0 ${readoutTime} \n0 -1 0 ${readoutTime}" > "$tmp"/singleecho_topupDataIn.txt
            Info "topup datain:\n$(cat "${tmp}"/singleecho_topupDataIn.txt)"
            Do_cmd fslmerge -t "${tmp}/singleecho_mergeForTopUp.nii.gz" "${tmp}/singleecho_mainPhaseAlignedMean.nii.gz" "${tmp}/singleecho_secondaryPhaseAlignedMean.nii.gz"
            Do_cmd topup --imain="${tmp}/singleecho_mergeForTopUp.nii.gz" --datain="${tmp}/singleecho_topupDataIn.txt" --config="${topupConfigFile}" --out="$tmp/singleecho_topup"
            Do_cmd applytopup --imain="${mainScan}" --inindex=1 --datain="${tmp}/singleecho_topupDataIn.txt" --topup="${tmp}/singleecho_topup" --method=jac --out="${singleecho}"

            # # rsfmri for registration
            # Do_cmd applytopup --imain="${tmp}/singleecho_mainPhaseAligned.nii.gz" --inindex=1 --datain="${tmp}/singleecho_topupDataIn.txt" --topup="${tmp}/singleecho_topup" --method=jac --out="${tmp}/singleecho_mainPhaseAlignedTopup.nii.gz"
            # Do_cmd fslmaths "${tmp}/singleecho_mainPhaseAlignedTopup.nii.gz" -Tmean "$rsfmri4reg"

            # Check if it worked
            if [[ ! -f "${singleecho}" ]]; then Error "Something went wrong with TOPUP check ${tmp} and log:\n\t\t${dir_logs}/proc_rsfmri.txt"; exit; fi; ((Nsteps++))
            export statusTopUp="YES"
        else
            Info "Subject ${id} has singleecho in fmrispace with TOPUP"; export statusTopUp="YES"; ((Nsteps++))
        fi
    fi
else
      Info "Subject ${id} has a singleecho_fmrispace processed"; Nsteps=$((Nsteps + 2))
fi
json_rsfmri "${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_clean.json"

#------------------------------------------------------------------------------#
Info "!!!!!  goin str8 to ICA-FIX yo  !!!!!"

fmri_mean="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_mean.nii.gz"
fmri_mask="${rsfmri_ICA}/mask.nii.gz"
fmri_HP="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_HP.nii.gz"
fmri_brain="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_brain.nii.gz"

if [[ ! -f "$fmri_mask" ]] || [[ ! -f "$fmri_brain" ]]; then
    Info "Generating a rsfMRI binary mask"
    # Calculates the mean rsfMRI volume
    Do_cmd fslmaths "$singleecho" -Tmean "$fmri_mean"

    # Creates a mask from the motion corrected time series
    Do_cmd bet "$fmri_mean ${rsfmri_ICA}/func.nii.gz" -m -n
    Do_cmd mv "${rsfmri_ICA}/func_mask.nii.gz" "$fmri_mask"

    # masked mean rsfMRI time series
    Do_cmd fslmaths "$fmri_mean" -mul "$fmri_mask" "$fmri_brain"
    if [[ -f "${fmri_mask}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a binary mask of the rsfMRI"; ((Nsteps++))
fi

# High-pass filter - Remove all frequencies EXCEPT those in the range
if [[ ! -f "$fmri_HP" ]]; then
    Info "High pass filter"
    Do_cmd 3dTproject -input "${singleecho}" -prefix "$fmri_HP" -passband 0.01 666
        if [[ -f "${fmri_HP}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has High-pass filter"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# run MELODIC for ICA-FIX
melodic_IC="${rsfmri_ICA}/filtered_func_data.ica/melodic_IC.nii.gz"
fmri_filtered="${rsfmri_ICA}/filtered_func_data.nii.gz"

# melodic will run ONLY if FIX is avaliable
if  [[ -f $(which fix) ]]; then
      if [[ ! -f "${melodic_IC}" ]]; then
          Info "Running melodic"
          Do_cmd cp "$fmri_HP" "$fmri_filtered"
          Do_cmd melodic --in="${fmri_filtered}" \
                          --tr="$RepetitionTime" \
                          --nobet \
                          --mask="${fmri_mask}" \
                          --bgthreshold=3 \
                          --mmthresh=0.5 \
                          --report \
                          --Oall \
                          --outdir="${rsfmri_ICA}/filtered_func_data.ica" \
                          --Omean="${rsfmri_ICA}/mean_func.nii.gz"
          if [[ -f "${melodic_IC}" ]]; then export statusMel="YES"; else export statusMel="FAILED"; fi
      else
          Info "Subject ${id} has MELODIC outputs"
      fi
fi

#------------------------------------------------------------------------------#
fmri_in_T1nativepro="${proc_struct}/${idBIDS}_space-nativepro_desc-rsfmri_bold.nii.gz"
T1nativepro_in_fmri="${rsfmri_ICA}/filtered_func_data.ica/t1w2fmri_brain.nii.gz"
str_rsfmri_affine="${dir_warp}/${idBIDS}_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_"
mat_rsfmri_affine="${str_rsfmri_affine}0GenericAffine.mat"
t1bold="${proc_struct}/${idBIDS}_space-nativepro_desc-t1wbold.nii.gz"

str_rsfmri_SyN="${dir_warp}/${idBIDS}_rsfmri_from-nativepro_to-nativepro_mode-image_desc-SyN_"
SyN_rsfmri_affine="${str_rsfmri_SyN}0GenericAffine.mat"
SyN_rsfmri_warp="${str_rsfmri_SyN}1Warp.nii.gz"
SyN_rsfmri_Invwarp="${str_rsfmri_SyN}1InverseWarp.nii.gz"

# Registration to native pro
if [[ ! -f "$mat_rsfmri_affine" ]] || [[ ! -f "$fmri_in_T1nativepro" ]]; then
    # if [[ -f "$rsfmri4reg" ]]; then
    #     Do_cmd fslmaths "$rsfmri4reg" -mul "$fmri_mask" "$fmri_brain"
    # fi

    Info "Creating a synthetic BOLD image for registration"
    # Inverse T1w
    Do_cmd ImageMath 3 "${tmp}/${id}_t1w_nativepro_NEG.nii.gz" Neg "$T1nativepro"
    # Dilate the T1-mask
    #Do_cmd ImageMath 3 "${tmp}/${id}_t1w_mask_dil-2.nii.gz" MD "$T1nativepro_mask" 2
    # Masked the inverted T1w
    Do_cmd ImageMath 3 "${tmp}/${id}_t1w_nativepro_NEG_brain.nii.gz" m "${tmp}/${id}_t1w_nativepro_NEG.nii.gz" "$T1nativepro_mask"
    # Match histograms values acording to rsfmri
    Do_cmd ImageMath 3 "${tmp}/${id}_t1w_nativepro_NEG-rescaled.nii.gz" HistogramMatch "${tmp}/${id}_t1w_nativepro_NEG_brain.nii.gz" "$fmri_brain"
    # Smoothing
    Do_cmd ImageMath 3 "$t1bold" G "${tmp}/${id}_t1w_nativepro_NEG-rescaled.nii.gz" 1

    Info "Registering fmri space to nativepro"
    # Affine from rsfMRI to t1-nativepro
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1nativepro_brain" -m "$fmri_brain" -o "$str_rsfmri_affine" -t a -n "$threads" -p d

    # SyN from T1_nativepro to t1-nativepro
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$t1bold" -m "${str_rsfmri_affine}Warped.nii.gz" -o "$str_rsfmri_SyN" -t s -n "$threads" -p d -i "$mat_rsfmri_affine"

    # fmri to t1-nativepro
    Do_cmd antsApplyTransforms -d 3 -i "$fmri_brain" -r "$T1nativepro_brain" -t "$SyN_rsfmri_warp" -t "$SyN_rsfmri_affine" -t "$mat_rsfmri_affine" -o "$fmri_in_T1nativepro" -v -u int

    # t1-nativepro to fmri
    Do_cmd antsApplyTransforms -d 3 -i "$T1nativepro" -r "$fmri_brain" -t ["$mat_rsfmri_affine",1] -t ["$SyN_rsfmri_affine",1] -t "$SyN_rsfmri_Invwarp" -o "${T1nativepro_in_fmri}" -v -u int
    Do_cmd cp "${T1nativepro_in_fmri}" "${rsfmri_volum}/${idBIDS}_space-rsfmri_t1w.nii.gz"
    if [[ -f "${SyN_rsfmri_Invwarp}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a rsfMRI volume and transformation matrix in T1nativepro space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Register rsfMRI to Freesurfer space with Freesurfer
fmri2fs_dat="${dir_warp}/${idBIDS}_from-rsfmri_to-fsnative_bbr.dat"
if [[ ! -f "${fmri2fs_dat}" ]] ; then
  Info "Registering fmri to FreeSurfer space"
    Do_cmd bbregister --s "$idBIDS" --mov "$fmri_mean" --reg "${fmri2fs_dat}" --o "${dir_warp}/${idBIDS}_from-rsfmri_to-fsnative_bbr_outbbreg_FIX.nii.gz" --init-fsl --bold
    if [[ -f "${fmri2fs_dat}" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a dat transformation matrix from fmri to Freesurfer space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# run ICA-FIX IF melodic ran succesfully
fix_output="${rsfmri_ICA}/filtered_func_data_clean.nii.gz"
fmri_processed="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_clean.nii.gz"

# Run if fmri_clean does not exist
if [[ $noFIX == 0 ]]; then
    if [[ ! -f "${fmri_processed}" ]] ; then
          if  [[ -f "${melodic_IC}" ]] && [[ -f $(which fix) ]]; then
              if [[ ! -f "${fix_output}" ]] ; then
                    Info "Getting ICA-FIX requirements"
                    # FIX requirements - https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide
                    if [ ! -d "${rsfmri_ICA}/mc" ]; then mkdir "${rsfmri_ICA}/mc"; fi
                    if [ ! -d "${rsfmri_ICA}/reg" ]; then mkdir "${rsfmri_ICA}/reg"; fi

                    # $fmri_filtered                                                                                 preprocessed 4D data
                    # $melodic_IC                                                                                    melodic (command-line program) full output directory
                    Do_cmd cp "${rsfmri_volum}/${idBIDS}_space-rsfmri_singleecho.1D" "${rsfmri_ICA}/mc/prefiltered_func_data_mcf.par"   # motion parameters created by mcflirt
                    # $fmri_mask                                                                                     valid mask relating to the 4D data
                    Do_cmd cp "${rsfmri_ICA}/filtered_func_data.ica/mean.nii.gz" "${rsfmri_ICA}/mean_func.nii.gz"      # temporal mean of 4D data
                    middleSlice=$(mrinfo "$fmri_filtered" -size | awk -F ' ' '{printf "%.0f\n", $4/2}')
                    Do_cmd fslroi "$fmri_filtered" "${rsfmri_ICA}/reg/example_func.nii.gz" "$middleSlice" 1          # example middle image from 4D data
                    Do_cmd cp "$T1nativepro_brain" "${rsfmri_ICA}/reg/highres.nii.gz"                                  # brain-extracted structural

                    # REQUIRED by FIX - reg/highres2example_func.mat                                               # FLIRT transform from structural to functional space
                    if [[ ! -f "${rsfmri_ICA}/reg/highres2example_func.mat" ]]; then
                        # Get transformation matrix T1native to rsfMRI space (ICA-FIX requirement)
                        Do_cmd antsApplyTransforms -v 1 -o Linear["$tmp/highres2example_func.mat",0] -t ["$mat_rsfmri_affine",1] -t ["$SyN_rsfmri_affine",1]
                        # Transform matrix: ANTs (itk binary) to text
                        Do_cmd ConvertTransformFile 3 "$tmp/highres2example_func.mat" "$tmp/highres2example_func.txt"

                        # Fixing the transformations incompatibility between ANTS and FSL
                        tmp_ants2fsl_mat="$tmp/itk2fsl_highres2example_func.mat"
                        # Transform matrix: ITK text to matrix (FSL format)
                        Do_cmd lta_convert --initk "$tmp/highres2example_func.txt" --outfsl "$tmp_ants2fsl_mat" --src "$T1nativepro" --trg "$fmri_brain"
                        # apply transformation with FSL
                        Do_cmd flirt -in "$T1nativepro" -out "$tmp/t1w2fmri_brain_ants2fsl.nii.gz" -ref "$fmri_brain" -applyxfm -init "$tmp_ants2fsl_mat"
                        # correct transformation matrix
                        Do_cmd flirt -in "$tmp/t1w2fmri_brain_ants2fsl.nii.gz" -ref "$T1nativepro_in_fmri" -omat "$tmp/ants2fsl_fixed.omat" -cost mutualinfo -searchcost mutualinfo -dof 6
                        # concatenate the matrices to fix the transformation matrix
                        Do_cmd convert_xfm -concat "$tmp"/ants2fsl_fixed.omat -omat "${rsfmri_ICA}/reg/highres2example_func.mat" "$tmp_ants2fsl_mat"
                    else Info "Subject ${id} has reg/highres2example_func.mat for ICA-FIX"; fi

                    Info "Running ICA-FIX"
                    Do_cmd fix "$rsfmri_ICA" "$icafixTraining" 20 -m -h 100

                    # Replace file if melodic ran correctly - Change single-echo files for clean ones
                    if [[ -f "$fix_output" ]]; then
                        yes | Do_cmd cp -rf "$fix_output" "$fmri_processed"
                        statusFIX="YES"
                    else
                        Error "FIX failed, but MELODIC ran log file:\n\t${dir_logs}/proc_rsfmri.txt"; exit
                    fi
              else
                    Info "Subject ${id} has filtered_func_data_clean from ICA-FIX already"
                    cp -rf "$fix_output" "$fmri_processed"; statusFIX="YES"
              fi
          else
              Warning "!!!!  Melodic Failed and/or FIX was not found, check the software installation !!!!
                             If you've installed FIX try to install required R packages and re-run:
                             'kernlab','ROCR','class','party','e1071','randomForest'"
              Do_cmd cp -rf "$fmri_HP" "$fmri_processed" # OR cp -rf  $singleecho $fmri_processed <<<<<<<<<<<<<<<<<<<<< NOT SURE YET
              statusFIX="$NO"
          fi
    else
        Info "Subject ${id} has singleecho_fmrispace_clean volume"
    fi
else
    # Skip FIX processing but rename variables anyways for simplicity
    Info "Further processing will be performed on distorsion corrected images."
    cp -rf "${singleecho}" "$fmri_processed"
fi
json_rsfmri "${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_clean.json"
#------------------------------------------------------------------------------#
global_signal="${rsfmri_volum}/${idBIDS}_space-rsfmri_global.txt"
if [[ ! -f "${global_signal}" ]] ; then
      Info "Calculating tissue-specific and global signals changes"
      tissues=(CSF GM WM)
      for idx in {0..2}; do
           tissue=${tissues[$idx]}
           tissuemap="${proc_struct}/${idBIDS}_space-nativepro_t1w_brain_pve_${idx}.nii.gz"
           tissue_series="${rsfmri_volum}/${idBIDS}_space-rsfmri_pve_${tissue}.txt"
           if [[ ! -f "${tissue_series}" ]] ; then
           Do_cmd antsApplyTransforms -d 3 -i "$tissuemap" -r "$fmri_mean" -t ["$mat_rsfmri_affine",1] -t ["$SyN_rsfmri_affine",1] -t "$SyN_rsfmri_Invwarp" -o "${tmp}/${idBIDS}_space-rsfmri_${tissue}.nii.gz" -v -u int
           Do_cmd fslmeants -i "$fmri_processed" -o "$tissue_series" -m "${tmp}/${idBIDS}_space-rsfmri_${tissue}.nii.gz" -w
         else
             Info "Subject ${idBIDS} has $tissue time-series"
         fi
      done
      Do_cmd fslmaths "${tmp}/${idBIDS}_space-rsfmri_WM.nii.gz" -add  "${tmp}/${idBIDS}_space-rsfmri_GM.nii.gz" -add  "${tmp}/${idBIDS}_space-rsfmri_CSF.nii.gz" "${tmp}/${idBIDS}_space-rsfmri_WB.nii.gz"
      Do_cmd fslmeants -i "$fmri_processed" -o "$global_signal" -m "${tmp}/${idBIDS}_space-rsfmri_WB.nii.gz" -w
      if [[ -f "${global_signal}" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has Global time-series"; ((Nsteps++))
fi

# Motion confound
spikeRegressors="${rsfmri_volum}/${idBIDS}_space-rsfmri_spikeRegressors_REFRMS.1D"
if [[ ! -f "$spikeRegressors" ]] ; then
    Do_cmd fsl_motion_outliers -i "$fmri_processed" -o "$spikeRegressors" -s "${rsfmri_volum}/${idBIDS}_space-rsfmri_metric_REFRMS.1D" --refmse --nomoco
    if [[ -f "$spikeRegressors" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a spike Regressors from fsl_motion_outliers"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Register to surface
# If three surfaces are found skipp this step
Nsurf=$(ls "${rsfmri_surf}/${idBIDS}"_rsfmri_space-fsnative_?h.mgh \
            "${rsfmri_surf}/${idBIDS}"_rsfmri_space-fsnative_?h_10mm.mgh \
            "${rsfmri_surf}/${idBIDS}"_rsfmri_space-fsaverage5_?h_10mm.mgh \
            "${rsfmri_surf}/${idBIDS}"_rsfmri_space-conte69-32k_?h_10mm.mgh 2>/dev/null | wc -l)

if [ "$Nsurf" -lt 8 ]; then
for hemisphere in lh rh; do
    HEMI=$(echo "${hemisphere/h/}" | tr [:lower:] [:upper:])
    Info "Mapping volumetric timeseries to native surface ${hemisphere}"
    vol2surfTS="${rsfmri_surf}/${idBIDS}"_rsfmri_space-fsnative_${hemisphere}.mgh
    if [[ ! -f "$vol2surfTS" ]] ; then

          # Map the non high-passed volumetric timeseries to the surface so we can compute tSNR
          Do_cmd mri_vol2surf \
              --mov "$singleecho" \
              --reg "$fmri2fs_dat" \
              --projfrac-avg 0.2 0.8 0.1 \
              --trgsubject "$idBIDS" \
              --interp trilinear \
              --hemi "${hemisphere}" \
              --out "${rsfmri_surf}/${idBIDS}"_rsfmri_space-fsnative_"${hemisphere}"_NoHP.mgh

          # Map processed timeseries to surface
          Do_cmd mri_vol2surf \
              --mov "$fmri_processed "\
              --reg "$fmri2fs_dat" \
              --projfrac-avg 0.2 0.8 0.1 \
              --trgsubject "$idBIDS" \
              --interp trilinear \
              --hemi "${hemisphere}" \
              --out "$vol2surfTS"

          if [[ -f "$vol2surfTS" ]] ; then ((Nsteps++)); fi
    else
        Info "Subject ${id} volumetric timeseries have been mapped to ${HEMI} cortical surface"; ((Nsteps++))
    fi

    # Convert native timeseries to gifti
    Do_cmd mri_convert "${rsfmri_surf}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.mgh" "${tmp}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.func.gii"

    # Apply smoothing on native surface
    out_surf_native="${rsfmri_surf}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}_10mm.mgh"
    if [[ ! -f "$out_surf_native" ]] ; then
          if [[ "$smooth" == 1 ]] ; then
            Do_cmd wb_command -metric-smoothing \
                "${dir_freesurfer}/surf/${hemisphere}.midthickness.surf.gii"  \
                "${tmp}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.func.gii" \
                10 \
                "${tmp}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}_10mm.func.gii"
            Do_cmd mri_convert "${tmp}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}_10mm.func.gii" "$out_surf_native"
          else
            Do_cmd mri_surf2surf \
                --hemi "${hemisphere}" \
                --srcsubject "$idBIDS" \
                --sval "${rsfmri_surf}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.mgh" \
                --trgsubject "$idBIDS" \
                --tval "$out_surf_native" \
                --fwhm-trg 10
          fi
    if [[ -f "$out_surf_native" ]] ; then ((Nsteps++)); fi
    else
        Info "Subject ${id} has native timeseries smoothed on ${HEMI} surface"; ((Nsteps++))
    fi

    # Register to fsa5 and smooth
    out_surf_fsa5="${rsfmri_surf}/${idBIDS}_rsfmri_space-fsaverage5_${hemisphere}.mgh"
    if [[ ! -f "$out_surf_fsa5" ]] ; then
         Do_cmd mri_surf2surf \
            --hemi "${hemisphere}" \
            --srcsubject "$idBIDS" \
            --sval "${rsfmri_surf}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.mgh" \
            --trgsubject fsaverage5 \
            --tval "$out_surf_fsa5"
         if [[ -f "$out_surf_fsa5" ]] ; then ((Nsteps++)); fi
    else
         Info "Subject ${id} has timeseries mapped to ${HEMI} fsa5"; ((Nsteps++))
    fi

    out_surf_fsa5_sm="${rsfmri_surf}/${idBIDS}_rsfmri_space-fsaverage5_${hemisphere}_10mm.mgh"
    if [[ ! -f "$out_surf_fsa5_sm" ]] ; then
         Do_cmd mri_surf2surf \
            --hemi "${hemisphere}" \
            --srcsubject "$idBIDS" \
            --sval "${rsfmri_surf}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.mgh" \
            --trgsubject fsaverage5 \
            --tval "$out_surf_fsa5_sm" \
            --fwhm-trg 10
         if [[ -f "$out_surf_fsa5_sm" ]] ; then ((Nsteps++)); fi
    else
         Info "Subject ${id} has smoothed timeseries mapped to ${HEMI} fsa5"; ((Nsteps++))
    fi

    # Register to conte69 and smooth
    out_surf="${rsfmri_surf}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}_10mm.mgh"
    if [[ ! -f "$out_surf" ]] ; then
          # Register to conte69
          Do_cmd wb_command -metric-resample \
              "${tmp}/${idBIDS}_rsfmri_space-fsnative_${hemisphere}.func.gii" \
              "${dir_conte69}/${idBIDS}_${hemisphere}_sphereReg.surf.gii" \
              "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMI}.sphere.32k_fs_LR.surf.gii" \
              ADAP_BARY_AREA \
              "${tmp}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}.func.gii" \
              -area-surfs \
              "${dir_freesurfer}/surf/${hemisphere}.midthickness.surf.gii" \
              "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemisphere}_midthickness.surf.gii"
          # Apply smooth on conte69
          Do_cmd wb_command -metric-smoothing \
              "${util_surface}/fsaverage.${HEMI}.midthickness_orig.32k_fs_LR.surf.gii" \
              "${tmp}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}.func.gii" \
              10 \
              "${tmp}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}_10mm.func.gii"

          Do_cmd mri_convert "${tmp}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}.func.gii" "${rsfmri_surf}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}.mgh"
          Do_cmd mri_convert "${tmp}/${idBIDS}_rsfmri_space-conte69-32k_${hemisphere}_10mm.func.gii" "${out_surf}"
          if [[ -f "$out_surf" ]] ; then ((Nsteps++)); fi
    else
          Info "Subject ${id} has a singleecho fmri2fs ${hemisphere} on conte69-32k 10mm surface"; ((Nsteps++))
    fi
done
else
  Info "Subject ${id} has a singleecho fmri2fs on all surfaces: native, native_fwhm-10mm, fsa5, fsa5_fwhm-10mm, and conte69fwhm_10mm"; Nsteps=$((Nsteps+10))
fi

#------------------------------------------------------------------------------#
#                           S U B C O R T E X
# Subcortical segmentation (nativepro) to rsfmri space
rsfmri_subcortex="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_subcortical.nii.gz"
timese_subcortex="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_timeseries_subcortical.txt"

if [[ ! -f "$timese_subcortex" ]] ; then
      Info "Getting subcortical timeseries"
      Do_cmd antsApplyTransforms -d 3 -i "$T1_seg_subcortex" -r "$fmri_mean" -n GenericLabel -t ["$mat_rsfmri_affine",1] -t ["$SyN_rsfmri_affine",1] -t "$SyN_rsfmri_Invwarp" -o "$rsfmri_subcortex" -v -u int
      # Extract subcortical timeseries
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i "$fmri_processed" --seg "$rsfmri_subcortex" --exclude 0 --exclude 16 --avgwf "$timese_subcortex"
      if [[ -f "$timese_subcortex" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has rsfmri subcortical time-series"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
#                           C E R E B E L L U M
rsfmri_cerebellum="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_cerebellum.nii.gz"
timese_cerebellum="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_timeseries_cerebellum.txt"
stats_cerebellum="${rsfmri_volum}/${idBIDS}_space-rsfmri_desc-singleecho_cerebellum_roi_stats.txt"

if [[ ! -f "$timese_cerebellum" ]] ; then
      Info "Getting cerebellar timeseries"
      Do_cmd antsApplyTransforms -d 3 -i "$T1_seg_cerebellum" -r "$fmri_mean" -n GenericLabel -t ["$mat_rsfmri_affine",1] -t ["$SyN_rsfmri_affine",1] -t "$SyN_rsfmri_Invwarp" -o "$rsfmri_cerebellum" -v -u int
      # Extract cerebellar timeseries (mean, one ts per segemented structure, exluding nuclei because many are too small for our resolution)
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i "$fmri_processed" --seg "$rsfmri_cerebellum" --exclude 0 --avgwf "$timese_cerebellum"
      3dROIstats -mask "$rsfmri_cerebellum" -nzsum "$rsfmri_cerebellum" > "$stats_cerebellum"
      if [[ -f "$timese_cerebellum" ]] ; then ((Nsteps++)); fi
else
      Info "Subject ${id} has rsfmri cerebellar time-series"; ((Nsteps++))
fi

# -----------------------------------------------------------------------------------------------
# QC: rsfmri processing Input files
QC_proc-rsfmri

#------------------------------------------------------------------------------#
# run post-rsfmri
cleanTS="${rsfmri_surf}/${idBIDS}_rsfmri_space-conte69-32k_desc-timeseries_clean.txt"
if [[ ! -f "$cleanTS" ]] ; then
    Info "Running rsfMRI post processing"
    labelDirectory="${dir_freesurfer}/label/"
    Do_cmd python "$MICAPIPE"/functions/03_FC.py "$idBIDS" "$proc_rsfmri" "$labelDirectory" "$util_parcelations" "$dir_volum" "$performNSR"
    if [[ -f "$cleanTS" ]] ; then ((Nsteps++)); fi
else
    Info "Subject ${id} has post-processed conte69 time-series"; ((Nsteps++))
fi


#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 21 ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "rsfMRI processing and post processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\tSteps completed : $(printf "%02d" "$Nsteps")/21
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/proc_rsfmri_*.txt)"
grep -v "${id}, ${SES/ses-/}, proc_rsfmri" "${out}/micapipe_processed_sub.csv" > tmpfile && mv tmpfile "${out}/micapipe_processed_sub.csv"
echo "${id}, ${SES/ses-/}, proc_rsfmri, ${status}, $(printf "%02d" "$Nsteps")/21, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" "$eri"), ${PROC}, ${Version}" >> "${out}/micapipe_processed_sub.csv"
cleanup "$tmp" "$nocleanup" "$here"
