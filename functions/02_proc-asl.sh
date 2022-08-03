#!/bin/bash
#
# ASL processing with bash:
#
# Preprocessing workflow for ASL MRI.
#
# This workflow makes use of FSL, BASIL, ANTs, AFNI
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#
umask 003
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
aslScanStr=$8
m0ScanStr=$9
PROC=${10}
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ]; then
    export MICAPIPE=/host/yeatman/local_raid/ango/micapipe/ #rcruces/git_here/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# mtrix configuration file
mrconf="${HOME}/.mrtrix.conf"
echo "BZeroThreshold: 61" > "$mrconf"

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

Info "Inputs of ASL:"
Note "tmpDir     : " "$tmpDir"
Note "aslScanStr : " "$aslScanStr"
Note "m0ScanStr  : " "$m0ScanStr"
Note "Processing : " "$PROC"

#------------------------------------------------------------------------------#
if [[ "$aslScanStr" == DEFAULT ]]; then
    # ASL scan
    N_aslScan=${#bids_asl[@]}
    if [ "$N_aslScan" -gt 1 ]; then
        if [[ "${aslScanStr}" == "DEFAULT" ]]; then
            Error "Multiple ASL runs found in BIDS rawdata directory!"; exit;
        fi
    else
        aslScan=${bids_asl[0]}
    fi
else
    Info "Using user provided ASL scan string: ${aslScanStr}"
    aslScan=;$(ls "${subject_bids}/perf/${idBIDS}_${aslScanStr}".nii* 2>/dev/null);
fi

if [[ "$m0ScanStr" == DEFAULT ]]; then
    # M0 calibration scan
    N_m0Scan=${#bids_m0[@]}
    if [ "$N_m0Scan" -gt 1 ]; then
        if [[ "${m0ScanStr}" == "DEFAULT" ]]; then
            Error "Multiple M0 runs found in BIDS rawdata directory!"; exit;
        fi
    else
        m0Scan=${bids_m0[0]}
    fi
else
    Info "Using user provided M0 scan string: ${m0ScanStr}"
    m0Scan=$(ls "${subject_bids}/perf/${idBIDS}_${m0ScanStr}".nii* 2>/dev/null);
fi

#------------------------------------------------------------------------------#
Title "Arterial Spin Labeling Imaging processing\n\t\tmicapipe $Version, $PROC"
micapipe_software
Note "Saving temporal dir     : " "$nocleanup"
Note "ANTs and MRtrix will use: " "$threads threads"

#	Timer
aloita=$(date +%s)
Nsteps=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_ASL_${id}"
Do_cmd mkdir -p "$tmp"
[[ ! -d "$dir_QC_png" ]] && Do_cmd mkdir -p "$dir_QC_png" && chmod -R 770 "$dir_QC_png"

# TRAP in case the script fails
trap 'rm $mrconf; cleanup $tmp $nocleanup $here' SIGINT SIGTERM

Do_cmd cd "$tmp"

# Check inputs
if [ ! -f "${aslScan}" ]; then Error "Couldn't find $id ASL scan : \n\t ls ${subject_bids}/perf/"; exit; fi
if [ ! -f "${m0Scan}" ]; then Error "Couldn't find $id m0 calibration scan : \n\t ls ${subject_bids}/perf/"; exit; fi
if [ ! -f "${T1nativepro_brain}" ]; then Error "Subject $id doesn't have T1_nativepro.\n\t\tRun -proc_structural"; exit; fi

# Create output directory
export proc_asl=$subject_dir/perf/
[[ ! -d "$proc_asl" ]] && Do_cmd mkdir -p "$proc_asl" && chmod -R 770 "$proc_asl"

#------------------------------------------------------------------------------#
#------------------------------- ASL processing -------------------------------#
#------------------------------------------------------------------------------#

#---------------------------------- ASL mask ----------------------------------#
asl_str="${idBIDS}_desc-asl"
asl_brain="${asl_str}_brain.nii.gz"
asl_mask="${asl_str}_brain_mask.nii.gz"
aslref="${asl_str}ref.nii.gz"

if [[ ! -f ${proc_asl}/${asl_mask} ]] || [[ ! -f ${proc_asl}/${aslref} ]]; then
    Info "Generating ASL binary mask"

    # Calculate mean ASL volume (reference image)
    Do_cmd fslmaths ${aslScan} -Tmean ${tmp}/${aslref}

    # Creat mask from reference image
    Do_cmd bet ${tmp}/${aslref} ${tmp}/${asl_brain} -B -m -n # ASL brain

    # Move important files to derivatives
    Do_cmd mv ${tmp}/${asl_mask} ${proc_asl}/${asl_mask}
    Do_cmd mv ${tmp}/${aslref} ${proc_asl}/${aslref}

    if [[ -f ${proc_asl}/${asl_mask} ]] && [[ -f ${proc_asl}/${aslref} ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a binary mask of the ASL scan"; ((Nsteps++))
fi

#------------ Preprocess ASL image using oxford_asl (from Basil) -------------#
m0_str="${idBIDS}_desc-m0"
m0_brain="${m0_str}_brain.nii.gz"
m0_mask="${m0_str}_brain_mask.nii.gz"
cbf_str="${idBIDS}_desc-preproc_cbf"
cbf="${cbf_str}.nii.gz"

if [[ ! -f ${proc_asl}/${cbf} ]]; then
    Info "Processing ASL scan with oxford_asl from the BASIL toolbox"

    # Mask ASl and M0 timeseries
    Do_cmd bet ${aslScan} ${tmp}/${asl_brain} -F -B -m # ASL brain
    Do_cmd bet ${m0Scan} ${tmp}/${m0_brain} -F -B -m # M0 brain

    Do_cmd fslmaths ${aslScan} -mul ${tmp}/${asl_mask} ${tmp}/${asl_brain}
    Do_cmd fslmaths ${m0Scan} -mul ${tmp}/${m0_mask} ${tmp}/${m0_brain}

    # Compute absolute CBF map using oxford_asl
    Do_cmd fsl_anat --clobber --nosubcortseg \
                    -i ${T1nativepro} \
                    -o ${tmp}/fsl_anat

    # Change flag to fit with ASL protocol
    Do_cmd oxford_asl -i ${tmp}/${asl_brain} \
                      -o ${tmp}/oxford_asl \
                      --iaf=tc --tis=3.25 --bolus 1.7 --slicedt=0.04571 --casl --wp --mc --pvcorr \
                      -c ${tmp}/${m0_brain} \
                      --tr 10 --fslanat=${tmp}/fsl_anat.anat

    # Move important files to derivatives
    Do_cmd mv ${tmp}/oxford_asl/native_space/perfusion_calib.nii.gz ${proc_asl}/${cbf} # CBF map

    if [[ -f ${proc_asl}/${cbf} ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has preprocessed ASL scan from BASIL"; ((Nsteps++))
fi

#----------------------- Register CBF to nativepro T1w ------------------------#
aslref_brain="${asl_str}ref_brain.nii.gz"
aslref_brain_LPI_orient="${asl_str}ref_brain_LPI.nii.gz"
cbf_LPI_orient="${cbf_str}_LPI.nii.gz"
str_asl_affine="${dir_warp}/${idBIDS}_from-asl_to-nativepro_mode-image_desc-affine_"
mat_asl_affine="${str_asl_affine}0GenericAffine.mat"
cbf_in_nativepro="${idBIDS}_space-nativepro_desc-preproc_cbf.nii.gz"

if [[ ! -f ${proc_asl}/${cbf_in_nativepro} ]]; then
    Info "Registering ASL MRI to nativepro T1w"

    # Mask ASL reference scan
    Do_cmd fslmaths ${proc_asl}/${aslref} -mul ${proc_asl}/${asl_mask} ${tmp}/${aslref_brain}

    # Resample/reorient asl reference/cbf to nativepro T1w orientation (LPI)
    Do_cmd 3dresample -orient LPI -dxyz 0.8 0.8 0.8 -prefix ${tmp}/${aslref_brain_LPI_orient} -inset ${tmp}/${aslref_brain}
    Do_cmd fslreorient2std ${tmp}/${aslref_brain_LPI_orient} ${tmp}/${aslref_brain_LPI_orient}
    Do_cmd 3dresample -orient LPI -dxyz 0.8 0.8 0.8 -prefix ${tmp}/${cbf_LPI_orient} -inset ${proc_asl}/${cbf}
    Do_cmd fslreorient2std ${tmp}/${cbf_LPI_orient} ${tmp}/${cbf_LPI_orient}

    # Generate affine from asl to t1-nativepro
    Do_cmd antsRegistrationSyN.sh -d 3 -f ${T1nativepro_brain} -m ${tmp}/${aslref_brain_LPI_orient} -o ${str_asl_affine} -t a -n "$threads" -p d

    # Register CBF to nativepro T1w
    Do_cmd antsApplyTransforms -d 3 -i ${tmp}/${cbf_LPI_orient} -r ${T1nativepro_brain} -t ${mat_asl_affine} -o ${tmp}/${cbf_in_nativepro} -v -u double -n NearestNeighbor
    Do_cmd fslmaths ${tmp}/${cbf_in_nativepro} -mul ${T1nativepro_mask} ${tmp}/${cbf_in_nativepro}

    # Move important files to derivatives
    Do_cmd mv ${tmp}/${cbf_in_nativepro} ${proc_asl}/${cbf_in_nativepro} # CBF in nativepro T1w space

    if [[ -f ${proc_asl}/${cbf_in_nativepro} ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has CBF and transformation matrix in T1nativepro space"; ((Nsteps++))
fi

#--------------- Register to fsa5/conte69 and apply 10mm smooth ---------------#
asl_fsa5

if
# -----------------------------------------------------------------------------------------------
# QC: Input files
# QC_proc-asl "${dir_QC}/micapipe_QC_proc-asl${asl_str_}.txt"

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
N=3
if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "ASL processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\t\tSteps completed: $(printf "%02d" "$Nsteps")/3
\tStatus          : ${status}
\tCheck logs:
$(ls "${dir_logs}"/proc_asl_*.txt)"
# Print QC stamp
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_asl" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_asl" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
cleanup "$tmp" "$nocleanup" "$here"
