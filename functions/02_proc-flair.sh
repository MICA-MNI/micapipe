#!/bin/bash
#
# T2-FLAIR processing:
#
# Performs rescaling and normalization of T2/FLAIR image and registers to nativepro
#
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
flairScanStr=$8
PROC=${9}
here=$(pwd)
export OMP_NUM_THREADS="${threads}"

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE"/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check dependencies Status: POST_STRUCTURAL
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"

#------------------------------------------------------------------------------#
# Check inputs: T2-FLAIR
if [[ "$flairScanStr" == DEFAULT ]]; then
    # T2/FLAIR scan
    N_flairScan=${#bids_flair[@]}
    if [ "$N_flairScan" -gt 1 ]; then
        if [[ "${flairScanStr}" == "DEFAULT" ]]; then
            Error "Multiple FLAIR runs found in BIDS rawdata directory! Please specify which one you want to process"; exit;
        fi
    else
        flairScan=${bids_flair[0]}
    fi
else
    Info "Using string to set the FLAIR scan: ${flairScanStr}"
    flairScan=$(realpath "${subject_bids}/anat/${idBIDS}_${flairScanStr}.nii"* 2>/dev/null)
fi
if [ ! -f "$flairScan" ]; then Error "T2-flair not found for Subject $idBIDS : ${flairScan}"; exit; fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_flair.json"
micapipe_check_json_status "${module_json}" "proc_flair"

#------------------------------------------------------------------------------#
Title "T2-FLAIR intensities\n\t\tmicapipe $Version, $PROC"
micapipe_software
Info "Module inputs:"
Note "wb_command threads  : " "${OMP_NUM_THREADS}"
Note "threads             : " "${threads}"
Note "Saving temporary dir: " "${nocleanup}"
Note "tmpDir              : " "${tmpDir}"
Note "flairScanStr        : " "${flairScanStr}"
Note "Processing          : " "${PROC}"

# Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_flair_${idBIDS}"
umask 000; mkdir -m 777 -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${subject_dir}/maps"

# json file
flair_json="${outDir}/${idBIDS}_space-nativepro_map-flair.json"

#------------------------------------------------------------------------------#
# Functions
function Image_threshold() {
  img_in=$1
  thrlo=$2
  thrhi=$3
  img_out=$4

  # Create upper and lower thresholded images
  rnd_str=${RANDOM}
  threshlo="${tmp}/${rnd_str}_threshlo.nii.gz"
  threshhi="${tmp}/${rnd_str}_threshhi.nii.gz"

  ThresholdImage 3 "${img_in}" "${threshlo}" threshlo "${thrlo}"
  ThresholdImage 3 "${img_in}" "${threshhi}" threshhi "${thrhi}"

  # Substract the thresholded images to get the unique thresholded img
  ImageMath 3 "${img_out}" - "${threshlo}" "${threshhi}"

  # Cleanup
  rm "${threshlo}" "${threshhi}"

}

# The the mode of the GM, WM and whole brain
function get_mode() {
  # This functions uses mrhistogram to calculate the mode with awk.
  # Using 1000 bins gets the same results as python
  MRI_img=$1
  MRI_mask=$2
  bins=$3
  # Create histogram file
  hist=${tmp}/hist_masked.txt
  mrhistogram -bins "${bins}" -mask "${MRI_mask}" -ignorezero -nthreads "${threads}" "${MRI_img}" "${hist}"
  # get the index with the max frecuency
  max_val=$(awk -F ',' 'NR==3 { m=$3; p=1; for(i=4;i<=NF;i++) { if ($i>m) { m=$i; p=i-2 } } printf "%d ",p }' "${hist}")
  # bash array of each intensity bin
  intensities=($(awk 'NR==2' "${hist}" | tr -s ',' ' '))
  # get the intensity of the maximun frecuency (aka mode)
  mode=${intensities[$((max_val-1))]}
  # remove tmp file
  rm "${hist}"
  # Print the mode
  echo "${mode}" | awk -F"e" 'BEGIN{OFMT="%10.10f"} {print $1 * (10 ^ $2)}'
}

#------------------------------------------------------------------------------#
# preproc FLAIR in nativepro space
flair_nativepro="${outDir}/${idBIDS}_space-nativepro_map-flair.nii.gz"
flair_preproc="${tmp}/${idBIDS}_flair_norm_brain.nii.gz"
flair_synthseg="${tmp}/${idBIDS}_flair_synthseg.nii.gz"
flair_resample="${tmp}/${idBIDS}_flair_synthseg_resampled.nii.gz"

if [[ ! -f "$flair_nativepro" ]]; then ((N++))
    #------------------------------------------------------------------------------#
    # 1 | FLAIR segmentation
    Do_cmd mri_synthseg --i "${flairScan}" --o "${flair_synthseg}" --robust --threads "$threads" --cpu --resample "${flair_resample}"

    # synthseg_resampled is only generated if the voxel sizes don't match; in this case, you will need to do an affine registration back to FLAIR native space
    if [ -f "${tmp}/${idBIDS}_flair_synthseg_resampled.nii.gz" ]; then
        # Affine registration from flair_synthseg_resampled to FLAIR native space
        flair_resample2orig="${tmp}/${idBIDS}_flair_synthresample_to_orig_"
        Do_cmd antsRegistrationSyN.sh -d 3 -f "${flairScan}" -m "${flair_resample}" -o "${flair_resample2orig}" -t a -n "$threads" -p d
        # Apply transformation to flair_synthseg
        flair_synthseg_orig="${tmp}/${idBIDS}_flair_synthseg_orig.nii.gz"
        Do_cmd antsApplyTransforms -d 3 -i "${flair_synthseg}" -r "${flairScan}" -t "${flair_resample2orig}"0GenericAffine.mat -o "${flair_synthseg_orig}" -v -u float -n GenericLabel
        # flair_synthseg in original FLAIR space
        flair_synthseg="${flair_synthseg_orig}"
    fi

    # 1.A | GM gray matter mask (lh=3, rh=42)
    flair_mask_gm_lh="${tmp}/${idBIDS}_hemi-L_label-gm_flair_mask.nii.gz"
    flair_mask_gm_rh="${tmp}/${idBIDS}_hemi-R_label-gm_flair_mask.nii.gz"
    flair_mask_gm="${tmp}/${idBIDS}_label-gm_flair_mask.nii.gz"
    Image_threshold "${flair_synthseg}" 3 2 "${flair_mask_gm_lh}"
    Image_threshold "${flair_synthseg}" 42 41 "${flair_mask_gm_rh}"
    ImageMath 3 "${flair_mask_gm}" + "${flair_mask_gm_lh}" "${flair_mask_gm_rh}"

    # 1.B | WM white matter mask (lh=2, rh=41)
    flair_mask_wm_lh="${tmp}/${idBIDS}_hemi-L_label-wm_flair_mask.nii.gz"
    flair_mask_wm_rh="${tmp}/${idBIDS}_hemi-R_label-wm_flair_mask.nii.gz"
    flair_mask_wm="${tmp}/${idBIDS}_label-wm_flair_mask.nii.gz"
    Image_threshold "${flair_synthseg}" 2 1 "${flair_mask_wm_lh}"
    Image_threshold "${flair_synthseg}" 41 40 "${flair_mask_wm_rh}"
    ImageMath 3 "${flair_mask_wm}" + "${flair_mask_wm_lh}" "${flair_mask_wm_rh}"

    #------------------------------------------------------------------------------#
    # 2 | Bias field correction weighted by white matter
    flair_N4="${tmp}/${idBIDS}_flairN4.nii.gz"
    Do_cmd N4BiasFieldCorrection -r -d 3 -w "${flair_mask_wm}" -i "${flairScan}" -o "${flair_N4}"

    #------------------------------------------------------------------------------#
    # 3 | Brain mask
    flair_mask="${tmp}/${idBIDS}_label-brain_flairN4_mask.nii.gz"
    Image_threshold "${flair_synthseg}" 60 1 "${flair_mask}"


    #------------------------------------------------------------------------------#
    # 4 | Get the mode for each tissue
    mode_gm=$(get_mode "${flair_N4}" "${flair_mask_gm}" 1000)
    mode_wm=$(get_mode "${flair_N4}" "${flair_mask_wm}" 1000)
    mode_brain=$(get_mode "${flair_N4}" "${flair_mask}" 1000)

    Note "mode_gm    :" "${mode_gm}"
    Note "mode_wm    :" "${mode_wm}"
    Note "mode_brain :" "${mode_brain}"

    #------------------------------------------------------------------------------#
    # 5 | Normalize intensities by peak of WM (mode).
    # This normalization will center the peak of the WM mode intensity at ZERO.
    # Mean mode between GM and WM | BG=(GM_mode+WM_mode)/2.0
    BG=$(echo "(${mode_gm}+${mode_wm})/2.0" | bc -l)
    # mode difference | mode_diff = np.abs(BG - WM_mode)
    mode_diff=$(echo "${BG}-${mode_wm}" | bc); mode_diff=$(echo ${mode_diff#-})
    # Normalize array | norm_wm = 100.0 * (array - WM_mode)/(mode_diff)
    flair_norm="${tmp}/${idBIDS}_flair_norm.nii.gz"
    Do_cmd mrcalc "${flair_N4}" "${mode_wm}" -subtract "${mode_diff}" -div 100 -mul "${flair_norm}"

    #------------------------------------------------------------------------------#
    # 6 | Mask only the brain of the normalized data
    Do_cmd ImageMath 3 "${flair_preproc}" m "${flair_mask}" "${flair_norm}"

    ((Nsteps++))
else
    Info "Subject ${id} T2-FLAIR has been processed"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
### FLAIR registration to nativepro ###
if [[ ! -f "$flair_nativepro" ]]; then ((N++))
    # Register nativepro and flair
    str_flair_affine="${dir_warp}/${idBIDS}_from-flair_to-nativepro_mode-image_desc-affine_"
    Info "Running label based affine registrations"
    flair_synth="${tmp}/flair_synthsegGM.nii.gz"
    T1_synth="${tmp}/T1w_synthsegGM.nii.gz"
    Do_cmd mri_synthseg --i "${T1nativepro}" --o "${tmp}/T1w_synthseg.nii.gz" --robust --threads "$threads" --cpu
    fslmaths "${tmp}/T1w_synthseg.nii.gz" -uthr 42 -thr 42 -bin -mul -39 -add "${tmp}/T1w_synthseg.nii.gz" "${T1_synth}"
    fslmaths "${flair_synthseg}" -uthr 42 -thr 42 -bin -mul -39 -add "${flair_synthseg}" "${flair_synth}"

    # Affine from func to t1-nativepro
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1_synth" -m "$flair_synth" -o "$str_flair_affine" -t a -n "$threads" -p d

    # Apply transformations
    Do_cmd antsApplyTransforms -d 3 -i "$flair_preproc" -r "$T1nativepro_brain" -t "$str_flair_affine"0GenericAffine.mat -o "$flair_nativepro" -v -u float
    ((Nsteps++))
else
    Info "Subject ${id} T2-FLAIR is registered to nativepro"; ((Nsteps++)); ((N++))
fi

# Write json file
json_nativepro_flair "$flair_nativepro" \
    "antsApplyTransforms -d 3 -i ${flair_preproc} -r ${T1nativepro_brain} -t ${str_flair_affine}0GenericAffine.mat -o ${flair_nativepro} -v -u float" \
    "$flair_json"

#------------------------------------------------------------------------------#
# Map to surface
Nmorph=$(ls "${dir_maps}/"*_flair.func.gii 2>/dev/null | wc -l)
if [[ "$Nmorph" -lt 16 ]]; then ((N++))
    Info "Mapping flair to fsLR-32k, fsLR-5k and fsaverage5"
    for HEMI in L R; do
        for label in midthickness white; do
            surf_fsnative="${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-${label}.surf.gii"
            # MAPPING metric to surfaces
            map_to-surfaces "${flair_nativepro}" "${surf_fsnative}" "${dir_maps}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-${label}_flair.func.gii" "${HEMI}" "${label}_flair" "${dir_maps}"
        done
    done
    Nmorph=$(ls "${dir_maps}/"*_flair.func.gii 2>/dev/null | wc -l)
    if [[ "$Nmorph" -eq 16 ]]; then ((Nsteps++)); fi
else
    Info "Subject ${idBIDS} has flair mapped to surfaces"; ((Nsteps++)); ((N++))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
module_name="proc_flair"
micapipe_completition_status "${module_name}"
micapipe_procStatus "${id}" "${SES/ses-/}" "${module_name}" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "${module_name}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
