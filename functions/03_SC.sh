#!/bin/bash
#
# DWI POST structural TRACTOGRAPHY processing with bash:
#
# POST processing workflow for diffusion MRI TRACTOGRAPHY.
#
# This workflow makes use of MRtrix3
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
tracts=$8
autoTract=$9
keep_tck=${10}
dwi_str=${11}
weighted_SC=${12}
tck_file=${13}
PROC=${14}
filter=SIFT2
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Update path for multiple acquisitions processing
if [[ "${dwi_str}" != "DEFAULT" ]]; then
  dwi_str="acq-${dwi_str/acq-/}"
  dwi_str_="_${dwi_str}"
  export proc_dwi=$subject_dir/dwi/"${dwi_str}"
  export dwi_cnntm=$proc_dwi/connectomes
  export autoTract_dir=$proc_dwi/auto_tract
else
  dwi_str=""; dwi_str_=""
fi

# Check inputs: DWI post TRACTOGRAPHY
fod_wmN="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
dwi_5tt="${proc_dwi}/${idBIDS}_space-dwi_desc-5tt.nii.gz"
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0.nii.gz"
dwi_mask="${proc_dwi}/${idBIDS}_space-dwi_desc-brain_mask.nii.gz"
str_dwi_affine="${dir_warp}/${idBIDS}_space-dwi_from-dwi${dwi_str_}_to-nativepro_mode-image_desc-affine_"
mat_dwi_affine="${str_dwi_affine}0GenericAffine.mat"
dwi_SyN_str="${dir_warp}/${idBIDS}_space-dwi_from-dwi${dwi_str_}_to-dwi_mode-image_desc-SyN_"
dwi_SyN_warp="${dwi_SyN_str}1Warp.nii.gz"
dwi_SyN_affine="${dwi_SyN_str}0GenericAffine.mat"
dti_FA="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.nii.gz"
lut_sc="${util_lut}/lut_subcortical-cerebellum_mics.csv"
# from proc_structural
dwi_cere="${proc_dwi}/${idBIDS}_space-dwi_atlas-cerebellum.nii.gz"
dwi_subc="${proc_dwi}/${idBIDS}_space-dwi_atlas-subcortical.nii.gz"
# TDI output
tdi="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tdi.nii.gz"

# Check inputs
micapipe_check_dependency "proc_structural" "${dir_QC}/${idBIDS}_module-proc_structural.json"
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"
micapipe_check_dependency "proc_dwi" "${dir_QC}/${idBIDS}_module-proc_dwi${dwi_str_}.json"
if [ "${weighted_SC}" != "FALSE" ]; then
  if [ ! -f "$weighted_SC" ]; then Error "The provided weighted map does NOT exist: ${weighted_SC}"; exit 1; fi
fi
if [ "${tck_file}" != "FALSE" ]; then
  if [ ! -f "$tck_file" ]; then Error "The provided tractography (tck) does NOT exist: ${tck_file}"; exit 1; fi
fi
if [ ! -f "$dwi_SyN_affine" ]; then Warning "Subject $id doesn't have an SyN registration, only AFFINE will be apply"; regAffine="TRUE"; else regAffine="FALSE"; fi

# -----------------------------------------------------------------------------------------------
# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-SC-${tracts}${dwi_str_}.json"
micapipe_check_json_status "${module_json}" "SC"

#------------------------------------------------------------------------------#
Title "Tractography and structural connectomes\n\t\tmicapipe $Version, $PROC"
micapipe_software
Note "Number of streamlines:" "${tracts}"
Note "Auto-tractograms     :" "${autoTract}"
Note "Saving tractography  :" "${keep_tck}"
Note "Saving temporal dir  :" "${nocleanup}"
Note "MRtrix will use      :" "${threads} threads"
Note "DWI acquisition      :" "${dwi_str}"
Note "Weighted SC          :" "${weighted_SC}"
Note "tck file             :" "${tck_file}"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_post-dwi_${id}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Create Connectomes directory for the outpust
[[ ! -d "$dwi_cnntm" ]] && Do_cmd mkdir -p "$dwi_cnntm" && chmod -R 770 "$dwi_cnntm"
Do_cmd cd "$tmp"

# -----------------------------------------------------------------------------------------------
# Create fsLR-5k rois volume
for HEMI in R L; do
  fsLR5k_indx="${util_surface}/fsLR-5k.${HEMI}.indices.shape.gii"
  fsLR5k_midt="${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsLR-5k_label-midthickness.surf.gii"
  fsLR5k_rois="${tmp}/${idBIDS}_fsLR-5k_hemi-${HEMI}_rois.nii.gz"
  Do_cmd wb_command -metric-to-volume-mapping "${fsLR5k_indx}" "${fsLR5k_midt}" "${T1nativepro}" "${fsLR5k_rois}" \
      -ribbon-constrained "${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsLR-5k_label-white.surf.gii" "${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsLR-5k_label-pial.surf.gii" \
      -greedy -voxel-subdiv 1
done
seg_fsLR5k="${dir_volum}/${idBIDS}_space-nativepro_T1w_surf-fsLR-5k.nii.gz"

# Threshold overlaping ROIs
Do_cmd fslmaths "${tmp}/${idBIDS}_fsLR-5k_hemi-L_rois.nii.gz" -bin "${tmp}/${idBIDS}_fsLR-5k_hemi-L_bin.nii.gz"
Do_cmd fslmaths "${tmp}/${idBIDS}_fsLR-5k_hemi-R_rois.nii.gz" -bin "${tmp}/${idBIDS}_fsLR-5k_hemi-R_bin.nii.gz"
Do_cmd fslmaths "${tmp}/${idBIDS}_fsLR-5k_hemi-L_bin.nii.gz" -add "${tmp}/${idBIDS}_fsLR-5k_hemi-R_bin.nii.gz" -thr 2 -binv "${tmp}/fsLR-5k_rois_bin.nii.gz"
Do_cmd fslmaths "${tmp}/${idBIDS}_fsLR-5k_hemi-L_rois.nii.gz" -add "${tmp}/${idBIDS}_fsLR-5k_hemi-R_rois.nii.gz" -mul "${tmp}/fsLR-5k_rois_bin.nii.gz" "${seg_fsLR5k}"

# -----------------------------------------------------------------------------------------------
# Prepare the segmentatons
parcellations=($(find "${dir_volum}" -name "*atlas*" ! -name "*cerebellum*" ! -name "*subcortical*"))
# Transformations from T1nativepro to DWI
if [[ ${regAffine}  == "FALSE" ]]; then
    trans_T12dwi="-t ${dwi_SyN_warp} -t ${dwi_SyN_affine} -t [${mat_dwi_affine},1]"
elif [[ ${regAffine}  == "TRUE" ]]; then
    trans_T12dwi="-t [${mat_dwi_affine},1]"
fi
# get wmFOD for registration
fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
Do_cmd mrconvert -coord 3 0 "$fod_wmN" "$fod"

# -----------------------------------------------------------------------------------------------
# Generate probabilistic tracts
weights=${tmp}/SIFT2_${tracts}.txt
if [ ${tck_file} != "FALSE" ]; then
  tck=${tck_file}
else
  tck="${tmp}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tractography.tck"
  Info "Building the ${tracts} streamlines connectome!!!"
  export tckjson="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tractography.json"
  Do_cmd tckgen -nthreads "$threads" \
          "$fod_wmN" \
          "$tck" \
          -act "$dwi_5tt" \
          -crop_at_gmwmi \
          -backtrack \
          -seed_dynamic "$fod_wmN" \
          -algorithm iFOD2 \
          -step 0.5 \
          -angle 22.5 \
          -cutoff 0.06 \
          -maxlength 400 \
          -minlength 10 \
          -select "$tracts"

  # Exit if tractography fails
  if [ ! -f "$tck" ]; then Error "Tractogram failed, check the logs: $(ls -Art "$dir_logs"/post-dwi_*.txt | tail -1)"; exit; fi

  # json file of tractogram
  tck_json iFOD2 0.5 22.5 0.06 400 10 seed_dynamic "$tck"
fi

# SIFT2
Do_cmd tcksift2 -nthreads "$threads" "$tck" "$fod_wmN" "$weights"

# TDI for QC
Info "Creating a Track Density Image (tdi) of the $tracts connectome for QC"
Do_cmd tckmap -template "${fod}" -dec -nthreads "$threads" "$tck" "$tdi" -force

# Map the weighted image to the whole brain tractography
if [ "${weighted_SC}" != "FALSE" ]; then Do_cmd tcksample "${tck}" "${weighted_SC}" "${tmp}"/mean_map_per_streamline.txt -stat_tck mean -nthreads "$threads"; fi

# -----------------------------------------------------------------------------------------------
# Build the Connectomes
function build_connectomes(){
	nodes=$1
	sc_file=$2
	# Build the weighted connectomes
    Do_cmd tck2connectome -nthreads "$threads" \
    	"${tck}" "${nodes}" "${sc_file}-connectome.txt" \
        -tck_weights_in "$weights" -quiet -force
    Do_cmd python "${MICAPIPE}"/functions/connectome_slicer.py --conn="${sc_file}-connectome.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"

    # Calculate the edge lenghts
    Do_cmd tck2connectome -nthreads "$threads" \
        "${tck}" "${nodes}" "${sc_file}-edgeLengths.txt" \
        -tck_weights_in "$weights" -scale_length -stat_edge mean -quiet -force
    Do_cmd python "${MICAPIPE}"/functions/connectome_slicer.py --conn="${sc_file}-edgeLengths.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"

    # Weighted connectome with a NIFTI map
	if [ ${weighted_SC} != "FALSE" ]; then
		Do_cmd tck2connectome -nthreads "$threads" \
			"${tck}" "${nodes}" "${sc_file}-weighted_connectome.txt" \
			-tck_weights_in "$weights" -scale_file "${tmp}"/mean_map_per_streamline.txt -stat_edge mean -quiet  -force
      Do_cmd python "${MICAPIPE}"/functions/connectome_slicer.py --conn="${sc_file}-weighted_connectome.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"
	fi
}

# Vertex-wise connectome
fs5k_cortex="${tmp}/${id}_fsLR5k-cor_dwi.nii.gz" # Segmentation in dwi space
lut="${util_lut}/lut_fsLR-5k_mics.csv"
Do_cmd antsApplyTransforms -d 3 -e 3 -i "$seg_fsLR5k" -r "${fod}" -n GenericLabel "$trans_T12dwi" -o "$fs5k_cortex" -u int
build_connectomes "$fs5k_cortex" "${dwi_cnntm}/${idBIDS}_surf-fsLR-5k_desc-iFOD2-${tracts}-${filter}_full"

# Create a connectome per parcellation
for seg in "${parcellations[@]}"; do
    parc_name=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    connectome_str="${dwi_cnntm}/${idBIDS}_space-dwi_atlas-${parc_name}_desc-iFOD2-${tracts}-${filter}"
    lut="${util_lut}/lut_${parc_name}_mics.csv"
    dwi_cortex="${tmp}/${id}_${parc_name}-cor_dwi.nii.gz" # Segmentation in dwi space
    dwi_cortexSub="${tmp}/${id}_${parc_name}-sub_dwi.nii.gz"
    dwi_all="${tmp}/${id}_${parc_name}-full_dwi.nii.gz"
    # -----------------------------------------------------------------------------------------------
    # Build the Full connectome (Cortical-Subcortical-Cerebellar)
    if [[ ! -f "${connectome_str}_full-connectome.shape.gii" ]]; then ((N++))
        Info "Building $parc_name cortical-subcortical-cerebellum connectome"
        # Take parcellation into DWI space
        Do_cmd antsApplyTransforms -d 3 -e 3 -i "$seg" -r "${fod}" -n GenericLabel "$trans_T12dwi" -o "$dwi_cortex" -v -u int
        # Remove the medial wall
        for i in 1000 2000; do Do_cmd fslmaths "$dwi_cortex" -thr "$i" -uthr "$i" -binv -mul "$dwi_cortex" "$dwi_cortex"; done
        Do_cmd fslmaths "$dwi_cortex" -binv -mul "$dwi_subc" -add "$dwi_cortex" "$dwi_cortexSub" -odt int # added the subcortical parcellation
        Do_cmd fslmaths "$dwi_cortex" -binv -mul "$dwi_cere" -add "$dwi_cortexSub" "$dwi_all" -odt int # added the cerebellar parcellation
        # Build the Cortical-Subcortical-Cerebellum connectomes
        build_connectomes "$dwi_all" "${connectome_str}_full"
        if [[ -f "${connectome_str}_full-connectome.shape.gii" ]]; then ((Nsteps++)); fi
    else
        ((N++)); ((Nsteps++))
    fi
done

# Change connectome permissions
chmod 770 -R "$dwi_cnntm"/* 2>/dev/null

# -----------------------------------------------------------------------------------------------
# Compute Auto-Tractography
if [ "$autoTract" == "TRUE" ]; then
    Info "Running Auto-tract"
    autoTract_dir="$proc_dwi"/auto_tract
    [[ ! -d "$autoTract_dir" ]] && Do_cmd mkdir -p "$autoTract_dir"
    echo -e "\033[38;5;118m\nCOMMAND -->  \033[38;5;122m03_auto_tracts.sh -tck ${tck} -outbase ${autoTract_dir}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}-${filter} -mask ${dwi_mask} -fa ${dti_FA} -weights ${weights} -tmpDir ${tmp} -keep_tmp  \033[0m"
    "$MICAPIPE"/functions/03_auto_tracts.sh -tck "$tck" -outbase "${autoTract_dir}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}-${filter}" -mask "$dwi_mask" -fa "$dti_FA" -weights "$weights" -tmpDir "$tmp" -keep_tmp
fi

# -----------------------------------------------------------------------------------------------
# save the tractogram and the SIFT2 weights
if [ "$keep_tck" == "TRUE" ]; then Do_cmd mv "$tck" "$proc_dwi"; Do_cmd mv "$weights" "${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tractography_weights.txt"; fi
#if [[ $nocleanup == "FALSE" ]]; then rm -Rf "$seg_fsLR5k"; else mv "${seg_fsLR5k}" ${tmp}; fi

# -----------------------------------------------------------------------------------------------
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status "SC-${tracts}${dwi_str_}"
micapipe_procStatus "${id}" "${SES/ses-/}" "SC-${tracts}${dwi_str_}" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "SC-${tracts}${dwi_str_}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
