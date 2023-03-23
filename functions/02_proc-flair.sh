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
PROC=$9
here=$(pwd)


#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source ${MICAPIPE}/functions/init.sh;
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check dependencies Status: POST_STRUCTURAL
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"

Info "Inputs of T2/FLAIR:"
Note "tmpDir     : " "$tmpDir"
Note "flairScanStr : " "$flairScanStr"
Note "Processing : " "$PROC"

#------------------------------------------------------------------------------#
# Check inputs: T2-FLAIR
if [ ! -f "$bids_flair" ]; then Error "T2-flair not found for Subject $id : ${subject_bids}/anat/${idBIDS}*FLAIR.nii.gz"; exit; fi

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
    Info "Using user provided FLAIR scan string: ${flairScanStr}"
    flairScan=;$(ls "${subject_bids}/anat/${idBIDS}_${flairScanStr}".nii* 2>/dev/null);
fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-proc_flair.json"
micapipe_check_json_status "${module_json}" "proc_flair"

#------------------------------------------------------------------------------#
Title "T2-FLAIR intensities\n\t\tmicapipe $Version, $PROC"
micapipe_software
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporary dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_flair_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${subject_dir}/maps"

# json file
flair_json="${outDir}/${idBIDS}_space-nativepro_map-flair.json"

#------------------------------------------------------------------------------#
### FLAIR intensity correction ###
# preproc FLAIR in nativepro space
flairNP="${outDir}/${idBIDS}_space-nativepro_map-flair.nii.gz"

# Bias field correction
flair_N4="${tmp}/${idBIDS}_space-flair_desc-flair_N4.nii.gz"
if [[ ! -f "$flairNP" ]]; then
    ((N++))
    Do_cmd N4BiasFieldCorrection -d 3 -i "$bids_flair" -r \
                                -o "$flair_N4"
    ((Nsteps++))
else
    Info "Subject ${id} T2-FLAIR is N4 bias corrected"; ((Nsteps++)); ((N++))
fi

# Clamp and rescale intensities
flair_clamp="${tmp}/${idBIDS}_space-flair_desc-flair_N4_clamp.nii.gz"
flair_rescale="${tmp}/${idBIDS}_space-flair_desc-flair_N4_rescale.nii.gz"
if [[ ! -f "$flairNP" ]]; then
    ((N++))

    # Clamp intensities
    Do_cmd ImageMath 3 "$flair_clamp" TruncateImageIntensity "$flair_N4" 0.01 0.99 75

    # Rescale intensity [0,100]
    Do_cmd ImageMath 3 "$flair_rescale" RescaleImage "$flair_clamp" 0 100

    ((Nsteps++))
else
    Info "Subject ${id} T2-FLAIR is intensity corrected"; ((Nsteps++)); ((N++))
fi

# Normalize intensities by GM/WM interface, uses 5ttgen
flair_preproc="${tmp}/${idBIDS}_space-flair_desc-flair_preproc.nii.gz"
if [[ ! -f "$flairNP" ]]; then ((N++))
    # Get gm/wm interface mask
    t1_gmwmi="${tmp}/${idBIDS}_space-nativepro_desc-gmwmi-mask.nii.gz"
    t1_5tt="${proc_struct}/${idBIDS}_space-nativepro_T1w_5tt.nii.gz"
    if [[ ! -f "$t1_gmwmi" ]]; then ((N++))
        Info "Calculating Gray matter White matter interface mask"
        Do_cmd 5tt2gmwmi "$t1_5tt" "$t1_gmwmi"
        ((Nsteps++))
    else
        Info "Subject ${id} has Gray matter White matter interface mask"; ((Nsteps++)); ((N++))
    fi

    # Register nativepro and flair
    str_flair_affine="${dir_warp}/${idBIDS}_from-flair_to-nativepro_mode-image_desc-affine_"
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1nativepro_brain" -m "$flair_rescale" -o "$str_flair_affine" -t a -n "$threads" -p d
    t1_gmwmi_in_flair="${tmp}/${idBIDS}_space-flair_desc-gmwmi-mask.nii.gz"
    Do_cmd antsApplyTransforms -d 3 -i "$t1_gmwmi" -r "$flair_rescale" -t ["$str_flair_affine"0GenericAffine.mat,1] -o "$t1_gmwmi_in_flair" -v -u float

    # binarize mask
    t1_gmwmi_in_flair_thr="${tmp}/${idBIDS}_space-flair_desc-gmwmi-thr.nii.gz"
    fslmaths "$t1_gmwmi_in_flair" -thr 0.5 -bin "$t1_gmwmi_in_flair_thr"

    # compute mean flair intensity in non-zero voxels
    gmwmi_mean=$(fslstats "$flair_rescale" -M -k "$t1_gmwmi_in_flair_thr")

    # Normalize flair
    fslmaths "$flair_rescale" -div $gmwmi_mean "$flair_preproc"

    ((Nsteps++))
else
    Info "Subject ${id} T2-FLAIR is normalized by GM/WM interface"; ((Nsteps++)); ((N++))
fi


#------------------------------------------------------------------------------#
### FLAIR registration to nativepro ###
if [[ ! -f "$flairNP" ]]; then
    Do_cmd antsApplyTransforms -d 3 -i "$flair_preproc" -r "$T1nativepro_brain" -t "$str_flair_affine"0GenericAffine.mat -o "$flairNP" -v -u float
    ((Nsteps++)); ((N++))
else
    Info "Subject ${id} T2-FLAIR is registered to nativepro"; ((Nsteps++)); ((N++))
fi

# Write json file
json_nativepro_flair "$flairNP" \
    "antsApplyTransforms -d 3 -i ${flair_preproc} -r ${T1nativepro_brain} -t ${str_flair_affine}0GenericAffine.mat -o ${flairNP} -v -u float" \
    "$flair_json"

#------------------------------------------------------------------------------#
# Map to surface
Nmorph=$(ls "${dir_maps}/"*flair*gii 2>/dev/null | wc -l)
if [[ "$Nmorph" -lt 8 ]]; then ((N++))
    Info "Mapping flair to fsLR-32k, fsLR-5k and fsaverage5"
    for HEMI in L R; do
        surf_fsnative="${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsnative_label-midthickness.surf.gii"
        # MAPPING metric to surfaces
        map_to-surfaces "${flairNP}" "${surf_fsnative}" "${dir_maps}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-midthickness_flair.func.gii" "${HEMI}" "midthickness_flair" "${dir_maps}"
    done
    Nmorph=$(ls "${dir_maps}/"*flair*gii 2>/dev/null | wc -l)
    if [[ "$Nmorph" -eq 8 ]]; then ((Nsteps++)); fi
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
micapipe_completition_status ${module_name}
micapipe_procStatus "${id}" "${SES/ses-/}" "${module_name}" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "${module_name}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
