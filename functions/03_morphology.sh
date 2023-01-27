#!/bin/bash
#
# Cortical morphology metrics processing:
#
# Generates vertexwise (native, fsa5, and conte69) outputs for:
#   Cortical Thickness
#   Mean Curvature
#
# This workflow makes use of freesurfer outputs and custom python scripts
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micapipe/tree/master/parcellations
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
PROC=$8
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check dependencies Status: POST_STRUCTURAL
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfaceProc ${post_struct_json} | awk -F '"' '{print $4}')
set_surface_directory "${recon}"

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-morphology.json"
micapipe_check_json_status "${module_json}" "morphology"

#------------------------------------------------------------------------------#
Title "Cortical morphology analysis\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Saving temporal dir: $nocleanup"

# Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_post-morpho_${idBIDS}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${proc_struct}/surf/morphology"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_subjsurf}/surf"

#------------------------------------------------------------------------------#
### Cortical Thickness ###
# Register to fsa5 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-fsaverage5_desc-rh_thickness_10mm.mgh" ]]; then
    for hemi in lh rh; do  ((N++))
        # Convert native file to mgh and save in output directory
        Do_cmd mri_convert "${dataDir}/${hemi}.thickness ${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_thickness.mgh"

        Do_cmd mri_surf2surf --hemi "$hemi" \
            --srcsubject "$idBIDS" \
            --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_thickness.mgh" \
            --trgsubject fsaverage5 \
            --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_thickness.mgh"

        Do_cmd mri_surf2surf --hemi "$hemi" \
            --fwhm-trg 10 \
            --srcsubject "$idBIDS" \
            --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_thickness.mgh" \
            --trgsubject fsaverage5 \
            --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_thickness_10mm.mgh"
            if [[ -f "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_thickness_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${id} cortical thickness is registered to fsa5"; Nsteps=$((Nsteps+2)); N=$((N+2))
fi

# Register to conte69 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-conte69-32k_desc-rh_thickness_10mm.mgh" ]]; then
    for hemi in lh rh; do ((N++))
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

        Do_cmd mri_convert "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_thickness.mgh" "${tmp}/${hemi}_thickness.func.gii"

        Do_cmd wb_command -metric-resample \
            "${tmp}/${hemi}_thickness.func.gii" \
            "${dir_conte69}/${idBIDS}_${hemi}_space-fsnative_desc-sphere.surf.gii" \
            "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
            ADAP_BARY_AREA \
            "${tmp}/${hemi}_thickness_c69-32k.func.gii" \
            -area-surfs \
            "${dir_subjsurf}/surf/${hemi}.midthickness.surf.gii" \
            "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii"

        Do_cmd mri_convert "${tmp}/${hemi}_thickness_c69-32k.func.gii" "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_thickness.mgh"

        # Smoothing
        Do_cmd wb_command -metric-smoothing \
            "${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii" \
            "${tmp}/${hemi}_thickness_c69-32k.func.gii" \
            10 \
            "${tmp}/${hemi}_thickness_10mm_c69-32k.func.gii"

        Do_cmd mri_convert "${tmp}/${hemi}_thickness_10mm_c69-32k.func.gii" "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_thickness_10mm.mgh"
        if [[ -f "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_thickness_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${idBIDS} cortical thickness is registered to conte69"; Nsteps=$((Nsteps+2)); N=$((N+2))
fi


#------------------------------------------------------------------------------#
### Curvature ###

# Register to fsa5 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-fsaverage5_desc-rh_curvature_10mm.mgh" ]]; then
    for hemi in lh rh; do ((N++))
        # Convert native file to mgh and save in output directory
        Do_cmd mri_convert "${dataDir}/${hemi}.curv ${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_curvature.mgh"

        Do_cmd mri_surf2surf --hemi "${hemi}" \
            --srcsubject "$idBIDS" \
            --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_curvature.mgh" \
            --trgsubject fsaverage5 \
            --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_curvature.mgh"

        Do_cmd mri_surf2surf --hemi "${hemi}" \
            --fwhm-trg 10 \
            --srcsubject "$idBIDS" \
            --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_curvature.mgh" \
            --trgsubject fsaverage5 \
            --trgsurfval "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_curvature_10mm.mgh"
        if [[ -f "${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_curvature_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${id} curvature is registered to fsa5"; Nsteps=$((Nsteps+2)); N=$((N+2))
fi

# Register to conte69 and apply 10mm smooth
if [[ ! -f "${outDir}/${idBIDS}_space-conte69-32k_desc-rh_curvature_10mm.mgh" ]]; then
    for hemi in lh rh; do ((N++))
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo "$hemisphere" | tr [:lower:] [:upper:])

        Do_cmd mri_convert "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_curvature.mgh" "${tmp}/${hemi}_curv.func.gii"

        Do_cmd wb_command -metric-resample \
            "${tmp}/${hemi}_curv.func.gii" \
            "${dir_conte69}/${idBIDS}_${hemi}_space-fsnative_desc-sphere.surf.gii" \
            "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
            ADAP_BARY_AREA \
            "${tmp}/${hemi}_curv_c69-32k.func.gii" \
            -area-surfs \
            "${dir_subjsurf}/surf/${hemi}.midthickness.surf.gii" \
            "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii"

        Do_cmd mri_convert "${tmp}/${hemi}_curv_c69-32k.func.gii" "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_curvature.mgh"

        # Smoothing
        Do_cmd wb_command -metric-smoothing \
            "${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii" \
            "${tmp}/${hemi}_curv_c69-32k.func.gii" \
            10 \
            "${tmp}/${hemi}_curv_10mm_c69-32k.func.gii"

        Do_cmd mri_convert "${tmp}/${hemi}_curv_10mm_c69-32k.func.gii" "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_curvature_10mm.mgh"
        if [[ -f "${outDir}/${idBIDS}_space-conte69-32k_desc-${hemi}_curvature_10mm.mgh" ]]; then ((Nsteps++)); fi
    done
else
    Info "Subject ${id} curvature is registered to conte69"; Nsteps=$((Nsteps+2)); N=$((N+2))
fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status morphology
micapipe_procStatus "${id}" "${SES/ses-/}" "morphology" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "morphology" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
