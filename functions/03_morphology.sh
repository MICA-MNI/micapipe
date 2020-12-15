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
PROC=$5
nocleanup=$6
threads=$7
export OMP_NUM_THREADS=$threads

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source ${MICAPIPE}/functions/init.sh;
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables $BIDS $id $out $SES

# Check inputs: Freesurfer space T1
if [ ! -f ${T1freesurfr} ]; then Error "T1 in freesurfer space not found for Subject $id : <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi

#------------------------------------------------------------------------------#
Title "Running MICA Morphology processing"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"
Info "Not erasing temporal dir: $nocleanup"

# Timer
aloita=$(date +%s)

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Temporary fsa5 directory
Do_cmd ln -s $FREESURFER_HOME/subjects/fsaverage5/ ${dir_surf}

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporary directory
tmp=${tmp}/${RANDOM}_micapipe_post-morpho_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# TRAP in case the script fails
trap cleanup INT TERM

# Make output directory
outDir="$dir_surf"/morphology/
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir=${dir_freesurfer}/surf/

#------------------------------------------------------------------------------#
### Cortical Thickness ###

# Register to fsa5 and apply 10mm smooth
if [[ ! -f ${outDir}/rh.thickness_10mm_fsa5.mgh ]]; then
    for hemi in lh rh; do
        # Convert native file to mgh and save in output directory
        Do_cmd mri_convert ${dataDir}/${hemi}.thickness ${outDir}/${hemi}_thickness.mgh

        Do_cmd mri_surf2surf --hemi ${hemi} \
            --srcsubject ${id} \
            --srcsurfval ${outDir}/${hemi}_thickness.mgh \
            --trgsubject fsaverage5 \
            --trgsurfval ${outDir}/${hemi}.thickness_fsa5.mgh

        Do_cmd mri_surf2surf --hemi ${hemi} \
            --fwhm-trg 10 \
            --srcsubject ${id} \
            --srcsurfval ${outDir}/${hemi}_thickness.mgh \
            --trgsubject fsaverage5 \
            --trgsurfval ${outDir}/${hemi}.thickness_10mm_fsa5.mgh
    done
else
    Info "Subject ${id} cortical thickness is registered to fsa5"
fi

# Register to conte69 and apply 10mm smooth
if [[ ! -f ${outDir}/rh_thickness_10mm_c69-32k.mgh ]]; then
    for hemi in lh rh; do
        [[ $hemi == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=`echo $hemisphere | tr [:lower:] [:upper:]`

        Do_cmd mri_convert ${outDir}/${hemi}_thickness.mgh ${tmp}/${hemi}_thickness.func.gii

        Do_cmd wb_command -metric-resample \
            ${tmp}/${hemi}_thickness.func.gii \
            ${dir_conte69}/${id}_${hemi}_sphereReg.surf.gii \
            ${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii \
            ADAP_BARY_AREA \
            ${tmp}/${hemi}_thickness_c69-32k.func.gii \
            -area-surfs \
            ${dir_surf}/${id}/surf/${hemi}.midthickness.surf.gii \
            ${dir_conte69}/${id}_${hemi}_midthickness_32k_fs_LR.surf.gii

        Do_cmd mri_convert ${tmp}/${hemi}_thickness.func.gii ${outDir}/${hemi}_thickness_c69-32k.mgh

        # Smoothing
        Do_cmd wb_command -metric-smoothing \
            ${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii \
            ${tmp}/${hemi}_thickness_c69-32k.func.gii \
            10 \
            ${tmp}/${hemi}_thickness_10mm_c69-32k.func.gii

        Do_cmd mri_convert ${tmp}/${hemi}_thickness_10mm_c69-32k.func.gii ${outDir}/${hemi}_thickness_10mm_c69-32k.mgh
    done
else
    Info "Subject ${id} cortical thickness is registered to conte69"
fi


#------------------------------------------------------------------------------#
### Curvature ###

# Register to fsa5 and apply 10mm smooth
if [[ ! -f ${outDir}/rh.curv_10mm_fsa5.mgh ]]; then
    for hemi in lh rh; do
        # Convert native file to mgh and save in output directory
        Do_cmd mri_convert ${dataDir}/${hemi}.curv ${outDir}/${hemi}_curv.mgh

        Do_cmd mri_surf2surf --hemi ${hemi} \
            --srcsubject ${id} \
            --srcsurfval ${outDir}/${hemi}_curv.mgh \
            --trgsubject fsaverage5 \
            --trgsurfval ${outDir}/${hemi}.curv_fsa5.mgh

        Do_cmd mri_surf2surf --hemi ${hemi} \
            --fwhm-trg 10 \
            --srcsubject ${id} \
            --srcsurfval ${outDir}/${hemi}_curv.mgh \
            --trgsubject fsaverage5 \
            --trgsurfval ${outDir}/${hemi}.curv_10mm_fsa5.mgh
    done
else
    Info "Subject ${id} curvature is registered to fsa5"
fi

# Register to conte69 and apply 10mm smooth
if [[ ! -f ${outDir}/rh_curv_10mm_c69-32k.mgh ]]; then
    for hemi in lh rh; do
        [[ $hemi == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=`echo $hemisphere | tr [:lower:] [:upper:]`

        Do_cmd mri_convert ${outDir}/${hemi}_curv.mgh ${tmp}/${hemi}_curv.func.gii

        Do_cmd wb_command -metric-resample \
            ${tmp}/${hemi}_thickness.func.gii \
            ${dir_conte69}/${id}_${hemi}_sphereReg.surf.gii \
            ${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii \
            ADAP_BARY_AREA \
            ${tmp}/${hemi}_curv_c69-32k.func.gii \
            -area-surfs \
            ${dir_surf}/${id}/surf/${hemi}.midthickness.surf.gii \
            ${dir_conte69}/${id}_${hemi}_midthickness_32k_fs_LR.surf.gii

        Do_cmd mri_convert ${tmp}/${hemi}_curv.func.gii ${outDir}/${hemi}_curv_c69-32k.mgh

        # Smoothing
        Do_cmd wb_command -metric-smoothing \
            ${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii \
            ${tmp}/${hemi}_curv_c69-32k.func.gii \
            10 \
            ${tmp}/${hemi}_curv_10mm_c69-32k.func.gii

        Do_cmd mri_convert ${tmp}/${hemi}_curv_10mm_c69-32k.func.gii ${outDir}/${hemi}_curv_10mm_c69-32k.mgh
    done
else
    Info "Subject ${id} curvature is registered to conte69"
fi

#------------------------------------------------------------------------------#
# Clean temporary directory and fsaverage5
if [[ $nocleanup == "FALSE" ]]; then Do_cmd rm -rf $tmp ${dir_surf}/fsaverage5; else Info "Mica-pipe tmp directory was not erased: \n\t\t\t${tmp}"; fi

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Post-Morphology processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\tlogs:
$dir_logs/post-morph_*.txt"
echo "${id}, post_morpho, ${status}, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
bids_variables_unset
