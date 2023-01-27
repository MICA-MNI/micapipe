#!/bin/bash
#
# Compute geodesic distance along cortical mesh:
#
# Generates geodesic distance matrices along native cortical surface using specified parcellations
#
# This workflow makes use of wb_command and custom python scripts
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
# nocleanup=$5
threads=$6
# tmpDir=$7
PROC=$8
export OMP_NUM_THREADS=$threads
export here=$(pwd)

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
module_json="${dir_QC}/${idBIDS}_module-GD.json"
micapipe_check_json_status "${module_json}" "GD"

# Define output directory
outPath="${proc_struct}/surf/geo_dist"

# wb_command
workbench_path=$(which wb_command)

# Check PARCELLATIONS
parcellations=($(find "${dir_volum}" -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"))
if [ "${#parcellations[*]}" -eq "0" ]; then Error "Subject $id doesn't have -post_structural processing"; exit; fi

#------------------------------------------------------------------------------#
Title "Geodesic distance analysis\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# creates output directory if it doesn't exist
[[ ! -d "$outPath" ]] && mkdir -p "$outPath"

#------------------------------------------------------------------------------#
# Compute geodesic distance on all parcellations
for seg in "${parcellations[@]}"; do ((N++))
    parc=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    lh_annot="${dir_subjsurf}/label/lh.${parc}_mics.annot"
    rh_annot="${dir_subjsurf}/label/rh.${parc}_mics.annot"
    outName="${outPath}/${idBIDS}_space-fsnative_atlas-${parc}_GD"
    if [ -f "${outName}.txt" ]; then
        Info "Geodesic Distance on $parc, already exists"; ((Nsteps++))
    else
        Info "Computing Geodesic Distance on $parc"
        Do_cmd python "$MICAPIPE"/functions/geoDistMapper.py "$lh_midsurf" "$rh_midsurf" "$outName" "$lh_annot" "$rh_annot" "$workbench_path"
        if [[ -f "${outName}.txt" ]]; then ((Nsteps++)); fi
    fi
done

Do_cmd rm -rf "${outPath}"/*.func.gii

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status GD
micapipe_procStatus "${id}" "${SES/ses-/}" "GD" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "GD" "${module_json}"
bids_variables_unset
