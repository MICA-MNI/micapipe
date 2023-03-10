#!/bin/bash
#!usr/bin/python
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
nocleanup=$5
threads=$6
tmpDir=$7
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

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-dist.json"
micapipe_check_json_status "${module_json}" "disty"

# Define output directory
outPath="${subject_dir}/dist"

# wb_command
workbench_path=$(which wb_command)

# Check PARCELLATIONS
parcellations=($(find "${dir_volum}" -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"))
if [ "${#parcellations[*]}" -eq "0" ]; then Error "Subject $id doesn't have -post_structural processing"; exit; fi

#------------------------------------------------------------------------------#
Title "Matrix distance analysis\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmpDir=${tmpDir}/${RANDOM}_micapipe_dist_${idBIDS}
Do_cmd mkdir -p "$tmpDir"

# creates output directory if it doesn't exist
[[ ! -d "$outPath" ]] && Do_cmd mkdir -p "$outPath"

# VARIABLES
T1_seg_cerebellum="${dir_volum}/${idBIDS}_space-nativepro_T1w_atlas-cerebellum.nii.gz"
T1_seg_subcortex="${dir_volum}/${idBIDS}_space-nativepro_T1w_atlas-subcortical.nii.gz"

# Set temporary files
T1w_subc=${tmpDir}/T1w_seg_subcortical.nii.gz
T1w_cere=${tmpDir}/T1w_seg_cderebellum.nii.gz

# Threshold unwanted ROIs
Do_cmd fslmaths "$T1_seg_cerebellum" -bin -mul 100 -add "$T1_seg_cerebellum" "$T1w_cere"
Do_cmd fslmaths "$T1_seg_subcortex" -thr 16 -uthr 16 -binv -mul "$T1_seg_subcortex" "$T1w_subc"

#------------------------------------------------------------------------------#
# Compute Geodesic distance matrix on fsLR-5k

for HEMI in L R; do
  outFile="${outPath}/${idBIDS}_hemi-${HEMI}_surf-fsLR-5k_desc-GD.txt"
  if [[ ! -f "$outFile" ]]; then ((N++))
    fsLR5k="${dir_conte69}/${idBIDS}_hemi-${HEMI}_space-nativepro_surf-fsLR-5k_label-pial.surf.gii"
    Info "Calculating GD on fsLR-5k hemi-${HEMI}"
    Do_cmd python "${MICAPIPE}"/functions/dist.py "$fsLR5k" "${outPath}/${idBIDS}_hemi-${HEMI}_surf-fsLR-5k_desc-GD.txt" "GD" "fsLR-5k"
    if [[ -f "$outFile" ]]; then ((Nsteps++)); fi
  else
    Info "GD connectome of fsLF-5k_hemi-${HEMI} exist"; ((Nsteps++)); ((N++))
  fi
done

#------------------------------------------------------------------------------#
# Compute Euclidean distance matrix on all parcellations
for seg in "${parcellations[@]}"; do ((N++))
    parc=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    lh_annot="${dir_subjsurf}/label/lh.${parc}_mics.annot"
    rh_annot="${dir_subjsurf}/label/rh.${parc}_mics.annot"
    outName="${outPath}/${idBIDS}_atlas-${parc}_EDM"
    if [ -f "${outName}.txt" ]; then
        Info "Distance on $parc, already exists"; ((Nsteps++))
    else
      # Remove the medial walldwi_all
      seg_ctx=${tmpDir}/${parc}.nii.gz
      seg_ctxsub=${tmpDir}/${parc}+sub.nii.gz
      seg_full=${tmpDir}/${parc}+sub+cereb.nii.gz
      cp "$parcellations" "${seg_ctx}"
      # remove the medial wall
      for i in 1000 2000; do Do_cmd fslmaths "$seg" -thr "$i" -uthr "$i" -binv -mul "$seg" "$seg_ctx"; done
        Do_cmd fslmaths "$seg_ctx" -binv -mul "$T1w_subc" -add "$seg_ctx" "$seg_ctxsub" -odt int # added the subcortical parcellation
        Do_cmd fslmaths "$seg_ctxsub" -binv -mul "$T1w_cere" -add "$seg_ctxsub" "${seg_full}" -odt int # added the cerebellar parcellation
        Info "Distance on $parc"
        # get the center of mass
        3dCM -all_rois "${seg_full}" > ${tmpDir}/${atlas}.txt
        # sort the output table
        sed -e 1d ${tmpDir}/${atlas}.txt | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/ /g' | sed -e 's/#ROI /\n/g' | sed -e 1d > ${tmpDir}/tmp.txt && mv ${tmpDir}/tmp.txt ${tmpDir}/${atlas}.txt
        # more sorting
        sed -e 's/  / /g' ${tmpDir}/${atlas}.txt | sed -e 's/ /,/g' | sed -e 's/.$//' > ${tmpDir}/tmp.txt && mv ${tmpDir}/tmp.txt ${tmpDir}/${atlas}.txt
        # Calculate the distance
        Do_cmd python "${MICAPIPE}"/functions/dist.py ${tmpDir}/${atlas}.txt "$outName" "EDM" "${parc}"
        if [[ -f "${outName}.txt" ]]; then ((Nsteps++)); fi
    fi
done

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status dist
micapipe_procStatus "${id}" "${SES/ses-/}" "dist" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "dist" "${module_json}"
bids_variables_unset
