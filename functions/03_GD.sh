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
    export MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Define output directory
outPath="${proc_struct}/surfaces/geo_dist"

# wb_command
workbench_path=$(which wb_command)

# Check inputs: freesurfer space T1 (to make sure freesurfer was run)
if [ ! -f "${T1freesurfr}" ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi
if [ ! -f "${lh_midsurf}" ]; then Error "Subject $id doesn't have left hemisphere midsurface gifti file"; exit; fi
if [ ! -f "${rh_midsurf}" ]; then Error "Subject $id doesn't have right hemisphere midsurface gifti file"; exit; fi

# Check PARCELLATIONS
parcellations=($(find "${dir_volum}" -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"))
if [ "${#parcellations[*]}" -eq "0" ]; then Error "Subject $id doesn't have -post_structural processing"; exit; fi

#------------------------------------------------------------------------------#
Title "Geodesic distance analysis\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"

# Check IF output exits and WARNING
N=$(ls "${outPath}"/* 2>/dev/null | wc -l)
if [ "$N" -gt 2 ]; then Warning "Existing connectomes will be skipped!! If you want to re-run -GD first clean the outpus:
          micapipe_cleanup -GD -sub ${id} -out ${out} -bids ${BIDS}"; fi

#	Timer
aloita=$(date +%s)
Nsteps=0
# creates output directory if it doesn't exist
[[ ! -d "$outPath" ]] && mkdir -p "$outPath"

#------------------------------------------------------------------------------#
# Compute geodesic distance on all parcellations
N=${#parcellations[*]}
for seg in "${parcellations[@]}"; do
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
if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "Post-GD processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\tSteps completed : $(printf "%02d" "$Nsteps")/$(printf "%02d" "$N")
\tStatus          : ${status}
\tCheck logs      : $(ls "${dir_logs}"/GD_*.txt)"
micapipe_procStatus "${id}" "${SES/ses-/}" "GD" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "GD" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
bids_variables_unset
