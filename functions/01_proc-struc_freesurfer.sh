#!/bin/bash
#
# T1w Structural processing with bash:
#
# Preprocessing workflow for structural T1w.
#
# This workflow makes use of FREESURFER, FSL (fslchfiletype)
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
# ONLY for scripting and debugging
# TEST=ON
# FastSurfer
# https://doi.org/10.1016/j.neuroimage.2020.117012
# For future implementation: https://github.com/Deep-MI/FastSurfer

BIDS=$1
id=$2
out=$3
SES=$4
PROC=$5
nocleanup=$6
here=`pwd`

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

# Stop if freesurfer has finished without errors
if grep -q "finished without error" ${dir_freesurfer}/scripts/recon-all.log; then
Warning "Subject ${id} has Freesurfer
                    > If you want to re-run for QC purposes try it manually
                    > If you want to run again this step first erase all the outputs with:
                      mica_cleanup -sub <subject_id> -out <derivatives> -bids <BIDS_dir> -proc_fresurfer";
Do_cmd exit
fi

#------------------------------------------------------------------------------#
Title "Running MICA structural processing: Freesurfer"
micapipe_software
# print the names on the terminal
bids_print.variables
Info "Not erasing temporal dir: $nocleanup"

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporal directory
tmp=${tmp}/${RANDOM}_micapipe_proc-freesurfer_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# BIDS T1w processing
N=${#bids_T1ws[@]} # total number of T1w

# End script if no T1 are found
if [ "$N" -lt 1 ]; then Error "Subject $id doesn't have T1 on: \n\t\t\t${subject_bids}/anat"; exit; fi

# Define SUBJECTS_DIR for freesurfer processing as a global variable
# Will work on a temporal directory
export SUBJECTS_DIR=${tmp}
if [ ! -d ${tmp}/nii ]; then Do_cmd mkdir ${tmp}/nii; fi

# Copy all the T1 from the BIDS directory to the TMP
# transform to NIFTI (if == NIFTI_GZ)
for t1 in ${bids_T1ws[@]}; do
    t1_name=`echo $t1 | awk -F 'anat/' '{print $2}'`
    if [[ $t1 == *'.gz'* ]]; then
          Do_cmd fslchfiletype NIFTI $t1 ${tmp}/nii/${t1_name/.gz/}
    else  Do_cmd cp $t1 ${tmp}/nii/
    fi
done

# List of Files for processing
fs_cmd=$(echo "-i $(echo ${tmp}/nii/*nii | sed 's: : -i :g')")

# Perform recon-all surface registration
Do_cmd recon-all -cm -all "$fs_cmd" -s "$id"

# Copy the recon-all log to our MICA-log Directory
Do_cmd cp -v ${tmp}/${id}/scripts/recon-all.log ${dir_logs}/recon-all.log

# Copy results to  freesurfer's SUBJECTS_DIR directory
Do_cmd cp -rv ${tmp}/${id} $dir_surf

# Remove temporal directory
if [[ $nocleanup == "FALSE" ]]; then Do_cmd rm -rf $tmp; else Info "Mica-pipe tmp directory was not erased: \n\t\t\t${tmp}"; fi

Info "Check log file:\n\t\t\t ${dir_logs}/recon-all.log"

# Notification of completition
if [ -f ${dir_freesurfer}/mri/T1.mgz ]; then status="COMPLETED"; else status="ERROR"; fi

Title "Freesurfer recon-all processing ended: ${status}\n\tlogs:
`ls ${dir_logs}/proc-freesurfer*.txt`"

echo "${id}, FREESURFER, ${status}, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
