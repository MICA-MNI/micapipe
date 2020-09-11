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
#   $4 : Temporal directory (default /tmp)
#
# ONLY for scripting and debugging
# TEST=ON
# FastSurfer
# https://doi.org/10.1016/j.neuroimage.2020.117012
# For future implementation: https://github.com/Deep-MI/FastSurfer

# source utilities
source $MICAPIPE/functions/utilities.sh

BIDS=$1
id=$2
out=$3
tmp=$4

#------------------------------------------------------------------------------#
Title "Running MICA structural processing: Freesurfer"

# Assigns variables names
bids_variables $BIDS $id $out
# print the names on the terminal
bids_print.variables

# Check tmp dir: temporary directory
if [ -z ${tmp} ]; then tmp=/tmp; fi
tmp=${tmp}/${RANDOM}_proc-freesurfer_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# BIDS T1w processing
N=${#bids_T1ws[@]} # total number of T1w

# End script if no T1 are found
if [ "$N" -lt 1 ]; then Error "Subject $id doesn't have T1 on: \n\t\t\t${subject_bids}/anat"; exit; fi

# Define SUBJECTS_DIR for freesurfer processing as a global variable
# Will work on a temporal directory
export SUBJECTS_DIR=${tmp}
if [ ! -d ${tmp}/nii ]; then mkdir ${tmp}/nii; fi

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
Do_cmd rm -rf $tmp

Info "Check log file:\n\t\t\t ${dir_logs}/recon-all.log"

# Notification of completition
if [ -f ${dir_freesurfer}/mri/T1.mgz ]; then status="COMPLETED"; else status="ERROR"; fi

Title "Freesurfer recon-all processing ended: ${status}\n\t\t\tlogs:${dir_logs}/proc_freesurfer.txt"
echo "${id}, FREESURFER, ${status}, `whoami`, $(date), `printf "%0.3f\n" ${eri}`" >> ${out}/brain-proc.csv
