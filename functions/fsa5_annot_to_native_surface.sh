#!/bin/bash
# Maps fsa5 template parcellation to native surface

source ~/.sge_profile
source ${MICAPIPE}/functions/utilities.sh

# Set parameters
subject=$1 # HC12
fsDir=$2 # /host/fladgate/local_raid/MICA-MTL/"$sub"/scan_session_01/proc_struct/surfaces/
annotPath=$3 # Full path to annot, normally in micapipe: /data_/mica1/03_projects/jessica/micapipe/parcellations/
annotName=$4 # e.g. schaefer_100

#####################
# Step 1: Map fsa5 annotation file to subject's native surface.
#####################

lhAnnot="$annotPath"/lh."$annotName"_fsa5.annot
rhAnnot="$annotPath"/rh."$annotName"_fsa5.annot

export FREESURFER_HOME=/data_/mica1/01_programs/Freesurfer-6.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR="$fsDir"

labelDir="$fsDir"/"$subject"/label/

# map annotation to subject space
mri_surf2surf --srcsubject fsaverage5 --trgsubject "$subject" --hemi lh \
    --sval-annot "$lhAnnot" \
    --tval       "$labelDir"/lh."$annotName"_mics.annot
mri_surf2surf --srcsubject fsaverage5 --trgsubject $subject --hemi rh \
    --sval-annot "$rhAnnot" \
    --tval       "$labelDir"/rh."$annotName"_mics.annot

#####################
# Step 2: Read subject annotation file and export to csv in same directory.
#####################

python ${MICAPIPE}/functions/annot2csv.py "$subject" "$labelDir" "$annotName"

