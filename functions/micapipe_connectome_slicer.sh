#!/bin/bash
#
# This script slices the connectomes using LUT
#
#
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
# ONLY for scripting and debugging
# TEST=ON

BIDS=$1
id=$2
out=$3
SES=$4
PROC=$5
nocleanup=$6
here=`pwd`

# Get the real path of the Inputs
out=`realpath $out`
BIDS=`realpath $BIDS`
id=${id/sub-/}

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

# Check inputs: DWI post TRACTOGRAPHY
lut_sc="${util_lut}/lut_subcortical-cerebellum_mics.csv"
tracts=10M # <<<<<<<<<<<<<<<<<< Number of streamlines

# Check inputs
if [ ! -f $fod ]; then Error "Subject $id doesn't have FOD:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $dwi_b0 ]; then Error "Subject $id doesn't have dwi_b0:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $mat_dwi_affine ]; then Error "Subject $id doesn't have an affine mat from T1nativepro to DWI space:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $dwi_5tt ]; then Error "Subject $id doesn't have 5tt in dwi space:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $T1_seg_cerebellum ]; then Error "Subject $id doesn't have cerebellar segmentation:\n\t\tRUN -post_structural"; exit; fi
if [ ! -f $T1_seg_subcortex ]; then Error "Subject $id doesn't have subcortical segmentation:\n\t\tRUN -post_structural"; exit; fi

#------------------------------------------------------------------------------#
Title "Running MICA connectome slicer"
micapipe_software

#	Timer
aloita=$(date +%s)
here=`pwd`
Nparc=0

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporal directory
tmp=${tmp}/${RANDOM}_micapipe_post-dwi_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

cd ${tmp}

# -----------------------------------------------------------------------------------------------
# Prepare the segmentatons
parcellations=`find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"`

# -----------------------------------------------------------------------------------------------
# Build the Connectomes
for seg in $parcellations; do
    parc_name=`echo ${seg/.nii.gz/} | awk -F 'nativepro_' '{print $2}'`
    nom=${dwi_cnntm}/${id}_${tracts}_${parc_name}
    lut="${util_lut}/lut_${parc_name}_mics.csv"
    dwi_cortex=$tmp/${id}_${parc_name}-cor_dwi.nii.gz # Segmentation in dwi space

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical connectomes
    Info "Slicing $parc_name cortical connectome"
    Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn="${nom}_cor-connectome.txt" --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
    # Calculate the edge lenghts
    Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn="${nom}_cor-edgeLengths.txt" --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical connectomes (-sub)
    Info "Slicing $parc_name cortical-subcortical connectome"
    dwi_cortexSub=$tmp/${id}_${parc_name}-sub_dwi.nii.gz
    Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn="${nom}_sub-connectome.txt" --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
    # Calculate the edge lenghts
    Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn="${nom}_sub-edgeLengths.txt" --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical-Cerebellar connectomes (-sub-cereb)
    Info "Slicing $parc_name cortical-subcortical-cerebellum connectome"
    Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn="${nom}_full-connectome.txt" --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
    Rscript ${MICAPIPE}/functions/connectome_slicer.R --conn="${nom}_full-edgeLengths.txt" --lut1=${lut_sc} --lut2=${lut} --mica=${MICAPIPE}
done

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
if [[ -z $nocleanup ]]; then Do_cmd rm -rf $tmp; else Info "tmp directory was not erased: ${tmp}"; fi
cd $here

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
if [ "$Nparc" -eq 54 ]; then status="DONE"; else status="ERROR missing a connectome: "; fi
Title "CONNECTOME SLICER processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:
\t\tNumber of connectomes: `printf "%02d" $Nparc`/56
\tlogs:
`ls ${dir_logs}/post-dwi_*.txt`"
