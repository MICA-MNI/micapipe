#!/bin/bash
#
# DWI POST structural TRACTOGRAPHY processing with bash:
#
# POST processing workflow for diffusion MRI TRACTOGRAPHY.
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
PROC=$4
nocleanup=$5
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
bids_variables $BIDS $id $out

# Check inputs: DWI post TRACTOGRAPHY
fod=$proc_dwi/${id}_wm_fod_norm.mif
dwi_5tt=${proc_dwi}/${id}_dwi_5tt.nii.gz
T1str_nat=${id}_t1w_${res}mm_nativepro
T1_seg_cerebellum=${dir_volum}/${T1str_nat}_cerebellum.nii.gz
T1_seg_subcortex=${dir_volum}/${T1str_nat}_subcortical.nii.gz
dwi_b0=${proc_dwi}/${id}_dwi_b0.nii.gz
mat_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_0GenericAffine.mat

# Check inputs
if [ ! -f $fod ]; then Error "Subject $id doesn't have FOD:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $dwi_b0 ]; then Error "Subject $id doesn't have dwi_b0:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $mat_dwi_affine ]; then Error "Subject $id doesn't have an affine mat from T1nativepro to DWI space:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $dwi_5tt ]; then Error "Subject $id doesn't have 5tt in dwi space:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f $T1_seg_cerebellum ]; then Error "Subject $id doesn't have cerebellar segmentation:\n\t\tRUN -post_structural"; exit; fi
if [ ! -f $T1_seg_subcortex ]; then Error "Subject $id doesn't have subcortical segmentation:\n\t\tRUN -post_structural"; exit; fi

#------------------------------------------------------------------------------#
Title "Running MICA POST-DWI processing (Tractography)"

#	Timer
aloita=$(date +%s)
here=`pwd`
Nparc=0
tracts=40M # <<<<<<<<<<<<<<<<<< Number of stremalines

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporal directory
tmp=${tmp}/${RANDOM}_micapipe_post-dwi_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Create Connectomes directory for the outpust
[[ ! -d ${dwi_cnntm} ]] && Do_cmd mkdir -p ${dwi_cnntm}
cd ${tmp}

# -----------------------------------------------------------------------------------------------
# Generate and optimize probabilistic tracts
Info "Building the 100 million streamlines connectome!!!"
tck=${tmp}/DWI_tractogram_${tracts}.tck
weights=${tmp}/SIFT2_${tracts}.txt
Do_cmd tckgen -nthreads $CORES \
    $fod \
    $tck \
    -act $dwi_5tt \
    -crop_at_gmwmi \
    -seed_dynamic $fod \
    -maxlength 300 \
    -minlength 10 \
    -angle 22.5 \
    -power 1.0 \
    -backtrack \
    -select ${tracts} \
    -step .5 \
    -cutoff 0.06 \
    -algorithm iFOD2

# SIFT2
Do_cmd tcksift2 -nthreads $CORES $tck $fod $weights
# TDI for QC
Do_cmd tckmap -vox 1,1,1 -dec -nthreads $CORES $tck $proc_dwi/${id}_tdi_iFOD2-${tracts}.mif

# -----------------------------------------------------------------------------------------------
# Prepare the segmentatons
parcellations=`find ${dir_volum} -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"`
T1_seg_cerebellum=${dir_volum}/${T1str_nat}_cerebellum.nii.gz
T1_seg_subcortex=${dir_volum}/${T1str_nat}_subcortical.nii.gz
dwi_cere=${proc_dwi}/${id}_dwi_cerebellum.nii.gz
dwi_subc=${proc_dwi}/${id}_dwi_subcortical.nii.gz

if [[ ! -f $dwi_cere ]]; then Info "Registering Cerebellar parcellation to DWI-b0 space"
      Do_cmd antsApplyTransforms -d 3 -e 3 -i $T1_seg_cerebellum -r $dwi_b0 -n linear -t [$mat_dwi_affine,1] -o $dwi_cere -v -u int
      if [[ -f $dwi_cere ]]; then ((Nparc++)); fi
else Info "Subject ${id} has a Cerebellar segmentation in DWI space"; ((Nsteps++)); fi

if [[ ! -f $dwi_subc ]]; then Info "Registering Subcortical parcellation to DWI-b0 space"
      Do_cmd antsApplyTransforms -d 3 -e 3 -i $T1_seg_subcortex -r $dwi_b0 -n linear -t [$mat_dwi_affine,1] -o $dwi_subc -v -u int
      if [[ -f $dwi_subc ]]; then ((Nparc++)); fi
else Info "Subject ${id} has a Subcortical segmentation in DWI space"; ((Nsteps++)); fi

# -----------------------------------------------------------------------------------------------
# Build the Connectomes
for seg in $parcellations; do
    parc_name=`echo ${seg/.nii.gz/} | awk -F 'nativepro_' '{print $2}'`
    nom=${dwi_cnntm}/${id}_${tracts}_${parc_name}
    dwi_cortex=$tmp/${id}_${parc_name}_dwi.nii.gz # Segmentation in dwi space
    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical connectomes
    Info "Building $parc_name cortical connectome"
    # Take parcellation into DWI space
    Do_cmd antsApplyTransforms -d 3 -e 3 -i $seg -r $dwi_b0 -n linear -t [$mat_dwi_affine,1] -o $dwi_cortex -v -u int
    # Build the Cortical connectomes
    Do_cmd tck2connectome -nthreads $CORES \
        $tck $dwi_cortex ${nom}.txt \
        -tck_weights_in $weights -out_assignments ${nom}_assignments.txt \

    Do_cmd tck2connectome -nthreads $CORES \
        $tck $dwi_cortex ${nom}_edgeLengths.txt \
        -tck_weights_in $weights -scale_length -stat_edge mean
    if [[ -f ${nom}.txt ]]; then ((Nparc++)); fi

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical connectomes (-sub)
    Info "Building $parc_name cortical-subcortical connectome"
    dwi_cortexSub=$tmp/${id}_${parc_name}-sub_dwi.nii.gz
    Do_cmd fslmaths $dwi_cortex -add $dwi_subc $dwi_cortexSub -odt int # added the subcortical parcellation

    Do_cmd tck2connectome -nthreads $CORES \
        $tck $dwi_cortexSub ${nom}-sub.txt \
        -tck_weights_in $weights -out_assignments ${nom}-sub_assignments.txt \

    Do_cmd tck2connectome -nthreads $CORES \
        $tck $dwi_cortexSub ${nom}-sub_edgeLengths.txt \
        -tck_weights_in $weights -scale_length -stat_edge mean
    if [[ -f ${nom}-sub.txt ]]; then ((Nparc++)); fi

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical-Cerebellar connectomes (-sub-cereb)
    Info "Building $parc_name cortical-subcortical-cerebellum connectome"
    dwi_all=$tmp/${id}_${parc_name}-all_dwi.nii.gz
    Do_cmd fslmaths $dwi_cere -add 100 -add $dwi_cortexSub $dwi_all -odt int # added the subcortical parcellation

    Do_cmd tck2connectome -nthreads $CORES \
        $tck $dwi_all ${nom}-all.txt \
        -tck_weights_in $weights -out_assignments ${nom}-all_assignments.txt \

    Do_cmd tck2connectome -nthreads $CORES \
        $tck $dwi_all ${nom}-all_edgeLengths.txt \
        -tck_weights_in $weights -scale_length -stat_edge mean
    if [[ -f ${nom}-all.txt ]]; then ((Nparc++)); fi

done

# -----------------------------------------------------------------------------------------------
# Compute Auto-Tractography (Future Release version)
# >>>>>> https://github.com/lconcha/auto_tracto

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
if [[ -z $nocleanup ]]; then Do_cmd rm -rf $tmp; fi
cd $here

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
if [ "$Nparc" -eq 36 ]; then status="DONE"; else status="ERROR missing a connectome: "; fi
Title "TEST-DWI-post TRACTOGRAPHY processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:
\t\tNumber of connectomes: `printf "%02d" $Nparc`/60
\tlogs:
`ls ${dir_logs}/post-dwi_*.txt`"
# Print QC stamp
echo "${id}, post_dwi, $status N=`printf "%02d" $Nparc`/36, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
