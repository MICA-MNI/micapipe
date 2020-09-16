#!/bin/bash
#
# DWI structural processing with bash:
#
# Preprocessing workflow for diffusion MRI.
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
#TEST=ON

BIDS=$1
id=$2
out=$3
PROC=$4
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

# Check inputs: DWI
if [ "${#bids_dwis[@]}" -lt 1 ]; then Error "Subject $id doesn't have DWIs:\n\t\t TRY <ls -l ${subject_bids}/dwi/>"; exit; fi

#------------------------------------------------------------------------------#
Title "Running MICA Diffusion Weighted Imaging processing"
# print the names on the terminal
bids_print.variables-dwi

# GLOBAL variables for this script
Info "ANTs will use $CORES CORES"

#	Timer
aloita=$(date +%s)
here=`pwd`

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporal directory
tmp=${tmp}/${RANDOM}_micapipe_proc-dwi_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Set basic parameters.
dwi2fodAlgorithm=msmt_csd  # msmt_csd or csd

#------------------------------------------------------------------------------#
# DWI processing
# Image denoising must be performed as the first step of the image-processing pipeline.
# Interpolation or smoothing in other processing steps, such as motion and distortion correction,
# may alter the noise characteristics and thus violate the assumptions upon which MP-PCA is based.
dwi_n4=${proc_dwi}/${id}_dwi_dnsN4.mif
dwi_cat=${tmp}/dwi_concatenate.mif
dwi_dns=${tmp}/dwi_concatenate_denoised.mif

if [[ ! -f ${dwi_n4} ]] ; then
# Concatenate shells -if only one shell then just convert to mif and rename.
      for dwi in ${bids_dwis[@]}; do
            dwi_nom=`echo $dwi | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}'`
            bids_dwi_str=`echo $dwi | awk -F . '{print $1}'`
            Do_cmd mrconvert $dwi -json_import ${bids_dwi_str}.json -fslgrad ${bids_dwi_str}.bvec ${bids_dwi_str}.bval ${tmp}/${dwi_nom}.mif
      done

      # Concatenate shells and convert to mif.
      Do_cmd mrcat ${tmp}/*.mif $dwi_cat -nthreads $CORES

      # Denoise DWI and calculate residuals
      Do_cmd dwidenoise $dwi_cat $dwi_dns -nthreads $CORES
      Do_cmd mrcalc $dwi_cat $dwi_dns -subtract ${proc_dwi}/${id}_dwi_residuals.mif -nthreads $CORES

      # Bias field correction DWI
      Do_cmd dwibiascorrect ants $dwi_dns $dwi_n4 -force -nthreads $CORES
else
      Info "Subject ${id} has DWI in mif, denoised and concatenaded"
fi

#------------------------------------------------------------------------------#
# dwifslpreproc and TOPUP preparations
dwi_corr=${proc_dwi}/${id}_dwi_corr.mif

if [[ ! -f $dwi_corr ]]; then
      # Get parameters
      ReadoutTime=`mrinfo $dwi_n4 -property TotalReadoutTime`
      pe_dir=`mrinfo $dwi_n4 -property PhaseEncodingDirection`

      # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
      dwi_4proc=${tmp}/dwi_n4_even.mif
      dim=`mrinfo $dwi_n4 -size`
      dimNew=($(echo $dim | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
      mrconvert $dwi_n4 $dwi_4proc -coord 0 0:${dimNew[0]} -coord 1 0:${dimNew[1]} -coord 2 0:${dimNew[2]} -coord 3 0:end -force

      # Get the mean b-zero (un-corrected)
      dwiextract -nthreads $CORES $dwi_n4 - -bzero | mrmath - mean ${tmp}/b0_meanMainPhase.mif -axis 3

      # Processing the reverse encoding b0
      if [ -f $dwi_reverse ]; then
            b0_pair_tmp=${tmp}/b0_pair_tmp.mif
            b0_pair=${tmp}/b0_pair.mif
            dwi_reverse_str=`echo $dwi_reverse | awk -F . '{print $1}'`

            Do_cmd mrconvert $dwi_reverse -json_import ${dwi_reverse_str}.json -fslgrad ${dwi_reverse_str}.bvec ${dwi_reverse_str}.bval ${tmp}/b0_ReversePhase.mif
            dwiextract ${tmp}/b0_ReversePhase.mif - -bzero | mrmath - mean ${tmp}/b0_meanReversePhase.mif -axis 3 -nthreads $CORES
            Do_cmd mrcat ${tmp}/b0_meanMainPhase.mif ${tmp}/b0_meanReversePhase.mif $b0_pair_tmp -nthreads $CORES

            # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
            dim=`mrinfo $b0_pair_tmp -size`
            dimNew=($(echo $dim | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
            mrconvert $b0_pair_tmp $b0_pair -coord 0 0:${dimNew[0]} -coord 1 0:${dimNew[1]} -coord 2 0:${dimNew[2]} -coord 3 0:end -force

            opt="-rpe_pair -se_epi ${b0_pair}"
      else
            opt='-rpe_none'
      fi

      Info "dwifslpreproc parameters:"
      Note "Shell values        :" "`mrinfo $dwi_4proc -shell_bvalues`"
      Note "DWI main dimensions :" "`mrinfo $dwi_4proc -size`"
      Note "DWI rpe dimensions  :" "`mrinfo $b0_pair -size`"
      Note "pe_dir              :" $pe_dir
      Note "Readout Time        :" $ReadoutTime

      # Preprocess each shell
      # DWIs all acquired with a single fixed phase encoding; but additionally a
      # pair of b=0 images with reversed phase encoding to estimate the inhomogeneity field:
      Do_cmd dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -align_seepi -eddy_options " --data_is_shelled" -nthreads $CORES -nocleanup -scratch $tmp

else
      Info "Subject ${id} has a DWI processed with dwifslpreproc"
fi

#------------------------------------------------------------------------------#
# Registration of corrected DWI-b0 to T1nativepro
dwi_mask=${proc_dwi}/${id}_dwi_mask.mif
dwi_b0=${proc_dwi}/${id}_dwi_b0.nii.gz # This should be a NIFTI for compatibility with ANTS
dwi_in_T1nativepro=${proc_struct}/${id}_dwi_nativepro.nii.gz
T1nativepro_in_dwi=${proc_dwi}/${id}_t1w_dwi.nii.gz
str_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_
mat_dwi_affine=${str_dwi_affine}0GenericAffine.mat

if [[ ! -f $T1nativepro_in_dwi ]]; then
      # Create a binary mask of the DWI
      Do_cmd dwi2mask -force -nthreads $CORES $dwi_corr $dwi_mask

      # Corrected DWI-b0s mean for registration
      dwiextract -force -nthreads $CORES $dwi_corr - -bzero | mrmath - mean $dwi_b0 -axis 3

      # Register DWI-b0 mean corrected to T1nativepro
      Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro -m $dwi_b0 -o $str_dwi_affine -t a -n $CORES -p d
      # Apply transformation DWI-b0 space to T1nativepro
      Do_cmd antsApplyTransforms -d 3 -i $dwi_b0 -r $T1nativepro -t $mat_dwi_affine -o $dwi_in_T1nativepro -v -u int
      # Apply inverse transformation T1nativepro to DWI-b0 space
      Do_cmd antsApplyTransforms -d 3 -i $T1nativepro -r $dwi_b0 -t [$mat_dwi_affine,1] -o $T1nativepro_in_dwi -v -u int
else
      Info "Subject ${id} has a T1nativepro in DWI-b0 space"
fi

#------------------------------------------------------------------------------#
# Get some basic metrics.
if [[ ! -f $proc_dwi/${id}_FA.mif ]]; then
      Do_cmd dwi2tensor -nthreads $CORES $dwi_corr $tmp/${id}_tensor.mif
      Do_cmd tensor2metric -nthreads $CORES -fa $proc_dwi/${id}_FA.mif -adc $proc_dwi/${id}_ADC.mif $tmp/${id}_tensor.mif
else
      Info "Subject ${id} has DWI tensor metrics"
fi


#------------------------------------------------------------------------------#
# Prepare tractography
if [[ ! -f $proc_dwi/${id}_FOD_WM_norm.mif ]]; then
      if [[ $dwi2fodAlgorithm == msmt_csd ]]; then
            Do_cmd dwi2response dhollander -nthreads $CORES $dwi_corr \
                ${proc_dwi}/${id}_dhollander_response_wm.txt \
                ${proc_dwi}/${id}_dhollander_response_gm.txt \
                ${proc_dwi}/${id}_dhollander_response_csf.txt \
                -mask $dwi_mask
            Do_cmd dwi2fod -nthreads $CORES $dwi2fodAlgorithm $dwi_corr \
                ${proc_dwi}/${id}_dhollander_response_wm.txt \
                ${tmp}/wmFOD.mif \
                -mask $dwi_mask
            Do_cmd mtnormalise -nthreads $CORES -mask $dwi_mask \
                ${tmp}/wmFOD.mif $proc_dwi/${id}_FOD_WM_norm.mif
      else
            Do_cmd dwi2response tournier -nthreads $CORES \
                $dwi_corr \
                ${proc_dwi}/${id}_tournier_response.txt \
                -mask $dwi_mask
            Do_cmd dwi2fod -nthreads $CORES $dwi2fodAlgorithm \
                $dwi_corr \
                ${proc_dwi}/${id}_tournier_response.txt  \
                ${tmp}/FOD.mif \
                -mask $dwi_mask
            Do_cmd mtnormalise -nthreads $CORES -mask $dwi_mask \
                ${tmp}/FOD.mif $proc_dwi/${id}_FOD_norm.mif
      fi
else
      Info "Subject ${id} has DWI metrics prepared for tractography"
fi


# -----------------------------------------------------------------------------------------------
# Clean temporal directory
Do_cmd rm -rf $tmp

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
status="DWI-proc"
Title "DWI processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\tlogs:
`ls ${dir_logs}/proc-dwi_*.txt`"

echo "${id}, proc_dwi, $fini N=`printf "%02d" $Nfiles`/21, `whoami`, ``,$(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
