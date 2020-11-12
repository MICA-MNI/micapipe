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

# Check inputs: DWI
if [ "${#bids_dwis[@]}" -lt 1 ]; then Error "Subject $id doesn't have DWIs:\n\t\t TRY <ls -l ${subject_bids}/dwi/>"; exit; fi
if [ ! -f ${T1_MNI152_InvWarp} ]; then Error "Subject $id doesn't have T1_nativepro warp to MNI152.\n\t\tRun -proc_structural"; exit; fi
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro.\n\t\tRun -proc_structural"; exit; fi
if [ ! -f ${T15ttgen} ]; then Error "Subject $id doesn't have a 5tt volume in nativepro space.\n\t\tRun -proc_structural"; exit; fi

#------------------------------------------------------------------------------#
Title "Running MICA Diffusion Weighted Imaging processing"
micapipe_software
# print the names on the terminal
bids_print.variables-dwi

# GLOBAL variables for this script
Info "ANTs will use $CORES CORES"

#	Timer
aloita=$(date +%s)
Nsteps=0

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporal directory
tmp=${tmp}/${RANDOM}_micapipe_proc-dwi_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

cd $tmp
#------------------------------------------------------------------------------#
# DWI processing
# Image denoising must be performed as the first step of the image-processing pipeline.
# Interpolation or smoothing in other processing steps, such as motion and distortion correction,
# may alter the noise characteristics and thus violate the assumptions upon which MP-PCA is based.
dwi_n4=${proc_dwi}/${id}_dwi_dnsN4.mif
dwi_cat=${tmp}/dwi_concatenate.mif
dwi_dns=${tmp}/dwi_concatenate_denoised.mif
dwi_corr=${proc_dwi}/${id}_dwi_corr.mif

if [[ ! -f $dwi_corr ]] && [[ ! -f ${dwi_n4} ]]; then
Info "DWI denoise, bias filed correction and concatenation"
# Concatenate shells -if only one shell then just convert to mif and rename.
      for dwi in ${bids_dwis[@]}; do
            dwi_nom=`echo $dwi | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}'`
            bids_dwi_str=`echo $dwi | awk -F . '{print $1}'`
            Do_cmd mrconvert $dwi -json_import ${bids_dwi_str}.json -fslgrad ${bids_dwi_str}.bvec ${bids_dwi_str}.bval ${tmp}/${dwi_nom}.mif
      done

      # Concatenate shells and convert to mif.
      if [ "${#bids_dwis[@]}" -eq 1 ]; then
        cp ${tmp}/*.mif $dwi_cat
      else
        Do_cmd mrcat ${tmp}/*.mif $dwi_cat -nthreads $CORES
      fi

      # Denoise DWI and calculate residuals
      Do_cmd dwidenoise $dwi_cat $dwi_dns -nthreads $CORES
      Do_cmd mrcalc $dwi_cat $dwi_dns -subtract ${proc_dwi}/${id}_dwi_residuals.mif -nthreads $CORES

      # Bias field correction DWI
      Do_cmd dwibiascorrect ants $dwi_dns $dwi_n4 -force -nthreads $CORES -scratch $tmp
      # Step QC
      if [[ -f ${dwi_n4} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has DWI in mif, denoised and concatenaded"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# dwifslpreproc and TOPUP preparations
if [[ ! -f $dwi_corr ]]; then
      Info "DWI dwifslpreproc"
      # Get parameters
      ReadoutTime=`mrinfo $dwi_n4 -property TotalReadoutTime`
      pe_dir=`mrinfo $dwi_n4 -property PhaseEncodingDirection`
      shells=(`mrinfo $dwi_n4 -shell_bvalues`)
      # Exclude shells with 0 value
      for i in "${!shells[@]}"; do if [[ ${shells[i]} = 0 ]]; then unset 'shells[i]'; fi; done

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
            opt="-rpe_pair -align_seepi -se_epi ${b0_pair}"
      else
            opt='-rpe_none'
      fi

      Info "dwifslpreproc parameters:"
      Note "Shell values        :" ${shells[@]}
      Note "DWI main dimensions :" "`mrinfo $dwi_4proc -size`"
      Note "DWI rpe dimensions  :" "`mrinfo $b0_pair -size`"
      Note "pe_dir              :" $pe_dir
      Note "Readout Time        :" $ReadoutTime

      # Preprocess each shell
      # DWIs all acquired with a single fixed phase encoding; but additionally a
      # pair of b=0 images with reversed phase encoding to estimate the inhomogeneity field:
      echo "COMMAND --> dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options " --data_is_shelled --slm=linear" -nthreads $CORES -nocleanup -scratch $tmp -force"
      dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options " --data_is_shelled --slm=linear" -nthreads $CORES -nocleanup -scratch $tmp -force
      # Step QC
      if [[ ! -f ${dwi_corr} ]]; then Error "dwifslpreproc failed, check the logs"; exit;
      else
          Do_cmd rm $dwi_n4; ((Nsteps++))
          # eddy_quad Quality Check
          eddy_quad $tmp/dwi_post_eddy -idx $tmp/eddy_indices.txt -par $tmp/eddy_config.txt -m $tmp/eddy_mask.nii -b $tmp/bvals -o ${proc_dwi}/eddy_QC
          # Copy eddy parameters
          eddy_DIR=${proc_dwi}/eddy
          if [ ! -d ${eddy_DIR} ]; then Do_cmd mkdir ${eddy_DIR}; fi
          Do_cmd cp -rf ${tmp}/dwifslpreproc*/*eddy.eddy* ${eddy_DIR}
      fi
else
      Info "Subject ${id} has a DWI processed with dwifslpreproc"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Registration of corrected DWI-b0 to T1nativepro
dwi_mask=${proc_dwi}/${id}_dwi_mask.nii.gz
dwi_5tt=${proc_dwi}/${id}_dwi_5tt.nii.gz
dwi_b0=${proc_dwi}/${id}_dwi_b0.nii.gz # This should be a NIFTI for compatibility with ANTS
dwi_in_T1nativepro=${proc_struct}/${id}_dwi_nativepro.nii.gz
T1nativepro_in_dwi=${proc_dwi}/${id}_t1w_dwi.nii.gz
str_dwi_affine=${dir_warp}/${id}_dwi_to_nativepro_
mat_dwi_affine=${str_dwi_affine}0GenericAffine.mat

if [[ ! -f $T1nativepro_in_dwi ]]; then
      Info "Registering DWI b0 to T1nativepro"
      # Corrected DWI-b0s mean for registration
      dwiextract -force -nthreads $CORES $dwi_corr - -bzero | mrmath - mean $dwi_b0 -axis 3

      # Register DWI-b0 mean corrected to T1nativepro
      Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro -m $dwi_b0 -o $str_dwi_affine -t a -n $CORES -p d
      # Apply transformation DWI-b0 space to T1nativepro
      Do_cmd antsApplyTransforms -d 3 -i $dwi_b0 -r $T1nativepro -t $mat_dwi_affine -o $dwi_in_T1nativepro -v -u int
      # Apply inverse transformation T1nativepro to DWI-b0 space
      Do_cmd antsApplyTransforms -d 3 -i $T1nativepro -r $dwi_b0 -t [$mat_dwi_affine,1] -o $T1nativepro_in_dwi -v -u int
      # Step QC
      if [[ -f ${T1nativepro_in_dwi} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has a T1nativepro in DWI-b0 space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
if [[ ! -f $dwi_mask ]]; then
      Info "Creating DWI binary mask of processed volumes"
      # Create a binary mask of the DWI
      Do_cmd antsApplyTransforms -d 3 -i $MNI152_mask \
                  -r ${dwi_b0} \
                  -n GenericLabel -t [$mat_dwi_affine,1] -t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp} \
                  -o ${tmp}/dwi_mask.nii.gz -v
      Do_cmd maskfilter ${tmp}/dwi_mask.nii.gz erode -npass 1 $dwi_mask
      # Step QC
      if [[ -f ${dwi_mask} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has a DWI-preproc binary mask"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# 5TT file in dwi space
if [[ ! -f $dwi_5tt ]]; then
      Info "Registering 5TT file to DWI-b0 space"
      Do_cmd antsApplyTransforms -d 3 -e 3 -i $T15ttgen -r ${dwi_b0} -n linear -t [$mat_dwi_affine,1] -o ${dwi_5tt} -v
      # Step QC
      if [[ -f ${dwi_5tt} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has a 5TT segmentation in DWI space"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Get some basic metrics.
dwi_FA=$proc_dwi/${id}_dti_FA.mif
if [[ ! -f $dwi_FA ]]; then
      Info "Calculating basic DTI metrics"
      dwi2tensor -mask $dwi_mask -nthreads $CORES $dwi_corr ${tmp}/${id}_dti.mif
      tensor2metric -nthreads $CORES -fa $proc_dwi/${id}_dti_FA.mif -adc $proc_dwi/${id}_dti_ADC.mif ${tmp}/${id}_dti.mif
      # Step QC
      if [[ -f ${dwi_FA} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has DWI tensor metrics"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Response function and Fiber Orientation Distribution
fod=$proc_dwi/${id}_wm_fod_norm.mif
fod_gmN=$proc_dwi/${id}_gm_fod_norm.mif
fod_csfN=$proc_dwi/${id}_csf_fod_norm.mif
if [[ ! -f $fod ]]; then
      Info "Calculating Multi-Shell Multi-Tissue, Response function and Fiber Orientation Distribution"
      # if [ "${#shells[@]}" -ge 2 ]; then
            rf=dhollander
            # Response function
            rf_wm=${tmp}/${id}_response_wm_${rf}.txt
            rf_gm=${tmp}/${id}_response_gm_${rf}.txt
            rf_csf=${tmp}/${id}_response_csf_${rf}.txt
            # Fiber Orientation Distriution
            fod_wm=${tmp}/${id}_wm_fod.mif
            fod_gm=${tmp}/${id}_gm_fod.mif
            fod_csf=${tmp}/${id}_csf_fod.mif

            Do_cmd dwi2response $rf -nthreads $CORES $dwi_corr $rf_wm $rf_gm $rf_csf -mask $dwi_mask
            Do_cmd dwi2fod -nthreads $CORES msmt_csd $dwi_corr \
                $rf_wm $fod_wm \
                $rf_gm $fod_gm \
                $rf_csf $fod_csf \
                -mask $dwi_mask
      if [ "${#shells[@]}" -ge 2 ]; then
            Do_cmd mtnormalise $fod_wm $fod $fod_gm $fod_gmN $fod_csf $fod_csfN -nthreads $CORES -mask $dwi_mask
      else
      #     Info "Calculating Single-Shell, Response function and Fiber Orientation Distribution"
      #
      #       Do_cmd dwi2response tournier -nthreads $CORES $dwi_corr \
      #           ${proc_dwi}/${id}_response_tournier.txt \
      #           -mask $dwi_mask
      #       Do_cmd dwi2fod -nthreads $CORES csd \
      #           $dwi_corr \
      #           ${proc_dwi}/${id}_response_tournier.txt  \
      #           ${tmp}/FOD.mif \
      #           -mask $dwi_mask
            Do_cmd mtnormalise -nthreads $CORES -mask $dwi_mask $fod_wm $fod
      fi
      # Step QC
      if [[ -f ${fod} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has Fiber Orientation Distribution files"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# QC of the tractography
tdi=${proc_dwi}/${id}_tdi_iFOD1-1M.mif
if [[ ! -f $tdi ]]; then
      tract=${tmp}/${id}_QC_iFOD1_1M.tck
      Do_cmd tckgen -nthreads $CORES \
          $fod \
          $tract \
          -act $dwi_5tt \
          -crop_at_gmwmi \
          -backtrack \
          -seed_dynamic $fod \
          -maxlength 300 \
          -angle 22.5 \
          -power 1.0 \
          -select 1M \
          -step .5 \
          -cutoff 0.06 \
          -algorithm iFOD1
      Do_cmd tckmap -vox 1,1,1 -dec -nthreads $CORES $tract $tdi
      # Step QC
      if [[ -f ${tdi} ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has a Tract Density Image for QC 1M streamlines"; ((Nsteps++))
fi

# -----------------------------------------------------------------------------------------------
# Clean temporal directory
cd $here
if [[ -z $nocleanup ]]; then Do_cmd rm -rf $tmp; else Info "tmp directory was not erased: ${tmp}"; fi

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
if [ "$Nsteps" -eq 8 ]; then status="COMPLETED"; else status="ERROR DWI is missing a processing step: "; fi
Title "DWI processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:
\t\tSteps completed: `printf "%02d" $Nsteps`/08
\tStatus          : $status
\tCheck logs:
`ls ${dir_logs}/proc-dwi_*.txt`"
# Print QC stamp
echo "${id}, proc_dwi, $status N=`printf "%02d" $Nsteps`/08, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
