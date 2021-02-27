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

BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
dwi_main=$8
dwi_rpe=$9
dwi_processed=${10}
PROC=${11}
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe
    source ${MICAPIPE}/functions/init.sh;
fi
# init file for local processing at MICA lab
if [ "$PROC" = "LOCAL-MICA" ]; then source ${MICAPIPE}/functions/init.sh; fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

Info "Inputs of proc_dwi:"
Note "tmpDir     :" "$tmpDir"
Note "dwi_main   :" "$dwi_main"
Note "dwi_rpe    :" "$dwi_rpe"
Note "Processing :" "$PROC"

# Manage manual inputs: DWI main image(s)
if [[ "$dwi_main" != "DEFAULT" ]]; then
    IFS=',' read -ra bids_dwis <<< "$dwi_main"
fi
# Manage manual inputs: DWI reverse phase encoding
if [[ "$dwi_rpe" != "DEFAULT" ]]; then dwi_reverse=$dwi_rpe; fi

# Manage manual inputs: DWI pre-processed
if [[ "$dwi_processed" != "FALSE" ]]; then
    if [[ ! -f "$dwi_processed" ]]; then Error "Provided dwi_processed image's path is wrong of file doesn't exist!!, check:\n\tls $dwi_processed"; exit; fi
    # Check mif format
    if [ ${dwi_processed: -4} != ".mif" ]; then Error "Provided dwi_processed image is not in mif format!!"; exit; fi
    # Check diffusion-weighting gradient table, as interpreted by MRtrix3
    Ndir=$(mrinfo -dwgrad "$dwi_processed" | wc -l)
    if [ "${Ndir}" -lt 12 ]; then Error "Provided dwi_processed image doesn't have a dwi gradient (dir=$Ndir) table encoded or enough directions (minimum 12)!!"; exit; fi
    Npar=$(mrinfo -dwgrad "$dwi_processed" | awk '{print NF}' | sort -nu | tail -n 1)
    if [ "${Npar}" -lt 4 ]; then Error "Provided dwi_processed image doesn't have enough number of parameters encoded (x,y,z,bval should be 4 is  $Npar)!!"; exit; fi
    Info "Provided dwi_processed with $Ndir encoded directions seems ok"
    dwi_corr="$dwi_processed"
fi

# Check inputs: DWI
if [ "${#bids_dwis[@]}" -lt 1 ]; then Error "Subject $id doesn't have DWIs:\n\t\t TRY <ls -l ${subject_bids}/dwi/>"; exit; fi
if [ ! -f "${T1_MNI152_InvWarp}" ]; then Error "Subject $id doesn't have T1_nativepro warp to MNI152.\n\t\tRun -proc_structural"; exit; fi
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro.\n\t\tRun -proc_structural"; exit; fi
if [ ! -f ${T15ttgen} ]; then Error "Subject $id doesn't have a 5tt volume in nativepro space.\n\t\tRun -proc_structural"; exit; fi

# CHECK if PhaseEncodingDirection and TotalReadoutTime exist
for i in ${bids_dwis[*]}; do
  json=$(echo $i awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1 ".json"}')
  ped=$(cat $json | grep PhaseEncodingDirection\": | awk -F "\"" '{print $4}')
  trt=$(cat $json | grep TotalReadoutTime | awk -F " " '{print $2}')
  if [[ -z "$ped" ]]; then Error "PhaseEncodingDirection is missing in $json"; exit; fi
  if [[ -z "$trt" ]]; then Error "TotalReadoutTime is missing in $json"; exit; fi
done

#------------------------------------------------------------------------------#
Title "Diffusion Weighted Imaging processing\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-dwi
Info "Saving temporal dir: $nocleanup"
Info "ANTs and MRtrix will use $threads threads"

#	Timer
aloita=$(date +%s)
Nsteps=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_proc-dwi_${id}"
Do_cmd mkdir -p $tmp

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

Do_cmd cd "$tmp"
#------------------------------------------------------------------------------#
# DWI processing
# Image denoising must be performed as the first step of the image-processing pipeline.
# Interpolation or smoothing in other processing steps, such as motion and distortion correction,
# may alter the noise characteristics and thus violate the assumptions upon which MP-PCA is based.
dwi_cat="${tmp}/dwi_concatenate.mif"
dwi_dns="${tmp}/${idBIDS}_space-dwi_desc-MP-PCA_dwi.mif"
dwi_n4="${proc_dwi}/${idBIDS}_space-dwi_desc-MP-PCA_N4_dwi.mif"
dwi_res="${proc_dwi}/${idBIDS}_space-dwi_desc-MP-PCA_residuals-dwi.mif"
dwi_corr="${proc_dwi}/${id}_dwi_corr.mif"
b0_refacq=$(echo ${bids_dwis[0]} | awk -F 'acq-' '{print $2}'| sed 's:_dwi.nii.gz::g')

if [[ "$dwi_processed" == "FALSE" ]]; then
    if [ ! -f "$dwi_res" ] || [ ! -f "$dwi_n4" ]; then
    Info "DWI denoise, bias filed correction and concatenation"
    # Concatenate shells -if only one shell then just convert to mif and rename.
          for dwi in ${bids_dwis[@]}; do
                dwi_nom=$(echo $dwi | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')
                bids_dwi_str=$(echo $dwi | awk -F . '{print $1}')
                Do_cmd mrconvert $dwi -json_import ${bids_dwi_str}.json -fslgrad ${bids_dwi_str}.bvec ${bids_dwi_str}.bval ${tmp}/${dwi_nom}.mif
                dwiextract ${tmp}/${dwi_nom}.mif - -bzero | mrmath - mean "${tmp}/${dwi_nom}_b0.nii.gz" -axis 3 -nthreads $threads
          done

          # Rigid registration between shells
          n=$((${#bids_dwis[*]} - 1))
          if [[ ${#bids_dwis[*]} -gt 1 ]]; then
            b0_ref=${tmp}/$(echo "${bids_dwis[0]}" | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')_b0.nii.gz
            for ((i=1; i<=$n; i++)); do
                dwi_nom=$(echo ${bids_dwis[i]} | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')
                bids_dwi_str=$(echo ${bids_dwis[i]} | awk -F . '{print $1}')
                b0_acq=$(echo ${bids_dwis[i]} | awk -F 'acq-' '{print $2}'| sed 's:_dwi.nii.gz::g')
                b0_nom="${tmp}/$(echo ${bids_dwis[i]} | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')_b0.nii.gz"
                b0_run="acq-${b0_acq}"
                b0mat_str="${dir_warp}/${idBIDS}_from-${b0_acq}_to-${b0_refacq}_mode-image_desc-rigid_"
                b0mat="${b0mat_str}0GenericAffine.mat"
                b0run_2_b0ref="${tmp}/${idBIDS}_from-${b0_acq}_to-${b0_refacq}.nii.gz"

                Info "Registering ${b0_acq} to ${b0_refacq}"
                Do_cmd antsRegistrationSyN.sh -d 3 -m "$b0_nom" -f "$b0_ref"  -o "$b0mat_str" -t r -n "$threads" -p d
                mrconvert ${tmp}/${dwi_nom}.mif ${tmp}/${dwi_nom}.nii.gz
                Do_cmd antsApplyTransforms -d 3 -e 3 -i "${tmp}/${dwi_nom}.nii.gz" -r "$b0_ref" -t "$b0mat" -o "${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz" -v -u int
                Do_cmd mrconvert ${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz -json_import ${bids_dwi_str}.json -fslgrad ${bids_dwi_str}.bvec ${bids_dwi_str}.bval ${tmp}/${dwi_nom}_Ralign.mif -force -quiet
            done
          fi

          Info "Concatenatenating shells"
          dwi_0=$(echo ${bids_dwis[0]} | awk -F "dwi/" '{print $2}' | awk -F ".nii" '{print $1}')
          if [ "${#bids_dwis[@]}" -eq 1 ]; then
            cp "${tmp}/${dwi_0}.mif" "$dwi_cat"
          else
            Do_cmd mrcat "${tmp}/${dwi_0}.mif" "${tmp}/*_Ralign.mif" "$dwi_cat" -nthreads $threads
          fi

          # Denoise DWI and calculate residuals
          Do_cmd dwidenoise $dwi_cat $dwi_dns -nthreads $threads
          Do_cmd mrcalc $dwi_cat $dwi_dns -subtract $dwi_res -nthreads $threads

          # Bias field correction DWI
          Do_cmd dwibiascorrect ants $dwi_dns $dwi_n4 -force -nthreads $threads -scratch $tmp
          # Step QC
          if [[ -f "$dwi_n4" ]]; then ((Nsteps++)); fi
    else
          Info "Subject ${id} has DWI in mif, denoised and concatenaded"; ((Nsteps++))
    fi
else
    Info "Subject ${id} has a DWI processed, skipping denoise and bias field correction"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# dwifslpreproc and TOPUP preparations
if [[ ! -f "$dwi_corr" ]]; then
      Info "DWI dwifslpreproc"
      # Get parameters
      ReadoutTime=$(mrinfo $dwi_n4 -property TotalReadoutTime)
      pe_dir=$(mrinfo $dwi_n4 -property PhaseEncodingDirection)
      shells=$(mrinfo $dwi_n4 -shell_bvalues)
      # Exclude shells with 0 value
      for i in "${!shells[@]}"; do if [[ ${shells[i]} = 0 ]]; then unset 'shells[i]'; fi; done

      # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
      dwi_4proc=${tmp}/dwi_n4_even.mif
      dim=$(mrinfo $dwi_n4 -size)
      dimNew=($(echo $dim | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
      mrconvert $dwi_n4 $dwi_4proc -coord 0 0:${dimNew[0]} -coord 1 0:${dimNew[1]} -coord 2 0:${dimNew[2]} -coord 3 0:end -force

      # Get the mean b-zero (un-corrected)
      dwiextract -nthreads $threads $dwi_n4 - -bzero | mrmath - mean ${tmp}/b0_meanMainPhase.mif -axis 3

      # Processing the reverse encoding b0
      Title "$dwi_reverse"
      if [ -f $dwi_reverse ]; then
            b0_pair_tmp="${tmp}/b0_pair_tmp.mif"
            b0_pair="${tmp}/b0_pair.mif"
            dwi_reverse_str=$(echo $dwi_reverse | awk -F . '{print $1}')

            Do_cmd mrconvert $dwi_reverse -json_import ${dwi_reverse_str}.json -fslgrad ${dwi_reverse_str}.bvec ${dwi_reverse_str}.bval ${tmp}/b0_ReversePhase.mif
            dwiextract ${tmp}/b0_ReversePhase.mif - -bzero | mrmath - mean ${tmp}/b0_meanReversePhase.nii.gz -axis 3 -nthreads $threads

            # Linear registration between both b0
            rpe=$(echo ${dwi_reverse} | awk -F 'acq-' '{print $2}'| sed 's:_dwi.nii.gz::g')
            rpemat_str="${dir_warp}/${idBIDS}_from-${rpe}_to-${b0_refacq}_mode-image_desc-rigid_"
            rpemat="${rpemat_str}0GenericAffine.mat"
            Do_cmd mrconvert ${tmp}/b0_meanMainPhase.mif ${tmp}/b0_meanMainPhase.nii.gz
            Do_cmd antsRegistrationSyN.sh -d 3 -m "${tmp}/b0_meanReversePhase.nii.gz" -f "${tmp}/b0_meanMainPhase.nii.gz"  -o "$rpemat_str" -t r -n "$threads" -p d
            Do_cmd antsApplyTransforms -d 3 -e 3 -i "${dwi_reverse}" -r "${tmp}/b0_meanMainPhase.nii.gz" -t "$rpemat" -o "${tmp}/b0_meanReversePhase-reg.nii.gz" -v -u int
            Do_cmd mrconvert "${tmp}/b0_meanReversePhase-reg.nii.gz" -json_import ${dwi_reverse_str}.json -fslgrad ${dwi_reverse_str}.bvec ${dwi_reverse_str}.bval ${tmp}/b0_ReversePhase.mif -force -quiet
            dwiextract ${tmp}/b0_ReversePhase.mif - -bzero | mrmath - mean ${tmp}/b0_meanReversePhase.mif -axis 3 -nthreads $threads

            # Concatenate the pe and rpe b0s
            Do_cmd mrcat "${tmp}/b0_meanMainPhase.mif" "${tmp}/b0_meanReversePhase.mif" "$b0_pair_tmp" -nthreads "$threads"

            # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
            dim=$(mrinfo $b0_pair_tmp -size)
            dimNew=($(echo $dim | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
            mrconvert $b0_pair_tmp $b0_pair -coord 0 0:${dimNew[0]} -coord 1 0:${dimNew[1]} -coord 2 0:${dimNew[2]} -coord 3 0:end -force
            opt="-rpe_pair -align_seepi -se_epi ${b0_pair}"
      else
            Info "Reverse phase encoding image was not found it will be omitted"
            opt='-rpe_none'
      fi

      Info "dwifslpreproc parameters:"
      Note "Shell values        :" "${shells[*]}"
      Note "DWI main dimensions :" "$(mrinfo $dwi_4proc -size)"
      if [ -f $dwi_reverse ]; then Note "DWI rpe dimensions  :" "$(mrinfo $b0_pair -size)"; fi
      Note "pe_dir              :" $pe_dir
      Note "Readout Time        :" $ReadoutTime

      # Preprocess each shell
      # DWIs all acquired with a single fixed phase encoding; but additionally a
      # pair of b=0 images with reversed phase encoding to estimate the inhomogeneity field:
      echo -e "COMMAND --> dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options \" --data_is_shelled --slm=linear\" -nthreads $threads -nocleanup -scratch $tmp -force"
      dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options " --data_is_shelled --slm=linear" -nthreads $threads -nocleanup -scratch $tmp -force
      # Step QC
      if [[ ! -f "$dwi_corr" ]]; then Error "dwifslpreproc failed, check the logs"; exit;
      else
          Do_cmd rm $dwi_n4; ((Nsteps++))
          # eddy_quad Quality Check
          Do_cmd cd $tmp/dwifslpreproc*
          Do_cmd eddy_quad dwi_post_eddy -idx eddy_indices.txt -par eddy_config.txt -m eddy_mask.nii -b bvals -o ${dir_QC}/eddy_QC
          Do_cmd cd $tmp
          # Copy eddy parameters
          eddy_DIR=${proc_dwi}/eddy
          if [ ! -d ${eddy_DIR} ]; then Do_cmd mkdir ${eddy_DIR}; fi
          Do_cmd cp -rf ${tmp}/dwifslpreproc*/*eddy.eddy* ${eddy_DIR}
          Do_cmd chmod 770 -R ${eddy_DIR}/*
      fi
else
      Info "Subject ${id} has a DWI processed"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Registration of corrected DWI-b0 to T1nativepro
dwi_mask="${proc_dwi}/${idBIDS}_space-dwi_desc-brain_mask.nii.gz"
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0.nii.gz" # This should be a NIFTI for compatibility with ANTS
str_dwi_affine="${dir_warp}/${idBIDS}_space-dwi_from-dwi_to-nativepro_mode-image_desc-"
mat_dwi_affine="${str_dwi_affine}0GenericAffine.mat"
T1nativepro_in_dwi="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro.nii.gz"

if [[ ! -f "$T1nativepro_in_dwi" ]]; then
      Info "Linear registration from DWI-b0 to T1nativepro"
      # Corrected DWI-b0s mean for registration
      dwiextract -force -nthreads $threads $dwi_corr - -bzero | mrmath - mean $dwi_b0 -axis 3 -force

      # Register DWI-b0 mean corrected to T1nativepro
      Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro_brain -m $dwi_b0 -o $str_dwi_affine -t a -n $threads -p d
      # Apply inverse transformation T1nativepro to DWI-b0 space
      Do_cmd antsApplyTransforms -d 3 -i $T1nativepro -r $dwi_b0 -t [$mat_dwi_affine,1] -o $T1nativepro_in_dwi -v -u int
      if [[ -f "$T1nativepro_in_dwi" ]]; then ((Nsteps++)); fi

      #------------------------------------------------------------------------------#
      Info "Creating DWI binary mask of processed volumes"
      # Create a binary mask of the DWI
      Do_cmd antsApplyTransforms -d 3 -i $MNI152_mask \
              -r ${dwi_b0} \
              -n GenericLabel -t [$mat_dwi_affine,1] -t [${T1_MNI152_affine},1] -t ${T1_MNI152_InvWarp} \
              -o ${tmp}/dwi_mask.nii.gz -v
      Do_cmd maskfilter ${tmp}/dwi_mask.nii.gz erode -npass 1 $dwi_mask
      if [[ -f "$dwi_mask" ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has an affine transformation from T1w to DWI-b0 space"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
# Get some basic metrics.
dwi_dti="${proc_dwi}/${idBIDS}_space-dwi_model-DTI.mif"
dti_FA="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.mif"
dti_ADC="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-ADC.mif"
if [[ ! -f "$dti_FA" ]]; then
      Info "Calculating basic DTI metrics"
      dwi2tensor -mask $dwi_mask -nthreads $threads $dwi_corr $dwi_dti
      tensor2metric -nthreads $threads -fa ${dti_FA} -adc ${dti_ADC} $dwi_dti
      if [[ -f "$dti_FA" ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has diffusion tensor metrics"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Response function and Fiber Orientation Distribution
fod_wmN="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.mif"
fod_gmN="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-gmNorm.mif"
fod_csfN="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-csfNorm.mif"
if [[ ! -f "$fod_wmN" ]]; then
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

            Do_cmd dwi2response $rf -nthreads $threads $dwi_corr $rf_wm $rf_gm $rf_csf -mask $dwi_mask
            Do_cmd dwi2fod -nthreads $threads msmt_csd $dwi_corr \
                $rf_wm $fod_wm \
                $rf_gm $fod_gm \
                $rf_csf $fod_csf \
                -mask $dwi_mask
      if [ "${#shells[@]}" -ge 2 ]; then
            Do_cmd mtnormalise $fod_wm $fod_wmN $fod_gm $fod_gmN $fod_csf $fod_csfN -nthreads $threads -mask $dwi_mask
      else
      #     Info "Calculating Single-Shell, Response function and Fiber Orientation Distribution"
      #
      #       Do_cmd dwi2response tournier -nthreads $threads $dwi_corr \
      #           ${proc_dwi}/${id}_response_tournier.txt \
      #           -mask $dwi_mask
      #       Do_cmd dwi2fod -nthreads $threads csd \
      #           $dwi_corr \
      #           ${proc_dwi}/${id}_response_tournier.txt  \
      #           ${tmp}/FOD.mif \
      #           -mask $dwi_mask
            Do_cmd mtnormalise -nthreads $threads -mask $dwi_mask $fod_wm $fod_wmN
      fi
      if [[ -f "$fod_wmN" ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has Fiber Orientation Distribution file"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Non-linear registration between DWI space and T1w
dwi_SyN_str="${dir_warp}/${idBIDS}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_"
dwi_SyN_warp="${dwi_SyN_str}1Warp.nii.gz"
dwi_SyN_Invwarp="${dwi_SyN_str}1InverseWarp.nii.gz"
dwi_SyN_affine="${dwi_SyN_str}0GenericAffine.mat"
dwi_5tt="${proc_dwi}/${idBIDS}_space-dwi_desc-5tt.nii.gz"

if [[ ! -f "$dwi_SyN_warp" ]]; then
    Info "Non-linear registration from T1w_dwi-space to DWI"
    dwi_in_T1nativepro="${proc_struct}/${idBIDS}_space-nativepro_desc-dwi.nii.gz" # Only for QC
    T1nativepro_in_dwi_brain="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro-brain.nii.gz"
    T1nativepro_in_dwi_NL="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro_NL.nii.gz"
    fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"

    Do_cmd fslmaths $T1nativepro_in_dwi -mul $dwi_mask $T1nativepro_in_dwi_brain
    Do_cmd mrconvert -coord 3 0 $fod_wmN $fod
    Do_cmd antsRegistrationSyN.sh -d 3 -m $T1nativepro_in_dwi_brain -f $fod -o $dwi_SyN_str -t s -n $threads
    if [[ -f "$dwi_SyN_warp" ]]; then ((Nsteps++)); fi
    Info "Registering T1w-nativepro and 5TT to DWI-b0 space, and DWI-b0 to T1w-nativepro"
    # Apply transformation DWI-b0 space to T1nativepro
    Do_cmd antsApplyTransforms -d 3 -r $T1nativepro_brain -i $dwi_b0 -r $T1nativepro_brain -t $mat_dwi_affine -t [$dwi_SyN_affine,1] -t $dwi_SyN_Invwarp -o $dwi_in_T1nativepro -v -u int
    # Apply transformation T1nativepro to DWI space
    Do_cmd antsApplyTransforms -d 3 -r $fod -i $T1nativepro -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $T1nativepro_in_dwi_NL -v -u int
    # Apply transformation 5TT to DWI space
    Do_cmd antsApplyTransforms -d 3 -r $fod -i $T15ttgen -t $dwi_SyN_warp -t $dwi_SyN_affine -t [$mat_dwi_affine,1] -o $dwi_5tt -v -e 3 -n linear
    if [[ -f "$dwi_5tt" ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a non-linear registration from T1w_dwi-space to DWI"; Nsteps=$((Nsteps + 2))
fi
#------------------------------------------------------------------------------#
# Gray matter White matter interface mask
dwi_gmwmi="${proc_dwi}/${idBIDS}_space-dwi_desc-gmwmi-mask.mif"
if [[ ! -f $dwi_gmwmi ]]; then
      Info "Calculating Gray matter White matter interface mask"
      5tt2gmwmi $dwi_5tt $dwi_gmwmi; ((Nsteps++))
else
      Info "Subject ${id} has Gray matter White matter interface mask"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# QC of the tractography
tracts=1M
tdi_1M="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD1-${tracts}_tdi.mif"
tckjson="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD1-${tracts}_tractography.json"
if [[ ! -f "$tdi_1M" ]]; then
  Info "Creating a track density image for quality check"
  tck_1M="${tmp}/${idBIDS}_space-dwi_desc-iFOD1-${tracts}_tractography.tck"
  Do_cmd tckgen -nthreads $threads \
      $fod_wmN \
      $tck_1M \
      -act $dwi_5tt \
      -crop_at_gmwmi \
      -backtrack \
      -seed_gmwmi $dwi_gmwmi \
      -algorithm iFOD1 \
      -step 0.5 \
      -angle 22.5 \
      -cutoff 0.06 \
      -maxlength 400 \
      -minlength 10 \
      -power 1.0 \
      -select ${tracts}

  Do_cmd tckmap -vox 1,1,1 -dec -nthreads $threads $tck_1M $tdi_1M
  if [[ -f "$tdi_1M" ]]; then ((Nsteps++)); fi
  tck_json iFOD1 0.5 22.5 0.06 400 10 seed_gmwmi $tck_1M
else
      Info "Subject ${id} has a Tract Density Image for QC 1M streamlines"; ((Nsteps++))
fi

# -----------------------------------------------------------------------------------------------
# QC: Input files
QC_proc-dwi

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print $eri/60 | perl)

# Notification of completition
if [ "$Nsteps" -eq 10 ]; then status="COMPLETED"; else status="ERROR DWI is missing a processing step: "; fi
Title "DWI processing ended in \033[38;5;220m $(printf "%0.3f\n" ${eri}) minutes \033[38;5;141m:
\t\tSteps completed: $(printf "%02d" $Nsteps)/10
\tStatus          : $status
\tCheck logs:
$(ls ${dir_logs}/proc-dwi_*.txt)"
# Print QC stamp
echo "${id}, proc_dwi, $status N=$(printf "%02d" $Nsteps)/10, $(whoami), $(uname -n), $(date), $(printf "%0.3f\n" ${eri}), $PROC" >> ${out}/brain-proc.csv
cleanup $tmp $nocleanup $here
