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
umask 003
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
rpe_all=${11}
regAffine=${12}
dwi_str=${13}
PROC=${14}
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/host/yeatman/local_raid/rcruces/git_here/micapipe
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# mtrix configuration file
mrconf="${HOME}/.mrtrix.conf"
echo "BZeroThreshold: 61" > "$mrconf"

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

Info "Inputs of proc_dwi:"
Note "tmpDir     : " "$tmpDir"
Note "dwi_main   : " "$dwi_main"
Note "dwi_rpe    : " "$dwi_rpe"
Note "rpe_all    : " "$rpe_all"
Note "dwi_acq    : " "$dwi_str"
Note "Affine only: " "$regAffine"
Note "Processing : " "$PROC"

# Manage manual inputs: DWI main image(s)
if [[ "$dwi_main" != "DEFAULT" ]]; then
    IFS=',' read -ra bids_dwis <<< "$dwi_main"
    bids_dwis=("${bids_dwis[@]}")
fi
# Manage manual inputs: DWI reverse phase encoding
if [[ "$dwi_rpe" != "DEFAULT" ]]; then
    IFS=',' read -ra dwi_reverse <<< "$dwi_rpe"
    dwi_reverse=("${dwi_reverse[@]}")
fi

# Manage manual inputs: DWI pre-processed
if [[ "$dwi_processed" != "FALSE" ]]; then
    if [[ ! -f "$dwi_processed" ]]; then Error "Provided dwi_processed image's path is wrong of file doesn't exist!!, check:\n\tls $dwi_processed"; exit; fi
    # Check mif format
    if [ "${dwi_processed: -4}" != ".mif" ]; then Error "Provided dwi_processed image is not in mif format!!"; exit; fi
    # Check diffusion-weighting gradient table, as interpreted by MRtrix3
    Ndir=$(mrinfo -dwgrad "$dwi_processed" | wc -l)
    if [ "${Ndir}" -lt 12 ]; then Error "Provided dwi_processed image doesn't have a dwi gradient (dir=$Ndir) table encoded or enough directions (minimum 12)!!"; exit; fi
    Npar=$(mrinfo -dwgrad "$dwi_processed" | awk '{print NF}' | sort -nu | tail -n 1)
    if [ "${Npar}" -lt 4 ]; then Error "Provided dwi_processed image doesn't have enough number of parameters encoded (x,y,z,bval should be 4 is  $Npar)!!"; exit; fi
    Info "Provided dwi_processed with $Ndir encoded directions seems ok"
    dwi_corr="$dwi_processed"
fi

# Check inputs: DWI
if [ "${#bids_dwis[@]}" -lt 1 ]; then Error "DWI string or path does not match the default:\n\t\t TRY to set it with -dwi_main <path to DWI.nii.gz>"; exit; fi
if [ ! -f "${bids_dwis[0]}" ]; then Error "Main DWI was not found:\n\t\t TRY to set it with -dwi_main <path to DWI.nii>"; exit; fi
if [ ! -f "${T1_MNI152_InvWarp}" ]; then Error "Subject $id doesn't have T1_nativepro warp to MNI152.\n\t\tRun -proc_structural"; exit; fi
if [ ! -f "${T1nativepro}" ]; then Error "Subject $id doesn't have T1_nativepro.\n\t\tRun -proc_structural"; exit; fi
if [ ! -f "${T15ttgen}" ]; then Error "Subject $id doesn't have a 5tt volume in nativepro space.\n\t\tRun -proc_structural"; exit; fi

# CHECK if PhaseEncodingDirection and TotalReadoutTime exist
for i in ${bids_dwis[*]}; do
  json=$(echo "${i}" | awk -F ".nii" '{print $1 ".json"}')
  ped=$(grep PhaseEncodingDirection\": "${json}" | awk -F "\"" '{print $4}')
  trt=$(grep TotalReadoutTime "${json}" | awk -F " " '{print $2}')
  if [[ -z "$ped" ]]; then Error "PhaseEncodingDirection is missing in $json"; exit; fi
  if [[ -z "$trt" ]]; then Error "TotalReadoutTime is missing in $json"; exit; fi
done

# Update path for multiple acquisitions processing
if [[ "${dwi_str}" != "DEFAULT" ]]; then
  dwi_str="acq-${dwi_str/acq-/}"
  dwi_str_="_${dwi_str}"
  export proc_dwi=$subject_dir/dwi/"${dwi_str}"
  [[ ! -d "$proc_dwi" ]] && Do_cmd mkdir -p "$proc_dwi" && chmod -R 770 "$proc_dwi"
else
  dwi_str=""; dwi_str_=""
fi

#------------------------------------------------------------------------------#
Title "Diffusion Weighted Imaging processing\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-dwi
Note "Saving temporal dir     :" "$nocleanup"
Note "ANTs and MRtrix will use: " "$threads threads"

#	Timer
aloita=$(date +%s)
Nsteps=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_proc-dwi_${id}"
Do_cmd mkdir -p "$tmp"
[[ ! -d "$dir_QC_png" ]] && Do_cmd mkdir -p "$dir_QC_png" && chmod -R 770 "$dir_QC_png"

# TRAP in case the script fails
trap 'rm $mrconf; cleanup $tmp $nocleanup $here' SIGINT SIGTERM

Do_cmd cd "$tmp"
#------------------------------------------------------------------------------#
# DWI processing
# Image denoising must be performed as the first step of the image-processing pipeline.
# Interpolation or smoothing in other processing steps, such as motion and distortion correction,
# may alter the noise characteristics and thus violate the assumptions upon which MP-PCA is based.
dwi_cat="${tmp}/dwi_concatenate.mif"
dwi_dns="${proc_dwi}/${idBIDS}_space-dwi_desc-MP-PCA_dwi.mif"
dwi_res="${proc_dwi}/${idBIDS}_space-dwi_desc-MP-PCA_residuals-dwi.mif"
dwi_resGibss="${proc_dwi}/${idBIDS}_space-dwi_desc-deGibbs_residuals-dwi.mif"
dwi_corr="${proc_dwi}/${idBIDS}_space-dwi_desc-dwi_preproc.mif"
b0_refacq=$(echo "${bids_dwis[0]##*/}" | awk -F ".nii" '{print $1}'); b0_refacq=$(echo "${b0_refacq/_dwi/}"); b0_refacq=$(echo "${b0_refacq/${idBIDS}_/}")

#------------------------------------------------------------------------------#
# Main phase encoding processing
if [[ "$dwi_processed" == "FALSE" ]] && [[ ! -f "$dwi_corr" ]]; then
    if [ ! -f "$dwi_res" ] || [ ! -f "$dwi_dns" ]; then
    Info "DWI denoise and concatenation"
    # Concatenate shells -if only one shell then just convert to mif and rename.
          for dwi in "${bids_dwis[@]}"; do
                dwi_nom=$(echo "${dwi##*/}" | awk -F ".nii" '{print $1}')
                bids_dwi_str=$(echo "$dwi" | awk -F . '{print $1}')
                Do_cmd mrconvert "$dwi" -json_import "${bids_dwi_str}.json" -fslgrad "${bids_dwi_str}.bvec" "${bids_dwi_str}.bval" "${tmp}/${dwi_nom}.mif"
                Do_cmd dwiextract "${tmp}/${dwi_nom}.mif" "${tmp}/${dwi_nom}_b0.mif" -bzero
                Do_cmd mrmath "${tmp}/${dwi_nom}_b0.mif" mean "${tmp}/${dwi_nom}_b0.nii.gz" -axis 3 -nthreads "$threads"
          done

          # Rigid registration between shells
          n=$((${#bids_dwis[*]} - 1))
          if [[ ${#bids_dwis[*]} -gt 1 ]]; then
            b0_ref=${tmp}/$(echo "${bids_dwis[0]##*/}" | awk -F ".nii" '{print $1}')_b0.nii.gz
            for ((i=1; i<=n; i++)); do
                dwi_nom=$(echo "${bids_dwis[i]##*/}" | awk -F ".nii" '{print $1}')
                bids_dwi_str=$(echo "${bids_dwis[i]}" | awk -F . '{print $1}')
                b0_nom="${tmp}/$(echo "${bids_dwis[i]##*/}" | awk -F ".nii" '{print $1}')_b0.nii.gz"
                b0_acq=$(echo "${dwi_nom/${idBIDS}_/}"); b0_acq=$(echo "${b0_acq/_dwi/}")
                b0mat_str="${dir_warp}/${idBIDS}_from-${b0_acq}_to-${b0_refacq}${dwi_str_}_mode-image_desc-rigid_"
                b0mat="${b0mat_str}0GenericAffine.mat"

                Info "Registering ${b0_acq} to ${b0_refacq}"
                Do_cmd antsRegistrationSyN.sh -d 3 -m "$b0_nom" -f "$b0_ref" -o "$b0mat_str" -t r -n "$threads" -p d
                mrconvert "${tmp}/${dwi_nom}.mif" "${tmp}/${dwi_nom}.nii.gz"
                Do_cmd antsApplyTransforms -d 3 -e 3 -i "${tmp}/${dwi_nom}.nii.gz" -r "$b0_ref" -t "$b0mat" -o "${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz" -v -u int
                Do_cmd mrconvert "${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz" -json_import "${bids_dwi_str}.json" -fslgrad "${bids_dwi_str}.bvec" "${bids_dwi_str}.bval" "${tmp}/${dwi_nom}__Ralign.mif" -force -quiet
            done
          fi

          Info "Concatenatenating shells"
          dwi_0=$(echo "${bids_dwis[0]##*/}" | awk -F ".nii" '{print $1}')
          if [ "${#bids_dwis[@]}" -eq 1 ]; then
            cp "${tmp}/${dwi_0}.mif" "$dwi_cat"
          else
            Do_cmd mrcat "${tmp}/${dwi_0}.mif" "${tmp}/*_Ralign.mif" "$dwi_cat" -nthreads "$threads"
          fi

          # Denoise DWI and calculate residuals
          Info "DWI MP-PCA denoising and Gibbs ringing correction"
          dwi_dns_tmp="${tmp}/MP-PCA_dwi.mif"
          Do_cmd dwidenoise "$dwi_cat" "$dwi_dns_tmp" -nthreads "$threads"
          mrcalc "$dwi_cat" "$dwi_dns_tmp" -subtract - -nthreads "$threads" | mrmath - mean "$dwi_res" -axis 3
          Do_cmd mrdegibbs "$dwi_dns_tmp" "$dwi_dns" -nthreads "$threads"
          mrcalc "$dwi_dns_tmp" "$dwi_dns" -subtract - -nthreads "$threads" | mrmath - mean "$dwi_resGibss" -axis 3
          Do_cmd mrinfo "$dwi_dns" -json_all "${dwi_dns/mif/json}"
          Do_cmd mrinfo "$dwi_res" -json_all "${dwi_res/mif/json}"
          Do_cmd mrinfo "$dwi_resGibss" -json_all "${dwi_resGibss/mif/json}"
          ((Nsteps++))
    else
          Info "Subject ${id} has DWI in mif, denoised and concatenaded"; ((Nsteps++))
    fi
else
    Info "Subject ${id} has a DWI processed, skipping denoise and concatenation"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
rpe_cat="${tmp}/dwi_rpe_concatenate.mif"
rpe_dns="${proc_dwi}/${idBIDS}_space-dwi_dir-rpe_desc-MP-PCA_dwi.mif"
if [[ "$dwi_processed" == "FALSE" ]] && [[ ! -f "$dwi_corr" ]]; then
    if [ ! -f "$rpe_dns" ] && [[ "${#dwi_reverse[@]}" -gt 0 ]]; then
    # Concatenate shells -if only one shell then just convert to mif and rename.
          for dwi in "${dwi_reverse[@]}"; do
                dwi_nom=$(echo "${dwi##*/}" | awk -F ".nii" '{print $1}')
                bids_dwi_str=$(echo "$dwi" | awk -F . '{print $1}')
                Info "DWI reverse phase encoding processing: $dwi_nom"
                if [[ -f "${bids_dwi_str}.bvec" ]] && [[ -f "${bids_dwi_str}.bval" ]]; then
                    Do_cmd mrconvert "$dwi" -json_import "${bids_dwi_str}.json" -fslgrad "${bids_dwi_str}.bvec" "${bids_dwi_str}.bval" "${tmp}/${dwi_nom}.mif"
                    Do_cmd dwiextract "${tmp}/${dwi_nom}.mif" "${tmp}/${dwi_nom}_b0.mif" -bzero
                else
                    Warning "No bval or bvecs were found the script will assumme that all the volumes are b0s!!!!"
                    Do_cmd mrconvert "$dwi" -json_import "${bids_dwi_str}.json" "${tmp}/${dwi_nom}.mif"
                    Do_cmd cp "${tmp}/${dwi_nom}.mif" "${tmp}/${dwi_nom}_b0.mif"
                fi

                # Get the reverse phase encoding mean b0
                rpe_dim=$(mrinfo "${tmp}/${dwi_nom}_b0.mif" -ndim)
                if [[ "$rpe_dim" -eq 3 ]]; then
                    Do_cmd mrconvert "${tmp}/${dwi_nom}.mif" "${tmp}/${dwi_nom}_b0.nii.gz"
                elif [[ "$rpe_dim" -gt 3 ]]; then
                    Do_cmd mrmath "${tmp}/${dwi_nom}_b0.mif" mean "${tmp}/${dwi_nom}_b0.nii.gz" -axis 3 -nthreads "$threads"
                fi
          done

          # Rigid registration between shells
          if [[ ${#dwi_reverse[*]} -gt 1 ]]; then
            b0_ref=${tmp}/dwi_pe_b0-mean.nii.gz
            dwiextract -nthreads "$threads" "$dwi_dns" - -bzero | mrmath - mean "$b0_ref" -axis 3
            for i in "${!dwi_reverse[@]}"; do
                dwi_nom=$(echo "${dwi_reverse[i]##*/}" | awk -F ".nii" '{print $1}')
                bids_dwi_str=$(echo "${dwi_reverse[i]}" | awk -F . '{print $1}')
                b0_nom="${tmp}/$(echo "${dwi_reverse[i]##*/}" | awk -F ".nii" '{print $1}')_b0.nii.gz"
                b0_acq=$(echo "${dwi_nom/${idBIDS}_/}"); b0_acq=$(echo "${b0_acq/_dwi/}")
                b0mat_str="${dir_warp}/${idBIDS}_from-${b0_acq}_to-${b0_refacq}${dwi_str_}_mode-image_desc-rigid_"
                b0mat="${b0mat_str}0GenericAffine.mat"

                Info "DWI rpe - Registering ${b0_acq} to ${b0_refacq}"
                Do_cmd antsRegistrationSyN.sh -d 3 -m "$b0_nom" -f "$b0_ref" -o "$b0mat_str" -t r -n "$threads" -p d
                mrconvert "${tmp}/${dwi_nom}.mif" "${tmp}/${dwi_nom}.nii.gz"
                Do_cmd antsApplyTransforms -d 3 -e 3 -i "${tmp}/${dwi_nom}.nii.gz" -r "$b0_ref" -t "$b0mat" -o "${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz" -v -u int
                Do_cmd mrconvert "${tmp}/${dwi_nom}_in-${b0_refacq}.nii.gz" -json_import "${bids_dwi_str}.json" -fslgrad "${bids_dwi_str}.bvec" "${bids_dwi_str}.bval" "${tmp}/${dwi_nom}_rpe_Ralign.mif" -force -quiet
            done
          fi

          Info "DWI-rpe: concatenatenating shells"
          dwi_0=$(echo "${dwi_reverse[0]##*/}" | awk -F ".nii" '{print $1}')
          rpe_dns_tmp="${tmp}/rpe-MP-PCA_dwi.mif"
          if [ "${#dwi_reverse[@]}" -eq 1 ]; then
            cp "${tmp}/${dwi_0}.mif" "$rpe_dns_tmp"
          else
            Do_cmd mrcat "${tmp}/*_rpe_Ralign.mif" "$rpe_cat" -nthreads "$threads"
            Do_cmd dwidenoise "$rpe_cat" "$rpe_dns_tmp" -nthreads "$threads"
          fi

          # Denoise DWI and calculate residuals
          Info "DWI-rpe: MP-PCA denoising and Gibbs ringing correction"
          Do_cmd mrdegibbs "$rpe_dns_tmp" "$rpe_dns" -nthreads "$threads"
          Do_cmd mrinfo "$rpe_dns" -json_all "${rpe_dns/mif/json}"
    else
          Info "Subject ${id} has DWI-rpe in mif, denoised and concatenaded"
    fi
else
    Info "Subject ${id} has a DWI processed, skipping denoise and concatenation"
fi

#------------------------------------------------------------------------------#
# dwifslpreproc and TOPUP preparations
if [[ ! -f "$dwi_corr" ]]; then
      Info "DWI dwifslpreproc"
      # Get parameters
      ReadoutTime=$(mrinfo "$dwi_dns" -property TotalReadoutTime)
      pe_dir=$(mrinfo "$dwi_dns" -property PhaseEncodingDirection)
      shells=($(mrinfo "$dwi_dns" -shell_bvalues))
      # Exclude shells with a threshold b-value lower than 15
      for i in "${!shells[@]}"; do if [ "${shells[i]%.*}" -le 15 ]; then unset 'shells[i]'; fi; done

      # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
      dwi_4proc=${tmp}/dwi_dns_even.mif
      dim=$(mrinfo "$dwi_dns" -size)
      dimNew=($(echo "$dim" | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i-1}}'))
      Do_cmd mrconvert "$dwi_dns" "$dwi_4proc" -coord 0 0:"${dimNew[0]}" -coord 1 0:"${dimNew[1]}" -coord 2 0:"${dimNew[2]}" -coord 3 0:end -force

      # Get the mean b-zero (un-corrected)
      dwiextract -nthreads "$threads" "$dwi_4proc" - -bzero | mrmath - mean "$tmp"/b0_meanMainPhase.mif -axis 3
      # Mean rpe QC image
      Do_cmd mrconvert "${tmp}/b0_meanMainPhase.mif" "${tmp}/b0_meanMainPhase.nii.gz"
      Do_cmd ${MICAPIPE}/functions/nifti_capture.py -img "${tmp}/b0_meanMainPhase.nii.gz" -out "${dir_QC_png}/${idBIDS}_space-dwi_pe${dwi_str_}.png"

      # Processing the reverse encoding b0
      if [[ -f "$rpe_dns" ]]; then
            rpe_dim=$(mrinfo "$rpe_dns" -ndim)
            dwgrad=$(mrinfo "$rpe_dns" -dwgrad | wc -l)

            # Mean reverse phase b0
            Info "Extracting the rpe b0(s)"
            Note "    rpe dwgrad :" "${dwgrad}"
            Note "    rpe ndim   :" "${rpe_dim}"
            if [[ "${dwgrad}" -eq 0 ]]; then
                Warning "No bval or bvecs were found the script will assumme that all the volumes are b0s!!!!"
                if [[ "$rpe_dim" -eq 3 ]]; then
                    Do_cmd mrconvert "$rpe_dns" "${tmp}/b0_meanReversePhase.nii.gz"
                elif [[ "$rpe_dim" -gt 3 ]]; then
                    mrmath "$rpe_dns" mean "${tmp}/b0_meanReversePhase.nii.gz" -axis 3 -nthreads "$threads"
                fi
            else
                if [[ "$rpe_dim" -eq 3 ]]; then
                    Do_cmd mrconvert "$rpe_dns" "${tmp}/b0_meanReversePhase.nii.gz"
                elif [[ "$rpe_dim" -gt 3 ]]; then
                    dwiextract "$rpe_dns" - -bzero | mrmath - mean "${tmp}/b0_meanReversePhase.nii.gz" -axis 3 -nthreads "$threads"
                fi
            fi
            Do_cmd mrconvert "${tmp}/b0_meanReversePhase.nii.gz" "${tmp}/b0_ReversePhase.nii.gz" -coord 0 0:"${dimNew[0]}" -coord 1 0:"${dimNew[1]}" -coord 2 0:"${dimNew[2]}" -force

            Do_cmd ${MICAPIPE}/functions/nifti_capture.py -img "${tmp}/b0_ReversePhase.nii.gz" -out "${dir_QC_png}/${idBIDS}_space-dwi_rpe${dwi_str_}.png"

            if [[ "$rpe_all" == TRUE ]]; then
                # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc-TOPUP).
                drpe_4proc=${tmp}/rpe_dns_even.mif
                dwi_all=${tmp}/all_dwi_dns_even.mif
                mrconvert "${rpe_dns}" "$drpe_4proc" -coord 0 0:"${dimNew[0]}" -coord 1 0:"${dimNew[1]}" -coord 2 0:"${dimNew[2]}" -coord 3 0:end -force
                Do_cmd mrcat "${dwi_4proc}" "${drpe_4proc}" "$dwi_all"
                dwi_4proc="$dwi_all"
                opt="-rpe_all"

            elif [[ "$rpe_all" == FALSE ]]; then
                b0_pair="${tmp}/b0_pair.mif"
                # Concatenate the pe and rpe b0s
                Do_cmd mrcat "${tmp}/b0_meanMainPhase.nii.gz" "${tmp}/b0_ReversePhase.nii.gz" "$b0_pair" -nthreads "$threads"
                opt="-rpe_pair -align_seepi -se_epi ${b0_pair}"
            fi
      else
            Warning "Reverse phase encoding image was not found it will be omitted"
            opt='-rpe_none'
      fi

      Info "dwifslpreproc parameters:"
      Note "Shell values        :" "${shells[*]}"
      Note "DWI main dimensions :" "$(mrinfo "$dwi_dns" -size)"
      if [ -f "${dwi_reverse[0]}" ]; then Note "DWI rpe dimensions  :" "$(mrinfo "$rpe_dns" -size)"; fi
      Note "DWI to process dim  :" "$(mrinfo "$dwi_4proc" -size)"
      Note "pe_dir              :" "$pe_dir"
      Note "Readout Time        :" "$ReadoutTime"
      Note "Options             :" "$opt"

      # Preprocess each shell
      # DWIs all acquired with a single fixed phase encoding; but additionally a
      # pair of b=0 images with reversed phase encoding to estimate the inhomogeneity field:
      echo -e "COMMAND --> dwifslpreproc $dwi_4proc $dwi_corr $opt -pe_dir $pe_dir -readout_time $ReadoutTime -eddy_options \" --data_is_shelled --slm=linear --repol\" -nthreads $threads -nocleanup -scratch $tmp -force"
      dwifslpreproc "$dwi_4proc" "$dwi_corr" $opt -pe_dir "$pe_dir" -readout_time "$ReadoutTime" -eddy_options " --data_is_shelled --slm=linear --repol" -nthreads "$threads" -nocleanup -scratch "$tmp" -force
      # Step QC
      if [[ ! -f "$dwi_corr" ]]; then Error "dwifslpreproc failed, check the logs"; exit;
      else
          # Bias field correction DWI
          dwi_n4="${tmp}/${idBIDS}_space-dwi_desc-dwi_preproc_N4.mif"
          Do_cmd dwibiascorrect ants "$dwi_corr" "$dwi_n4" -force -nthreads "$threads" -scratch "$tmp"
          Do_cmd mv "$dwi_n4" "$dwi_corr"
          Do_cmd mrinfo "$dwi_corr" -json_all "${dwi_corr/mif/json}"
          json_dwipreproc $dwi_dns $dwi_corr ${proc_dwi}/${idBIDS}_desc-preproc_dwi.json ${rpe_dns}
          Do_cmd rm "$dwi_dns"
          # eddy_quad Quality Check
          Do_cmd cd "$tmp"/dwifslpreproc*
          Do_cmd eddy_quad dwi_post_eddy -idx eddy_indices.txt -par eddy_config.txt -m eddy_mask.nii -b bvals -o "${dir_QC}/eddy_QC${dwi_str_}"
          Do_cmd cd "$tmp"

          # Copy eddy parameters
          eddy_DIR=${proc_dwi}/eddy
          if [ ! -d "$eddy_DIR" ]; then Do_cmd mkdir "$eddy_DIR"; fi
          Do_cmd cp -rf "${tmp}"/dwifslpreproc*/*eddy.eddy* "$eddy_DIR"
          Do_cmd chmod 770 -R "$eddy_DIR"/*
          ((Nsteps++))
      fi
else
      Info "Subject ${id} has a DWI processed"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Registration of corrected DWI-b0 to T1nativepro

dwi_mask="${proc_dwi}/${idBIDS}_space-dwi_desc-brain_mask.nii.gz"
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0.nii.gz" # This should be a NIFTI for compatibility with ANTS
str_dwi_affine="${dir_warp}/${idBIDS}_space-dwi_from-dwi${dwi_str_}_to-nativepro_mode-image_desc-affine_"
mat_dwi_affine="${str_dwi_affine}0GenericAffine.mat"
T1nativepro_in_dwi="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro.nii.gz"

if [[ ! -f "$T1nativepro_in_dwi" ]] || [[ ! -f "$mat_dwi_affine" ]]; then
      Info "Affine registration from DWI-b0 to T1nativepro"
      # Corrected DWI-b0s mean for registration
      dwiextract -force -nthreads "$threads" "$dwi_corr" - -bzero | mrmath - mean "$dwi_b0" -axis 3 -force

      # Register DWI-b0 mean corrected to T1nativepro
      Do_cmd antsRegistrationSyN.sh -d 3 -f "$T1nativepro_brain" -m "$dwi_b0" -o "$str_dwi_affine" -t a -n "$threads" -p d
      # Apply inverse transformation T1nativepro to DWI-b0 space
      Do_cmd antsApplyTransforms -d 3 -i "$T1nativepro" -r "$dwi_b0" -t ["$mat_dwi_affine",1] -o "$T1nativepro_in_dwi" -v -u int
      if [[ -f "$T1nativepro_in_dwi" ]]; then ((Nsteps++)); fi

      #------------------------------------------------------------------------------#
      Info "Creating DWI binary mask of processed volumes"
      # Create a binary mask of the DWI
      Do_cmd antsApplyTransforms -d 3 -i "$MNI152_mask" \
              -r "$dwi_b0" \
              -n GenericLabel -t ["$mat_dwi_affine",1] -t ["$T1_MNI152_affine",1] -t "$T1_MNI152_InvWarp" \
              -o "${tmp}"/dwi_mask.nii.gz -v
      Do_cmd maskfilter "${tmp}"/dwi_mask.nii.gz erode -npass 1 "$dwi_mask"
      # Do_cmd dwi2mask "$dwi_corr" "$dwi_mask" -nthreads "$threads"
      if [[ -f "$dwi_mask" ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has an affine transformation from T1w to DWI-b0 space"; Nsteps=$((Nsteps + 2))
fi

#------------------------------------------------------------------------------#
# Get some basic metrics.
dwi_dti="${proc_dwi}/${idBIDS}_space-dwi_model-DTI.mif"
dti_FA="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.mif"
dti_AD="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-AD.mif"
dti_RD="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-RD.mif"
dti_ADC="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-ADC.mif"
if [[ ! -f "$dti_FA" ]]; then
      Info "Calculating basic DTI metrics"
      dwi2tensor -mask "$dwi_mask" -nthreads "$threads" "$dwi_corr" "$dwi_dti"
      tensor2metric -nthreads "$threads" -fa "$dti_FA" -adc "$dti_ADC" -ad "$dti_AD" -rd "$dti_RD" "$dwi_dti"
      Do_cmd mrinfo "$dwi_dti" -json_all "${dwi_dti/mif/json}"
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
            rf=dhollander
            # Response function
            rf_wm="${tmp}/${id}_response_wm_${rf}.txt"
            rf_gm="${tmp}/${id}_response_gm_${rf}.txt"
            rf_csf="${tmp}/${id}_response_csf_${rf}.txt"
            # Fiber Orientation Distriution
            fod_wm="${tmp}/${id}_wm_fod.mif"
            fod_gm="${tmp}/${id}_gm_fod.mif"
            fod_csf="${tmp}/${id}_csf_fod.mif"

            Do_cmd dwi2response "$rf" -nthreads "$threads" "$dwi_corr" "$rf_wm" "$rf_gm" "$rf_csf" -mask "$dwi_mask"
            Do_cmd dwi2fod -nthreads "$threads" msmt_csd "$dwi_corr" \
                "$rf_wm" "$fod_wm" \
                "$rf_gm" "$fod_gm" \
                "$rf_csf" "$fod_csf" \
                -mask "$dwi_mask"
      if [ "${#shells[@]}" -ge 2 ]; then
            Do_cmd mtnormalise "$fod_wm" "$fod_wmN" "$fod_gm" "$fod_gmN" "$fod_csf" "$fod_csfN" -nthreads "$threads" -mask "$dwi_mask"
            Do_cmd mrinfo "$fod_wmN" -json_all "${fod_wmN/mif/json}"
            Do_cmd mrinfo "$fod_gmN" -json_all "${fod_gmN/mif/json}"
            Do_cmd mrinfo "$fod_csfN" -json_all "${fod_csfN/mif/json}"
      else
            Do_cmd mtnormalise -nthreads "$threads" -mask "$dwi_mask" "$fod_wm" "$fod_wmN"
            Do_cmd mrinfo "$fod_wmN" -json_all "${fod_wmN/mif/json}"
      fi
      if [[ -f "$fod_wmN" ]]; then ((Nsteps++)); fi
else
      Info "Subject ${id} has Fiber Orientation Distribution file"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# Non-linear registration between DWI space and T1w
dwi_SyN_str="${dir_warp}/${idBIDS}_space-dwi_from-dwi${dwi_str_}_to-dwi_mode-image_desc-SyN_"
dwi_SyN_warp="${dwi_SyN_str}1Warp.nii.gz"
dwi_SyN_Invwarp="${dwi_SyN_str}1InverseWarp.nii.gz"
dwi_SyN_affine="${dwi_SyN_str}0GenericAffine.mat"
dwi_5tt="${proc_dwi}/${idBIDS}_space-dwi_desc-5tt.nii.gz"


if [[ ! -f "$dwi_SyN_warp" ]] || [[ ! -f "$dwi_5tt" ]]; then
    dwi_in_T1nativepro="${proc_struct}/${idBIDS}_space-nativepro_desc-dwi.nii.gz" # Only for QC
    T1nativepro_in_dwi_brain="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro-brain.nii.gz"
    fod="${tmp}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz"
    Do_cmd fslmaths "$T1nativepro_in_dwi" -mul "$dwi_mask" "$T1nativepro_in_dwi_brain"
    Do_cmd mrconvert -coord 3 0 "$fod_wmN" "$fod"

    if [[ ${regAffine}  == "FALSE" ]]; then
        Info "Non-linear registration from T1w_dwi-space to DWI"
        T1nativepro_in_dwi_NL="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro_SyN.nii.gz"
        Do_cmd antsRegistrationSyN.sh -d 3 -m "$T1nativepro_in_dwi_brain" -f "$fod" -o "$dwi_SyN_str" -t s -n "$threads"
        export reg="Affine+SyN"
        trans_T12dwi="-t ${dwi_SyN_warp} -t ${dwi_SyN_affine} -t [${mat_dwi_affine},1]" # T1nativepro to DWI
        trans_dwi2T1="-t ${mat_dwi_affine} -t [${dwi_SyN_affine},1] -t ${dwi_SyN_Invwarp}"  # DWI to T1nativepro
        if [[ -f "$dwi_SyN_warp" ]]; then ((Nsteps++)); fi

    elif [[ ${regAffine}  == "TRUE" ]]; then
        Info "Only affine registration from T1w_dwi-space to DWI"; ((Nsteps++))
        T1nativepro_in_dwi_NL="${proc_dwi}/${idBIDS}_space-dwi_desc-t1w_nativepro_Affine.nii.gz"
        trans_T12dwi="-t [${mat_dwi_affine},1]"
        trans_dwi2T1="-t ${mat_dwi_affine}"
    fi

    Info "Registering T1w-nativepro and 5TT to DWI-b0 space, and DWI-b0 to T1w-nativepro"
    # Apply transformation DWI-b0 space to T1nativepro
    Do_cmd antsApplyTransforms -d 3 -r "$T1nativepro_brain" -i "$dwi_b0" -r "$T1nativepro_brain" "$trans_dwi2T1" -o "$dwi_in_T1nativepro" -v -u int
    # Apply transformation of each DTI derived map to T1nativepro
    for metric in FA AD RD ADC; do
        dti_map="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-${metric}.mif"
        dti_map_nativepro="${proc_dwi}/${idBIDS}_space-nativepro_model-DTI_map-${metric}.nii.gz"
        Do_cmd mrconvert "$dti_map" "${dti_map/.mif/.nii.gz}"
        Do_cmd antsApplyTransforms -d 3 -r "$T1nativepro_brain" -i "${dti_map/.mif/.nii.gz}" -r "$T1nativepro_brain" "$trans_dwi2T1" -o "$dti_map_nativepro" -v -u int -n NearestNeighbor
    done
    # Apply transformation T1nativepro to DWI space
    Do_cmd antsApplyTransforms -d 3 -r "$fod" -i "$T1nativepro" "$trans_T12dwi" -o "$T1nativepro_in_dwi_NL" -v -u int
    # Apply transformation 5TT to DWI space
    Do_cmd antsApplyTransforms -d 3 -r "$fod" -i "$T15ttgen" "$trans_T12dwi" -o "$dwi_5tt" -v -e 3 -n linear
    if [[ -f "$dwi_5tt" ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a registration from T1w_dwi-space to DWI"; Nsteps=$((Nsteps + 2))
fi
# Remove unused warped files
Do_cmd rm -rf ${dir_warp}/*Warped.nii.gz 2>/dev/null

#------------------------------------------------------------------------------#
# Gray matter White matter interface mask
dwi_gmwmi="${proc_dwi}/${idBIDS}_space-dwi_desc-gmwmi-mask.mif"
if [[ ! -f "$dwi_gmwmi" ]]; then
      Info "Calculating Gray matter White matter interface mask"
      Do_cmd 5tt2gmwmi "$dwi_5tt" "$dwi_gmwmi"
      Do_cmd mrinfo "$dwi_gmwmi" -json_all "${dwi_gmwmi/mif/json}"; ((Nsteps++))
else
      Info "Subject ${id} has Gray matter White matter interface mask"; ((Nsteps++))
fi

#------------------------------------------------------------------------------#
# QC of the tractography
tracts=1M
tdi_1M="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD1-${tracts}_tdi.mif"
export tckjson="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD1-${tracts}_tractography.json"
if [[ ! -f "$tdi_1M" ]]; then
  Info "Creating a track density image for quality check"
  tck_1M="${tmp}/${idBIDS}_space-dwi_desc-iFOD1-${tracts}_tractography.tck"
  Do_cmd tckgen -nthreads "$threads" \
      "$fod_wmN" \
      "$tck_1M" \
      -act "$dwi_5tt" \
      -crop_at_gmwmi \
      -backtrack \
      -seed_gmwmi "$dwi_gmwmi" \
      -algorithm iFOD1 \
      -step 0.5 \
      -angle 22.5 \
      -cutoff 0.06 \
      -maxlength 400 \
      -minlength 10 \
      -power 1.0 \
      -select "$tracts"

  Do_cmd tckmap -vox 1,1,1 -dec -nthreads "$threads" "$tck_1M" "$tdi_1M"
  if [[ -f "$tdi_1M" ]]; then ((Nsteps++)); fi
  tck_json iFOD1 0.5 22.5 0.06 400 10 seed_gmwmi "$tck_1M"
else
      Info "Subject ${id} has a Tract Density Image for QC 1M streamlines"; ((Nsteps++))
fi

# -----------------------------------------------------------------------------------------------
# QC: Input files
QC_proc-dwi "${dir_QC}/micapipe_QC_proc-dwi${dwi_str_}.txt"

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
N=10
if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "DWI processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\t\tSteps completed: $(printf "%02d" "$Nsteps")/10
\tStatus          : ${status}
\tCheck logs:
$(ls "${dir_logs}"/proc_dwi_*.txt)"
# Print QC stamp
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_dwi${dwi_str_}" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "proc_dwi${dwi_str_}" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
cleanup "$tmp" "$nocleanup" "$here"
