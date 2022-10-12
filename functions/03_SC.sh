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
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
tracts=$8
autoTract=$9
keep_tck=${10}
dwi_str=${11}
filter=SIFT2
PROC=${12}
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE/functions/utilities.sh"

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Update path for multiple acquisitions processing
if [[ "${dwi_str}" != "DEFAULT" ]]; then
  dwi_str="acq-${dwi_str/acq-/}"
  dwi_str_="_${dwi_str}"
  export proc_dwi=$subject_dir/dwi/"${dwi_str}"
  export dwi_cnntm=$proc_dwi/connectomes
  export autoTract_dir=$proc_dwi/auto_tract
else
  dwi_str=""; dwi_str_=""
fi

# Check inputs: DWI post TRACTOGRAPHY
fod_wmN="${proc_dwi}/${idBIDS}_space-dwi_model-CSD_map-FOD_desc-wmNorm.mif"
dwi_5tt="${proc_dwi}/${idBIDS}_space-dwi_desc-5tt.nii.gz"
dwi_b0="${proc_dwi}/${idBIDS}_space-dwi_desc-b0.nii.gz"
dwi_mask="${proc_dwi}/${idBIDS}_space-dwi_desc-brain_mask.nii.gz"
str_dwi_affine="${dir_warp}/${idBIDS}_space-dwi_from-dwi${dwi_str_}_to-nativepro_mode-image_desc-affine_"
mat_dwi_affine="${str_dwi_affine}0GenericAffine.mat"
dwi_SyN_str="${dir_warp}/${idBIDS}_space-dwi_from-dwi${dwi_str_}_to-dwi_mode-image_desc-SyN_"
dwi_SyN_warp="${dwi_SyN_str}1Warp.nii.gz"
dwi_SyN_affine="${dwi_SyN_str}0GenericAffine.mat"
dti_FA="${proc_dwi}/${idBIDS}_space-dwi_model-DTI_map-FA.mif"
lut_sc="${util_lut}/lut_subcortical-cerebellum_mics.csv"
# from proc_structural
T1str_nat="${idBIDS}_space-nativepro_t1w_atlas"
T1_seg_cerebellum="${dir_volum}/${T1str_nat}-cerebellum.nii.gz"
T1_seg_subcortex="${dir_volum}/${T1str_nat}-subcortical.nii.gz"
# TDI output
tdi="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tdi.mif"

# Check inputs
if [ ! -f "$fod_wmN"  ]; then Error "Subject $id doesn't have FOD:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f "$dwi_b0" ]; then Error "Subject $id doesn't have dwi_b0:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f "$mat_dwi_affine" ]; then Error "Subject $id doesn't have an affine mat from T1nativepro to DWI space:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f "$dwi_5tt" ]; then Error "Subject $id doesn't have 5tt in dwi space:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f "$T1_seg_cerebellum" ]; then Error "Subject $id doesn't have cerebellar segmentation:\n\t\tRUN -post_structural"; exit; fi
if [ ! -f "$T1_seg_subcortex" ]; then Error "Subject $id doesn't have subcortical segmentation:\n\t\tRUN -post_structural"; exit; fi
if [ ! -f "$dwi_mask" ]; then Error "Subject $id doesn't have DWI binary mask:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f "$dti_FA" ]; then Error "Subject $id doesn't have a FA:\n\t\tRUN -proc_dwi"; exit; fi
if [ ! -f "$dwi_SyN_affine" ]; then Warning "Subject $id doesn't have an SyN registration, only AFFINE will be apply"; regAffine="TRUE"; else regAffine="FALSE"; fi

# -----------------------------------------------------------------------------------------------
# Check IF output exits and WARNING
N=$(ls "${dwi_cnntm}"/"${idBIDS}"_space-dwi_atlas-*_desc-iFOD2-"${tracts}"-"${filter}"_full-connectome.txt 2>/dev/null | wc -l)
if [ "$N" -gt 3 ]; then Warning " Connectomes with $tracts streamlines already exist!!
  If you want to re-run the $tracts tractogram or add parcellations first clean the outpus:
    micapipe_cleanup -SC -sub $id -out ${out/"/micapipe"/} -bids $BIDS -tracts ${tracts}"; fi

#------------------------------------------------------------------------------#
Title "Tractography and structural connectomes\n\t\tmicapipe $Version, $PROC"
micapipe_software
Note "Number of streamlines:" "${tracts}"
Note "Auto-tractograms     :" "${autoTract}"
Note "Saving tractography  :" "${keep_tck}"
Note "Saving temporal dir  :" "${nocleanup}"
Note "MRtrix will use      :" "${threads} threads"
Note "DWi acquisition      :" "${dwi_str}"

#	Timer
aloita=$(date +%s)
Nsteps=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_post-dwi_${id}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Create Connectomes directory for the outpust
[[ ! -d "$dwi_cnntm" ]] && Do_cmd mkdir -p "$dwi_cnntm" && chmod -R 770 "$dwi_cnntm"
Do_cmd cd "$tmp"

# -----------------------------------------------------------------------------------------------
# Prepare the segmentatons
parcellations=($(find "${dir_volum}" -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"))
dwi_cere="${proc_dwi}/${idBIDS}_space-dwi_atlas-cerebellum.nii.gz"
dwi_subc="${proc_dwi}/${idBIDS}_space-dwi_atlas-subcortical.nii.gz"

# Transformations from T1nativepro to DWI
if [[ ${regAffine}  == "FALSE" ]]; then
    trans_T12dwi="-t ${dwi_SyN_warp} -t ${dwi_SyN_affine} -t [${mat_dwi_affine},1]"
elif [[ ${regAffine}  == "TRUE" ]]; then
    trans_T12dwi="-t [${mat_dwi_affine},1]"
fi

if [[ ! -f "$dwi_cere" ]]; then Info "Registering Cerebellar parcellation to DWI-b0 space"
      Do_cmd antsApplyTransforms -d 3 -e 3 -i "$T1_seg_cerebellum" -r "$dwi_b0" -n GenericLabel "$trans_T12dwi" -o "$dwi_cere" -v -u int
      if [[ -f "$dwi_cere" ]]; then ((Nsteps++)); fi
      # Threshold cerebellar nuclei (29,30,31,32,33,34) and add 100
      # Do_cmd fslmaths $dwi_cere -uthr 28 $dwi_cere
      Do_cmd fslmaths "$dwi_cere" -bin -mul 100 -add "$dwi_cere" "$dwi_cere"
else Info "Subject ${id} has a Cerebellar segmentation in DWI space"; ((Nsteps++)); fi

if [[ ! -f "$dwi_subc" ]]; then Info "Registering Subcortical parcellation to DWI-b0 space"
    Do_cmd antsApplyTransforms -d 3 -e 3 -i "$T1_seg_subcortex" -r "$dwi_b0" -n GenericLabel "$trans_T12dwi" -o "$dwi_subc" -v -u int
    # Remove brain-stem (label 16)
    Do_cmd fslmaths "$dwi_subc" -thr 16 -uthr 16 -binv -mul "$dwi_subc" "$dwi_subc"
    if [[ -f "$dwi_subc" ]]; then ((Nsteps++)); fi
else Info "Subject ${id} has a Subcortical segmentation in DWI space"; ((Nsteps++)); fi

# -----------------------------------------------------------------------------------------------
# Generate probabilistic tracts
tck="${tmp}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tractography.tck"
if [ ! -f "$tdi" ]; then
    Info "Building the ${tracts} streamlines connectome!!!"
    export tckjson="${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tractography.json"
    weights=${tmp}/SIFT2_${tracts}.txt
    Do_cmd tckgen -nthreads "$threads" \
        "$fod_wmN" \
        "$tck" \
        -act "$dwi_5tt" \
        -crop_at_gmwmi \
        -backtrack \
        -seed_dynamic "$fod_wmN" \
        -algorithm iFOD2 \
        -step 0.5 \
        -angle 22.5 \
        -cutoff 0.06 \
        -maxlength 400 \
        -minlength 10 \
        -select "$tracts"

    # Exit if tractography fails
    if [ ! -f "$tck" ]; then Error "Tractogram failed, check the logs: $(ls -Art "$dir_logs"/post-dwi_*.txt | tail -1)"; exit; fi

    # json file of tractogram
    tck_json iFOD2 0.5 22.5 0.06 400 10 seed_dynamic "$tck"

    # SIFT2
    Do_cmd tcksift2 -nthreads "$threads" "$tck" "$fod_wmN" "$weights"

    # TDI for QC
    Info "Creating a Track Density Image (tdi) of the $tracts connectome for QC"
    Do_cmd tckmap -vox 1,1,1 -dec -nthreads "$threads" "$tck" "$tdi" -force
    ((Nsteps++))
else
    Warning "SC has been processed for Subject $id: TDI of ${tracts} was found"; ((Nsteps++))
fi


# -----------------------------------------------------------------------------------------------
# Build the Connectomes
for seg in "${parcellations[@]}"; do
    parc_name=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    connectome_str="${dwi_cnntm}/${idBIDS}_space-dwi_atlas-${parc_name}_desc-iFOD2-${tracts}-${filter}"
    lut="${util_lut}/lut_${parc_name}_mics.csv"
    dwi_cortex="${tmp}/${id}_${parc_name}-cor_dwi.nii.gz" # Segmentation in dwi space

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical connectomes
    if [[ ! -f "${connectome_str}_cor-connectome.txt" ]]; then
        Info "Building $parc_name cortical connectome"
        # Take parcellation into DWI space
        Do_cmd antsApplyTransforms -d 3 -e 3 -i "$seg" -r "$dwi_b0" -n GenericLabel "$trans_T12dwi" -o "$dwi_cortex" -v -u int
        # Remove the medial wall
        for i in 1000 2000; do Do_cmd fslmaths "$dwi_cortex" -thr "$i" -uthr "$i" -binv -mul "$dwi_cortex" "$dwi_cortex"; done

        # Build the Cortical connectomes
        Do_cmd tck2connectome -nthreads "$threads" \
            "$tck" "$dwi_cortex" "${connectome_str}_cor-connectome.txt" \
            -tck_weights_in "$weights" -quiet
        Do_cmd Rscript "$MICAPIPE"/functions/connectome_slicer.R --conn="${connectome_str}_cor-connectome.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"

        # Calculate the edge lenghts
        Do_cmd tck2connectome -nthreads "$threads" \
            "$tck" "$dwi_cortex" "${connectome_str}_cor-edgeLengths.txt" \
            -tck_weights_in "$weights" -scale_length -stat_edge mean -quiet
        Do_cmd Rscript "$MICAPIPE"/functions/connectome_slicer.R --conn="${connectome_str}_cor-edgeLengths.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"
        if [[ -f "${connectome_str}_cor-connectome.txt" ]]; then ((Nsteps++)); fi
    else
        ((Nsteps++))
    fi

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical connectomes (-sub)
    if [[ ! -f "${connectome_str}_sub-connectome.txt" ]]; then
        Info "Building $parc_name cortical-subcortical connectome"
        dwi_cortexSub="${tmp}/${id}_${parc_name}-sub_dwi.nii.gz"
        Do_cmd fslmaths "$dwi_cortex" -binv -mul "$dwi_subc" -add "$dwi_cortex" "$dwi_cortexSub" -odt int # added the subcortical parcellation

        # Build the Cortical-Subcortical connectomes
        Do_cmd tck2connectome -nthreads "$threads" \
            "$tck" "$dwi_cortexSub" "${connectome_str}_sub-connectome.txt" \
            -tck_weights_in "$weights" -quiet
        Do_cmd Rscript "$MICAPIPE"/functions/connectome_slicer.R --conn="${connectome_str}_sub-connectome.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"

        # Calculate the edge lenghts
        Do_cmd tck2connectome -nthreads "$threads" \
            "$tck" "$dwi_cortexSub" "${connectome_str}_sub-edgeLengths.txt" \
            -tck_weights_in "$weights" -scale_length -stat_edge mean -quiet
        Do_cmd Rscript "$MICAPIPE"/functions/connectome_slicer.R --conn="${connectome_str}_sub-edgeLengths.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"
        if [[ -f "${connectome_str}_sub-connectome.txt" ]]; then ((Nsteps++)); fi
    else
        ((Nsteps++))
    fi

    # -----------------------------------------------------------------------------------------------
    # Build the Cortical-Subcortical-Cerebellar connectomes (-sub-cereb)
    if [[ ! -f "${connectome_str}_full-connectome.txt" ]]; then
        Info "Building $parc_name cortical-subcortical-cerebellum connectome"
        dwi_all="${tmp}/${id}_${parc_name}-full_dwi.nii.gz"
        Do_cmd fslmaths "$dwi_cortex" -binv -mul "$dwi_cere" -add "$dwi_cortexSub" "$dwi_all" -odt int # added the cerebellar parcellation

        # Build the Cortical-Subcortical-Cerebellum connectomes
        Do_cmd tck2connectome -nthreads "$threads" \
            "$tck" "$dwi_all" "${connectome_str}_full-connectome.txt" \
            -tck_weights_in "$weights" -quiet
        Do_cmd Rscript "$MICAPIPE"/functions/connectome_slicer.R --conn="${connectome_str}_full-connectome.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"

        # Calculate the edge lenghts
        Do_cmd tck2connectome -nthreads "$threads" \
            "$tck" "$dwi_all" "${connectome_str}_full-edgeLengths.txt" \
            -tck_weights_in "$weights" -scale_length -stat_edge mean -quiet
        Do_cmd Rscript "$MICAPIPE"/functions/connectome_slicer.R --conn="${connectome_str}_full-edgeLengths.txt" --lut1="$lut_sc" --lut2="$lut" --mica="$MICAPIPE"
        if [[ -f "${connectome_str}_full-connectome.txt" ]]; then ((Nsteps++)); fi
    else
        ((Nsteps++))
    fi
done

# Change connectome permissions
chmod 770 -R "$dwi_cnntm"/* 2>/dev/null

# -----------------------------------------------------------------------------------------------
# Compute Auto-Tractography
if [ "$autoTract" == "TRUE" ]; then
    Info "Running Auto-tract"
    autoTract_dir="$proc_dwi"/auto_tract
    [[ ! -d "$autoTract_dir" ]] && Do_cmd mkdir -p "$autoTract_dir"
    fa_niigz="${tmp}/${id}_dti_FA.nii.gz"
    Do_cmd mrconvert "$dti_FA" "$fa_niigz"
    echo -e "\033[38;5;118m\nCOMMAND -->  \033[38;5;122m03_auto_tracts.sh -tck $tck -outbase $autoTract_dir/${id} -mask $dwi_mask -fa $fa_niigz -tmpDir $tmp -keep_tmp  \033[0m"
    "$MICAPIPE"/functions/03_auto_tracts.sh -tck "$tck" -outbase "${autoTract_dir}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}-${filter}" -mask "$dwi_mask" -fa "$fa_niigz" -weights "$weights" -tmpDir "$tmp" -keep_tmp
fi

# -----------------------------------------------------------------------------------------------
# save the tractogram and the SIFT2 weights
if [ "$keep_tck" == "TRUE" ]; then Do_cmd mv "$tck" "$proc_dwi"; Do_cmd mv "$weights" "${proc_dwi}/${idBIDS}_space-dwi_desc-iFOD2-${tracts}_tractography_weights.txt"; fi

# -----------------------------------------------------------------------------------------------
# QC notification of completition
QC_SC "${dir_QC}/micapipe_QC_SC${dwi_str_}.txt"
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
N="$(( 3 + ${#parcellations[*]} * 3))"
if [ "$Nsteps" -eq "$N" ]; then status="COMPLETED"; else status="INCOMPLETE"; fi
Title "DWI-post TRACTOGRAPHY processing ended in \033[38;5;220m $(printf "%0.3f\n" "$eri") minutes \033[38;5;141m:
\tSteps completed : $(printf "%02d" "$Nsteps")/$(printf "%02d" "$N")
\tStatus          : ${status}
\tCheck logs      : $(ls "$dir_logs"/SC_*.txt)"
# Print QC stamp
micapipe_procStatus "${id}" "${SES/ses-/}" "SC-${tracts}${dwi_str_}" "${out}/micapipe_processed_sub.csv"
micapipe_procStatus "${id}" "${SES/ses-/}" "SC-${tracts}${dwi_str_}" "${dir_QC}/${idBIDS}_micapipe_processed.csv"
cleanup "$tmp" "$nocleanup" "$here"
