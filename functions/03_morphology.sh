#!/bin/bash
#
# Cortical morphology metrics processing:
#
# Generates vertexwise (native, fsa5, and conte69) outputs for:
#   Cortical Thickness
#   Mean Curvature
#
# This workflow makes use of freesurfer outputs and custom python scripts
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micapipe/tree/master/parcellations
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out Directory
#
BIDS=$1
id=$2
out=$3
SES=$4
nocleanup=$5
threads=$6
tmpDir=$7
PROC=$8
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ];then
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfaceProc ${post_struct_json} | awk -F '"' '{print $4}')
set_surface_directory "${recon}"

#------------------------------------------------------------------------------#
Title "Cortical morphology mapping"

# Freesurfer SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Create script specific temp directory
tmp_morph="${tmpDir}/${RANDOM}_micapipe_post-morpho_${idBIDS}"
Do_cmd mkdir -p "$tmp_morph"

# TRAP in case the script fails
trap 'cleanup $tmp_morph $nocleanup $here' SIGINT SIGTERM

# Make output directory
outDir="${dir_maps}"
[[ ! -d "$outDir" ]] && Do_cmd mkdir -p "$outDir"

# Data location
dataDir="${dir_subjsurf}/surf"

#------------------------------------------------------------------------------#
### Register morphology to surfaces ###
function reg_surfaces(){
  morph_data=$1
  smooth=$2
  Info "Mapping ${morph_data}"
  # Register to fsa5 and apply 10mm smooth
  if [[ ! -f "${outDir}/${idBIDS}_space-fsaverage5_desc-rh_${morph_data}_10mm.mgh" ]]; then
      for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
          # Convert native file to mgh and save in output directory
          Do_cmd mri_convert "${dataDir}/${hemi}.${morph_data} ${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsnative_label-${morph_data}.mgh"

          Do_cmd mri_surf2surf --hemi "$hemi" \
              --srcsubject "$idBIDS" \
              --srcsurfval "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsnative_label-${morph_data}.mgh" \
              --trgsubject fsaverage5 \
              --trgsurfval "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsaverage5_label-${morph_data}.mgh"

          Do_cmd mri_surf2surf --hemi "$hemi" \
              --fwhm-trg "${smooth}" \
              --srcsubject "$idBIDS" \
              --srcsurfval "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsnative_label-${morph_data}.mgh" \
              --trgsubject fsaverage5 \
              --trgsurfval "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsaverage5_label-${morph_data}_${smooth}mm.mgh"
              if [[ -f "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsaverage5_label-${morph_data}_${smooth}mm.mgh" ]]; then ((Nsteps++)); fi
      done
  else
      Info "Subject ${id} cortical ${morph_data} is registered to fsa5"
  fi

  # Register to conte69 and apply 10mm smooth
  if [[ ! -f "${outDir}/${idBIDS}_space-conte69-32k_desc-rh_${morph_data}_${smooth}mm.mgh" ]]; then
      for hemi in lh rh; do
          [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
          HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])

          Do_cmd mri_convert "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-fsnative_label-${morph_data}.mgh" "${tmp_morph}/${hemi}_${morph_data}.func.gii"

          Do_cmd wb_command -metric-resample \
              "${tmp_morph}/${hemi}_${morph_data}.func.gii" \
              "${dir_conte69}/${idBIDS}_hemi-${HEMICAP}_space-fsnative_label-sphere.surf.gii" \
              "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
              ADAP_BARY_AREA \
              "${tmp_morph}/${hemi}_${morph_data}_c69-32k.func.gii" \
              -area-surfs \
              "${dir_subjsurf}/surf/${hemi}.midthickness.surf.gii" \
              "${dir_conte69}/${idBIDS}_hemi-${HEMICAP}_space-conte69-32k_label-midthickness.surf.gii"

          Do_cmd mri_convert "${tmp_morph}/${hemi}_${morph_data}_c69-32k.func.gii" "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-conte69-32k_label-${morph_data}.mgh"

          # Smoothing
          Do_cmd wb_command -metric-smoothing \
              "${util_surface}/fsaverage.${HEMICAP}.midthickness_orig.32k_fs_LR.surf.gii" \
              "${tmp_morph}/${hemi}_${morph_data}_c69-32k.func.gii" \
              10 \
              "${tmp_morph}/${hemi}_${morph_data}_${smooth}mm_c69-32k.func.gii"

          Do_cmd mri_convert "${tmp_morph}/${hemi}_${morph_data}_${smooth}mm_c69-32k.func.gii" "${outDir}/${idBIDS}_hemi-${HEMICAP}_space-conte69-32k_label-${morph_data}_${smooth}mm.mgh"
      done
  else
      Info "Subject ${idBIDS} cortical ${morph_data} is registered to conte69"
  fi
}

reg_surfaces "thickness" 10
reg_surfaces "curv" 10

#------------------------------------------------------------------------------#
# End if module has been processed
cleanup "$tmp_morph" "$nocleanup" "$here"
