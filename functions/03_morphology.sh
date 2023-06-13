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
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source "$MICAPIPE"/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfRecon "${post_struct_json}" | awk -F '"' '{print $4}')
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
  Info "Mapping ${morph_data}"
  # Register to fsa5 and apply
  if [[ ! -f "${outDir}/${idBIDS}_hemi-R_surf-fsaverage5_label-${morph_data}.func.gii" ]]; then
      for hemi in lh rh; do
        [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=$(echo "$hemisphere" | tr [:lower:] [:upper:])
        surf_id="${idBIDS}_hemi-${HEMICAP}_surf"
          # Convert native file to mgh and save in output directory
          Do_cmd mri_convert "${dataDir}/${hemi}.${morph_data} ${tmp_morph}/${surf_id}-fsnative_label-${morph_data}.mgh"
          Do_cmd mri_convert "${dataDir}/${hemi}.${morph_data} ${outDir}/${surf_id}-fsnative_label-${morph_data}.func.gii"

          Do_cmd mri_surf2surf --hemi "$hemi" \
              --srcsubject "$idBIDS" \
              --srcsurfval "${tmp_morph}/${surf_id}-fsnative_label-${morph_data}.mgh" \
              --trgsubject fsaverage5 \
              --trgsurfval "${tmp_morph}/${surf_id}-fsaverage5_label-${morph_data}.mgh"
          Do_cmd mri_convert "${tmp_morph}/${surf_id}-fsaverage5_label-${morph_data}.mgh" "${outDir}/${surf_id}-fsaverage5_label-${morph_data}.func.gii"

          if [[ -f "${outDir}/${surf_id}-fsaverage5_label-${morph_data}.func.gii" ]]; then ((Nsteps++)); fi
      done
  else
      Info "Subject ${id} cortical ${morph_data} is registered to fsa5"
  fi

  # Register to fsLR-32k and fsLR-5k
  if [[ ! -f "${outDir}/${idBIDS}_hemi-R_surf-fsLR-5k_label-${morph_data}.func.gii" ]]; then
    for Surf in "fsLR-32k" "fsLR-5k"; do
      Info "Resampling ${morph_data} to ${Surf}"
      for hemi in lh rh; do
          [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
          HEMICAP=$(echo "$hemisphere" | tr [:lower:] [:upper:])
          surf_id="${idBIDS}_hemi-${HEMICAP}"_surf
          Do_cmd wb_command -metric-resample \
              "${outDir}/${surf_id}-fsnative_label-${morph_data}.func.gii" \
              "${dir_conte69}/${surf_id}-fsnative_label-sphere.surf.gii" \
              "${util_surface}/${Surf}.${HEMICAP}.sphere.reg.surf.gii" \
              BARYCENTRIC \
              "${outDir}/${surf_id}-${Surf}_label-${morph_data}.func.gii" #\
      done
    done
  else
      Info "Subject ${idBIDS} cortical ${morph_data} is registered to all surfaces"
  fi
}

reg_surfaces "thickness"
reg_surfaces "curv"

#------------------------------------------------------------------------------#
# End if module has been processed
cleanup "$tmp_morph" "$nocleanup" "$here"
