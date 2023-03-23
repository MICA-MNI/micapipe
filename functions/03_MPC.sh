#!/bin/bash
#
# Microstructural imaging processing:
#
# Preprocessing workflow for qT1.
# Generates microstructural profiles and mpc matrices on specified parcellations
#
# This workflow makes use of freesurfer and custom python scripts
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
input_im=$8
mpc_reg=$9
mpc_str=${10}
PROC=${11}
export OMP_NUM_THREADS=$threads
here=$(pwd)

#------------------------------------------------------------------------------#
# qsub configuration
if [ "$PROC" = "qsub-MICA" ] || [ "$PROC" = "qsub-all.q" ] || [ "$PROC" = "LOCAL-MICA" ]; then
    export MICAPIPE=/data_/mica1/01_programs/micapipe-v0.2.0
    source "${MICAPIPE}/functions/init.sh" "$threads"
fi

# source utilities
source $MICAPIPE/functions/utilities.sh

# Assigns variables names
bids_variables "$BIDS" "$id" "$out" "$SES"

# Check dependencies Status: POST_STRUCTURAL
micapipe_check_dependency "post_structural" "${dir_QC}/${idBIDS}_module-post_structural.json"

# Setting Surface Directory from post_structural
post_struct_json="${proc_struct}/${idBIDS}_post_structural.json"
recon=$(grep SurfaceProc ${post_struct_json} | awk -F '"' '{print $4}')
set_surface_directory "${recon}"

# Variables naming for multiple acquisitions
if [[ "${mpc_str}" == DEFAULT ]]; then
  mpc_str="qMRI"
  mpc_p="acq-qMRI"
else
  mpc_p="acq-${mpc_str}"
fi

# End if module has been processed
module_json="${dir_QC}/${idBIDS}_module-MPC-${mpc_str}.json"
micapipe_check_json_status "${module_json}" "MPC"

# Check microstructural image input flag and set parameters accordingly
if [[ "$input_im" == "DEFAULT" ]]; then microImage="$bids_T1map"; else microImage="${input_im}"; fi
Note "Microstructural image =" "$microImage"

# Check microstructural image to registrer
if [[ "$mpc_reg" == "DEFAULT" ]]; then regImage="${bids_inv1}"; else regImage="${mpc_reg}"; fi
Note "Microstructural image for registration =" "$regImage"

# Exit if microImage or Registration image do not exists
if [ ! -f "${microImage}" ]; then Error "Image for MPC was not found or the path is wrong!!!"; exit; fi
if [ ! -f "${regImage}" ]; then Error "Image for MPC registration was not found or the path is wrong!!!"; exit; fi

#------------------------------------------------------------------------------#
Title "Microstructural Profiles Covariance\n\t\tmicapipe $Version, $PROC"
micapipe_software
bids_print.variables-post
Note "Saving temporal dir:" "$nocleanup"
Note "Temporal dir:" "${tmp}"
Note "Parallel processing:" "$threads threads"

#	Timer
aloita=$(date +%s)
Nsteps=0
N=0

# Create script specific temp directory
tmp="${tmpDir}/${RANDOM}_micapipe_post-MPC_${id}"
Do_cmd mkdir -p "$tmp"

# TRAP in case the script fails
trap 'cleanup $tmp $nocleanup $here' SIGINT SIGTERM

# Freesurface SUBJECTs directory
export SUBJECTS_DIR="$dir_surf"
outDir="${subject_dir}/mpc/${mpc_p}"
Note "acqMRI:" "${mpc_str}"

#------------------------------------------------------------------------------#
# Affine registration between both images
T1_in_fs=${tmp}/orig.mgz
qT1_fsnative=${proc_struct}/${idBIDS}_space-fsnative_${mpc_str}.nii.gz
mat_fsnative_affine=${dir_warp}/${idBIDS}_from-fsnative_to_${mpc_str}_
qT1_fsnative_affine=${mat_fsnative_affine}0GenericAffine.mat

if [[ ! -f "$qT1_fsnative" ]] || [[ ! -f "$qT1_fsnative_affine" ]]; then ((N++))
    Do_cmd mrconvert "$T1surf" "$T1_in_fs"
    Do_cmd antsRegistrationSyN.sh -d 3 -f "$regImage" -m "$T1_in_fs" -o "$mat_fsnative_affine" -t a -n "$threads" -p d
    Do_cmd antsApplyTransforms -d 3 -i "$regImage" -r "$T1_in_fs" -t ["${qT1_fsnative_affine}",1] -o "$qT1_fsnative" -v -u int
    if [[ -f ${qT1_fsnative} ]]; then ((Nsteps++)); fi
else
    Info "Subject ${id} has a ${mpc_str} on Surface space"; ((Nsteps++)); ((N++))
fi

# Convert the ANTs transformation file for wb_command
wb_affine="${tmp}/${idBIDS}_from-fsnative_to_qMRI_wb.mat"
Do_cmd c3d_affine_tool -itk $qT1_fsnative_affine -inv -o ${wb_affine}

##------------------------------------------------------------------------------#
## Register qT1 intensity to surface
num_surfs=14
[[ ! -d "$outDir" ]] && mkdir -p "$outDir" && chmod -R 770 "$outDir"
json_mpc "$microImage" "${outDir}/${idBIDS}_MPC-${mpc_str}.json"

MPC_fsLR5k="${outDir}/${idBIDS}_surf-fsLR-5k_desc-MPC.shape.gii"
if [[ ! -f "${MPC_fsLR5k}" ]]; then ((N++))
    for hemi in lh rh ; do
        [[ "$hemi" == lh ]] && HEMI=L || HEMI=R
        unset LD_LIBRARY_PATH
        tot_surfs=$((num_surfs + 2))
        Do_cmd python "$MICAPIPE"/functions/generate_equivolumetric_surfaces.py \
            "${dir_subjsurf}/surf/${hemi}.pial" \
            "${dir_subjsurf}/surf/${hemi}.white" \
            "$tot_surfs" \
            "${outDir}/${hemi}.${num_surfs}surfs" \
            "$tmp" \
            --software freesurfer --subject_id "$idBIDS"

        # remove top and bottom surface
        Do_cmd rm -rfv "${outDir}/${hemi}.${num_surfs}surfs0.0.pial" "${outDir}/${hemi}.${num_surfs}surfs1.0.pial"

        # find all equivolumetric surfaces and list by creation time
        x=$(ls -t "$outDir"/"$hemi".${num_surfs}surfs*)
        for n in $(seq 1 1 "$num_surfs") ; do
            which_surf=$(sed -n "$n"p <<< "$x")
            surf_gii="${tmp}/${hemi}.${n}by${num_surf}_space-fsnative.surf.gii"
            surf_tmp="${tmp}/${hemi}.${n}by${num_surf}_no_offset.surf.gii"
            out_surf="${tmp}/${hemi}.${n}by${num_surf}_space-qMRI.surf.gii"
            out_feat="${outDir}/${idBIDS}_hemi-${HEMI}_surf-fsnative_label-MPC-${n}.func.gii"
            # Register surface to qMRI space
            Do_cmd mris_convert "$which_surf" "${surf_gii}"
            # Remove offset-to-origin from any gifti surface derived from FS
            Do_cmd python "$MICAPIPE"/functions/removeFSoffset.py "${surf_gii}" "${surf_tmp}"
            # Apply transformation to register surface to nativepro
            Do_cmd wb_command -surface-apply-affine "${surf_tmp}" "${wb_affine}" "${out_surf}"
            # Sample intensity and resample to other surfaces
            map_to-surfaces "${microImage}" "${out_surf}" "${out_feat}" "${HEMI}" "MPC-${n}" "${outDir}"
            # remove tmp surfaces
            rm "${surf_tmp}" "${which_surf}"
        done
    done
    ((Nsteps++))
else
    Info "Subject ${id} has microstructural intensities mapped to native surface";((Nsteps++)); ((N++));
fi

#------------------------------------------------------------------------------#
# Create MPC connectomes and Intensity profiles per parcellations
parcellations=($(find "$dir_volum" -name "*atlas*" ! -name "*cerebellum*" ! -name "*subcortical*"))
for seg in "${parcellations[@]}"; do
    parc=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    parc_annot="${parc}_mics.annot"
    MPC_int="${outDir}/${idBIDS}_atlas-${parc}_desc-intensity_profiles.txt"
    if [[ ! -f "$MPC_int" ]]; then ((N++))
        Info "Running MPC on $parc"
        Do_cmd python $MICAPIPE/functions/surf2mpc.py "$out" "$id" "$SES" "$num_surfs" "$parc_annot" "$dir_subjsurf" "${mpc_p}"
        if [[ -f "$MPC_int" ]]; then ((Nsteps++)); fi
    else Info "Subject ${id} has MPC connectome and intensity profile on ${parc}"; ((Nsteps++)); ((N++)); fi
done

#------------------------------------------------------------------------------#
# Create vertex-wise MPC connectome and directory cleanup
if [[ ! -f "${MPC_fsLR5k}" ]]; then ((N++))
  Info "Running MPC vertex-wise on fsLR-5k"
  Do_cmd python $MICAPIPE/functions/build_mpc-vertex.py "$out" "$id" "$SES" "${mpc_p}"
  ((Nsteps++))
else Info "Subject ${id} has MPC vertex-wise on fsLR-5k"; ((Nsteps++)); ((N++)); fi

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print "$eri"/60 | perl)

# Notification of completition
micapipe_completition_status "MPC"
micapipe_procStatus "${id}" "${SES/ses-/}" "MPC-${mpc_str}" "${out}/micapipe_processed_sub.csv"
Do_cmd micapipe_procStatus_json "${id}" "${SES/ses-/}" "MPC-${mpc_str}" "${module_json}"
cleanup "$tmp" "$nocleanup" "$here"
