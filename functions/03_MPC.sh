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
input_dat=$9
mpc_reg=${10}
mpc_str=${11}
PROC=${12}
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
if [[ "$input_im" == "DEFAULT" ]]; then
    Warning "MPC processing will be performed from default input image: qT1"
    microImage="$bids_T1map"
else
    Warning "MPC processing will be performed from provided input file"
    microImage="${input_im}"
fi
Note "Microstructural image =" "$microImage"

# Check microstructural image to registrer
if [[ "$mpc_reg" != "DEFAULT" ]]; then
    Warning "MPC will register the provided input image to the surface: \n\t$mpc_reg"
    regImage="$mpc_reg"
else
    if [[ "$input_im" == "DEFAULT" ]]; then
        Warning "MPC will register the default input image to the surface: inv1_T1map"
        regImage="$bids_inv1"
    else
        Warning "MPC will register the provided input file"
        microImage="$input_im"
    fi
fi
Note "micro2fs registration image =" "$regImage"

# Check .dat file input flag and set parameters accordingly
if [[ "$input_dat" == "DEFAULT" ]]; then
    Warning "Registration to surface space will be performed within script"
else
    Warning "Applying provided .dat file to perform registration to native surface space"
    Note "micro2fs transform =" "$input_dat"
    if [ -z "${input_dat}" ] || [ ! -f "${input_dat}" ]; then Error "Provided registration for MPC was not found or the path is wrong!!!"; exit; fi
fi

# Exit if microImage or Registration image do not exists
if [ -z "${microImage}" ] || [ ! -f "${microImage}" ]; then Error "Image for MPC was not found or the path is wrong!!!"; exit; fi
if [ -z "${regImage}" ] || [ ! -f "${regImage}" ]; then Error "Image for MPC registration was not found or the path is wrong!!!"; exit; fi

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
# If no lta specified by user, register to surface space using T1w as intermediate volume
T1_fsnative=${proc_struct}/${idBIDS}_space-fsnative_T1w.nii.gz
if [[ ${input_dat} == "DEFAULT" ]]; then
    fs_transform="${dir_warp}/${idBIDS}_from-${mpc_str}_to-fsnative_bbr.dat"
    if [[ ! -f "${dir_warp}/${idBIDS}_from-${mpc_str}_to-fsnative_bbr.dat" ]]; then ((N++))
        Info "Running microstructural -> Surface registration for Subject ${id}" #            --int "$T1_fsnative" \
        Do_cmd bbregister --s "$idBIDS" \
            --mov "$regImage" \
            --init-rr \
            --reg "$fs_transform" \
            --o "${subject_dir}/anat/${idBIDS}_space-fsnative_desc-${mpc_str}.nii.gz" \
            --t1
        if [[ -f "${subject_dir}/anat/${idBIDS}_space-fsnative_desc-${mpc_str}.nii.gz" ]]; then ((Nsteps++)); fi
    else
        Info "Subject ${id} already has a microstructural -> Surface space transformation"; ((Nsteps++)); ((N++))
    fi
else
    Info "Using provided input .dat for vol2surf"
    fs_transform="$input_dat"; ((Nsteps++)); ((N++))
fi

##------------------------------------------------------------------------------#
## Register qT1 intensity to surface
num_surfs=14
[[ ! -d "$outDir" ]] && mkdir -p "$outDir" && chmod -R 770 "$outDir"
json_mpc "$microImage" "$fs_transform" "${outDir}/${idBIDS}_MPC-${mpc_str}.json"

if [[ ! -f "${outDir}/${idBIDS}_space-fsnative_desc-rh_MPC-14.mgh" ]]; then
    for hemi in lh rh ; do
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
        for n in $(seq 1 1 "$num_surfs") ; do ((N++))
            which_surf=$(sed -n "$n"p <<< "$x")
            cp "$which_surf" "${dir_subjsurf}/surf/${hemi}.${n}by${num_surf}surf"
            # sample intensity
            Do_cmd mri_vol2surf \
                --mov "$microImage" \
                --reg "$fs_transform" \
                --hemi "$hemi" \
                --out_type mgh \
                --interp trilinear \
                --out "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_MPC-${n}.mgh" \
                --surf "${n}by${num_surf}surf"

            #Remove surfaces used by vol2surf
            Do_cmd rm -rfv "$which_surf" "${dir_subjsurf}/surf/${hemi}.${n}by${num_surf}surf"
            if [[ -f "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_MPC-${n}.mgh" ]]; then ((Nsteps++)); fi
        done
    done
else
    Info "Subject ${id} has microstructural intensities mapped to native surface"; Nsteps=$((Nsteps+28)); N=$((N+28))
fi

# Register to fsa5
for hemi in lh rh; do
    for n in $(seq 1 1 "$num_surfs"); do ((N++))
        MPC_fs5="${outDir}/${idBIDS}_space-fsaverage5_desc-${hemi}_MPC-${n}.mgh"
        if [[ ! -f "$MPC_fs5" ]]; then
            Do_cmd mri_surf2surf --hemi "$hemi" \
                --srcsubject "$idBIDS" \
                --srcsurfval "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_MPC-${n}.mgh" \
                --trgsubject fsaverage5 \
                --trgsurfval "$MPC_fs5"
            if [[ -f "$MPC_fs5" ]]; then ((Nsteps++)); fi
        else
            Info "Subject ${id} ${hemi}_MPC-${n} is registered to fsaverage5"; ((Nsteps++))
        fi
    done
done

# Register to conte69
for hemi in lh rh; do
    [[ "$hemi" == lh ]] && hemisphere=l || hemisphere=r
    HEMICAP=$(echo $hemisphere | tr [:lower:] [:upper:])
    for n in $(seq 1 1 "$num_surfs"); do
        MPC_c69="$outDir/${idBIDS}_space-conte69-32k_desc-${hemi}_MPC-${n}.mgh"
        if [[ ! -f "$MPC_c69" ]]; then ((N++))
            Do_cmd mri_convert "${outDir}/${idBIDS}_space-fsnative_desc-${hemi}_MPC-${n}.mgh" "${tmp}/${hemi}_${n}.func.gii"
            Do_cmd wb_command -metric-resample \
                    "${tmp}/${hemi}_${n}.func.gii" \
                    "${dir_conte69}/${idBIDS}_${hemi}_space-fsnative_desc-sphere.surf.gii" \
                    "${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii" \
                    ADAP_BARY_AREA \
                    "${tmp}/${hemi}_${n}_c69-32k.func.gii" \
                    -area-surfs \
                    "${dir_subjsurf}/surf/${hemi}.midthickness.surf.gii" \
                    "${dir_conte69}/${idBIDS}_space-conte69-32k_desc-${hemi}_midthickness.surf.gii"
            Do_cmd mri_convert "${tmp}/${hemi}_${n}_c69-32k.func.gii" "$MPC_c69"
            if [[ -f "$MPC_c69" ]]; then ((Nsteps++)); fi
        else
            Info "Subject ${id} ${hemi}_MPC-${n} is registered to conte69"; ((Nsteps++)); ((N++))
        fi
    done
done

#------------------------------------------------------------------------------#
# Create MPC connectomes and Intensity profiles per parcellations
parcellations=($(find "$dir_volum" -name "*.nii.gz" ! -name "*cerebellum*" ! -name "*subcortical*"))
for seg in "${parcellations[@]}"; do
    parc=$(echo "${seg/.nii.gz/}" | awk -F 'atlas-' '{print $2}')
    parc_annot="${parc}_mics.annot"
    MPC_int="${outDir}/${idBIDS}_space-fsnative_atlas-${parc}_desc-intensity_profiles.txt"
    if [[ ! -f "$MPC_int" ]]; then ((N++))
        Info "Running MPC on $parc"
        Do_cmd python $MICAPIPE/functions/surf2mpc.py "$out" "$id" "$SES" "$num_surfs" "$parc_annot" "$dir_subjsurf" "${mpc_p}"
        if [[ -f "$MPC_int" ]]; then ((Nsteps++)); fi
    else Info "Subject ${id} MPC connectome and intensity profile on ${parc}"; ((Nsteps++)); ((N++)); fi
done

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
