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
PROC=$5
input_im=$6
input_lta=$7
nocleanup=$8
threads=$9
export OMP_NUM_THREADS=$threads

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

# Check inputs: Freesurfer space T1
if [ ! -f ${T1freesurfr} ]; then Error "T1 in freesurfer space not found for Subject $id : <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi

# Check microstructural image input flag and set parameters accordingly
if [[ ${input_im} == "DEFAULT" ]]; then
    Warning "MPC processing will be performed from default input image: qT1"
    Note "qT1 =" ${bids_T1map}
    microImage=${bids_T1map}
    regImage=${bids_inv1}
else
    Warning "MPC processing will be performed from provided input file"
    Note "Microstructural image =" ${input_im}
    microImage=${input_im}
fi

# Check .lta file input flag and set parameters accordingly
if [[ ${input_lta} == "DEFAULT" ]]; then
    Warning "Registration to freesurfer space will be performed within script"
else
    Warning "Applying provided .lta file to perform registration to native freesurfer space"
    Note "micro2fs transform =" ${input_lta}
fi

#------------------------------------------------------------------------------#
Title "Running MICA MPC processing"
micapipe_software
bids_print.variables-post
Info "Not erasing temporal dir: $nocleanup"
Info "wb_command will use $OMP_NUM_THREADS threads"

#	Timer
aloita=$(date +%s)

# if temporary directory is empty
if [ -z ${tmp} ]; then tmp=/tmp; fi
# Create temporary directory
tmp=${tmp}/${RANDOM}_micapipe_post-MPC_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# TRAP in case the script fails
trap cleanup INT TERM

# Freesurface SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Temporary fsa5 directory
ln -s $FREESURFER_HOME/subjects/fsaverage5/ ${dir_surf}


#------------------------------------------------------------------------------#
# If no lta specified by user, register to Freesurfer space using T1w as intermediate volume

origImage=${bids_T1ws[0]}
if [[ ${input_lta} == "DEFAULT" ]]; then
    fs_transform="$dir_warp"/"$id"_micro2fs_bbr.lta
    if [[ ! -f "$dir_warp"/"$id"_micro2fs_bbr.lta ]]; then
        Info "Running microstructural -> freesurfer registration for Subject ${id}"
        Do_cmd bbregister --s "$id" \
            --mov "$regImage" \
            --int "$origImage" \
            --reg "$fs_transform" \
            --o "$subject_dir"/proc_struct/"$id"_micro2fsspace.nii.gz \
            --init-header --t1
    else
        Info "Subject ${id} already has a microstructural -> freesurfer transformation"
    fi
else
    Info "Using provided input .lta for vol2surf"
    fs_transform=${input_lta}
fi


##------------------------------------------------------------------------------#
## Register qT1 intensity to surface

num_surfs=14
outDir="$subject_dir"/proc_struct/surfaces/micro_profiles/
[[ ! -d "$outDir" ]] && mkdir -p "$outDir"

if [[ ! -f ${outDir}/rh.14.mgh ]]; then
    for hemi in lh rh ; do

        unset LD_LIBRARY_PATH
        tot_surfs=$((num_surfs + 2))
        Do_cmd python $MICAPIPE/functions/generate_equivolumetric_surfaces.py \
            ${SUBJECTS_DIR}/${id}/surf/${hemi}.pial \
            ${SUBJECTS_DIR}/${id}/surf/${hemi}.white \
            "$tot_surfs" \
            $outDir/${hemi}.${num_surfs}surfs \
            ${tmp} \
            --software freesurfer --subject_id $id

        # remove top and bottom surface
        Do_cmd rm -rfv ${outDir}/${hemi}.${num_surfs}surfs0.0.pial ${outDir}/${hemi}.${num_surfs}surfs1.0.pial

        # find all equivolumetric surfaces and list by creation time
        x=$(ls -t "$outDir"/"$hemi".${num_surfs}surfs*)
        for n in $(seq 1 1 $num_surfs) ; do
            which_surf=$(sed -n "$n"p <<< "$x")
            cp "$which_surf" "$SUBJECTS_DIR"/"$id"/surf/"$hemi"."$n"by"$num_surf"surf
            # sample intensity
            Do_cmd mri_vol2surf \
                --mov "$microImage" \
                --reg "$fs_transform" \
                --hemi "$hemi" \
                --out_type mgh \
                --interp trilinear \
                --out "$outDir"/"$hemi"."$n".mgh \
                --surf "$n"by"$num_surf"surf

            #Remove surfaces used by vol2surf
            Do_cmd rm -rfv "$which_surf" "$SUBJECTS_DIR"/"$id"/surf/"$hemi"."$n"by"$num_surf"surf
        done

    done
else
    Info "Subject ${id} has microstructural intensities mapped to native surface"
fi

# Register to fsa5
if [[ ! -f ${outDir}/rh.14_fsa5.mgh ]]; then
    for hemi in lh rh; do
        for n in $(seq 1 1 $num_surfs); do
            Do_cmd mri_surf2surf --hemi "$hemi" \
                --srcsubject "$id" \
                --srcsurfval "$outDir"/"$hemi"."$n".mgh \
                --trgsubject fsaverage5 \
                --trgsurfval "$outDir"/"$hemi"."$n"_fsa5.mgh
        done
    done
else
    Info "Subject ${id} microstructural intensities are registered to fsa5"
fi

# Register to conte69
if [[ ! -f ${outDir}/rh_14_c69-32k.mgh ]]; then
    for hemi in lh rh; do
        [[ $hemi == lh ]] && hemisphere=l || hemisphere=r
        HEMICAP=`echo $hemisphere | tr [:lower:] [:upper:]`
        for n in $(seq 1 1 $num_surfs); do
            Do_cmd mri_convert "$outDir"/"$hemi"."$n".mgh "$tmp"/"$hemi"_"$n".func.gii

            Do_cmd wb_command -metric-resample \
                ${tmp}/${hemi}_"$n".func.gii \
                ${dir_conte69}/${id}_${hemi}_sphereReg.surf.gii \
                ${util_surface}/fs_LR-deformed_to-fsaverage.${HEMICAP}.sphere.32k_fs_LR.surf.gii \
                ADAP_BARY_AREA \
                ${tmp}/${hemi}_"$n"_c69-32k.func.gii \
                -area-surfs \
                ${dir_surf}/${id}/surf/${hemi}.midthickness.surf.gii \
                ${dir_conte69}/${id}_${hemi}_midthickness_32k_fs_LR.surf.gii

            Do_cmd mri_convert ${tmp}/${hemi}_"$n"_c69-32k.func.gii "$outDir"/${hemi}_"$n"_c69-32k.mgh
         done
    done
else
    Info "Subject ${id} microstructural intensities are registered to conte69"
fi

#------------------------------------------------------------------------------#
# run  mpc on native surface
all_parcellations='vosdewael-100 vosdewael-200 vosdewael-300 vosdewael-400
schaefer-100 schaefer-200 schaefer-300 schaefer-400 schaefer-500 schaefer-600 schaefer-700 schaefer-800 schaefer-900 schaefer-1000
glasser-360
economo
aparc
aparc-a2009s'
for parc in ${all_parcellations}; do
    parc_annot=${parc}_mics.annot
    Do_cmd python $MICAPIPE/functions/surf2mpc.py "$out" "$id" "$SES" "$num_surfs" "$parc_annot"
    echo completed "$parc"
done

#------------------------------------------------------------------------------#
# Clean temporary directory and fsaverage5
if [[ $nocleanup == "FALSE" ]]; then Do_cmd rm -rf $tmp ${dir_surf}/fsaverage5; else Info "Mica-pipe tmp directory was not erased: \n\t\t\t${tmp}"; fi

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Post-MPC processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\tlogs:
$dir_logs/post-mpc_*.txt"
echo "${id}, post_mpc, ${status}, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
bids_variables_unset
