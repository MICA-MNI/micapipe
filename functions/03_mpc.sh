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
#   $4 : Temporal directory (default /tmp)
#
# ONLY for scripting and debugging:
# TEST=ON
# source utilities
source $MICAPIPE/functions/utilities.sh

BIDS=$1
id=$2
out=$3
tmp=$4

#------------------------------------------------------------------------------#
Title "Running MICA MPC processing"

# Assigns variables names
bids_variables $BIDS $id $out
# print the names on the terminal
bids_print.variables-post

# GLOBAL variables for this script
Info "wb_command will use $OMP_NUM_THREADS threads"

# Check inputs: mp2rage
if [ ! -f ${qT1} ]; then Error "Subject $id doesn't have qT1"; exit; fi

# Check inputs: freesurfer space T1
if [ ! -f ${T1freesurfr} ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi

#	Timer
aloita=$(date +%s)
here=`pwd`

# Check tmp dir: temporary directory
if [ -z ${tmp} ]; then tmp=/tmp; fi
tmp=${tmp}/${RANDOM}_post-MPC_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Freesurface SUBJECTs directory
export SUBJECTS_DIR=${dir_surf}

# Temporary fsa5 directory
ln -s $FREESURFER_HOME/subjects/fsaverage5/ ${dir_surf}

#------------------------------------------------------------------------------#
# Set up parameters
num_surfs=14
origImage=${bids_T1ws[0]}
microImage=${bids_T1map[0]}
invImage=${bids_inv1[0]}

#------------------------------------------------------------------------------#
# Register to Freesurfer space

if [[ ! -f "$dir_warp"/"$id"_T1map2fs_bbr.lta ]] ; then
	Do_cmd bbregister --s "$id" \
        --mov "$invImage" \
        --int "$origImage" \
        --reg "$dir_warp"/"$id"_T1map2fs_bbr.lta \
        --o "$subject_dir"/proc_struct/"$id"_inv1_T1map_fsspace.nii \
        --init-header --t1
fi


#------------------------------------------------------------------------------#
# Register qT1 intensity to surface

outDir="$subject_dir"/proc_struct/surfaces/micro_profiles/
[[ ! -d "$outDir" ]] && mkdir -p "$outDir"

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
	x=$(ls "$outDir"/"$hemi".${num_surfs}surfs* | sort)
	for n in $(seq 1 1 $num_surfs) ; do
		which_surf=$(sed -n "$n"p <<< "$x")
		cp "$which_surf" "$SUBJECTS_DIR"/"$id"/surf/"$hemi"."$n"by"$num_surf"surf
		# sample intensity
		Do_cmd mri_vol2surf \
			--mov "$microImage" \
			--reg "$dir_warp"/"$id"_T1map2fs_bbr.lta \
			--hemi "$hemi" \
			--out_type mgh \
			--interp trilinear \
			--out "$outDir"/"$hemi"."$n".mgh \
			--surf "$n"by"$num_surf"surf

        #Remove surfaces used by vol2surf
        Do_cmd rm -rfv "$which_surf" "$SUBJECTS_DIR"/"$id"/surf/"$hemi"."$n"by"$num_surf"surf
	done

done

# Register to fsa5
for hemi in lh rh; do
    for n in $(seq 1 1 $num_surfs); do
        Do_cmd mri_surf2surf --hemi "$hemi" \
            --srcsubject "$id" \
            --srcsurfval "$outDir"/"$hemi"."$n".mgh \
            --trgsubject fsaverage5 \
            --trgsurfval "$outDir"/"$hemi"."$n"_fsa5.mgh
    done
done

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
    Do_cmd python $MICAPIPE/functions/surf2mpc.py "$out" "$id" "$num_surfs" "$parc_annot"
    echo completed "$parc"
done

#------------------------------------------------------------------------------#
# Clean temporary directory and fsaverage5
Do_cmd rm -rfv $tmp ${dir_surf}/fsaverage5

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Post-MPC processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\t\t\tlogs:${dir_logs}/post_MPC.txt"
echo "${id}, post_MPC, TEST, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
