#!/bin/bash
#
# Compute geodesic distance along cortical mesh:
#
# Generates geodesic distance matrices along native cortical surface using specified parcellations
#
# This workflow makes use of wb_command and custom python scripts
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
nocleanup=$6
threads=$7
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

# Check inputs: freesurfer space T1 (to make sure freesurfer was run)
if [ ! -f ${T1freesurfr} ]; then Error "Subject $id doesn't have a T1 in freesurfer space: <SUBJECTS_DIR>/${id}/mri/T1.mgz"; exit; fi
if [ ! -f ${lh_midsurf} ]; then Error "Subject $id doesn't have left hemisphere midsurface gifti file"; exit; fi
if [ ! -f ${rh_midsurf} ]; then Error "Subject $id doesn't have right hemisphere midsurface gifti file"; exit; fi


#------------------------------------------------------------------------------#
Title "Running MICA Geodesic distance processing"
micapipe_software
bids_print.variables-post
Info "wb_command will use $OMP_NUM_THREADS threads"

#	Timer
aloita=$(date +%s)

#------------------------------------------------------------------------------#
# Set up parameters

# Define output directory
outPath=${dir_surf}/geo_dist/
[[ ! -d "$outPath" ]] && mkdir -p "$outPath"

# wb_command
workbench_path=`which wb_command`

#------------------------------------------------------------------------------#
# Compute geodesic distance on all parcellations
all_parcellations='vosdewael-100 vosdewael-200 vosdewael-300 vosdewael-400
schaefer-100 schaefer-200 schaefer-300 schaefer-400 schaefer-500 schaefer-600 schaefer-700 schaefer-800 schaefer-900 schaefer-1000
glasser-360
economo
aparc
aparc-a2009s'
for parc in ${all_parcellations}; do
    lh_annot=${dir_surf}/${id}/label/lh.${parc}_mics.annot
    rh_annot=${dir_surf}/${id}/label/rh.${parc}_mics.annot
    outName=${outPath}/${id}_${parc}
    Do_cmd python $MICAPIPE/functions/geoDistMapper.py "$lh_midsurf" "$rh_midsurf" "$outName" "$lh_annot" "$rh_annot" "$workbench_path"
    echo completed "$parc"
done

Do_cmd rm -rf ${outPath}/*.func.gii

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Post-GD processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\tlogs:
$dir_logs/post-gd_*.txt"
echo "${id}, post_gd, ${status}, `whoami`, `uname -n`, $(date), `printf "%0.3f\n" ${eri}`, $PROC" >> ${out}/brain-proc.csv
bids_variables_unset
