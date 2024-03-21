#/bin/bash

##################################
### Set variables beneath this ###
##################################

# Permissions
umask 002

# Save OLD PATH
export OLD_PATH=$PATH

#------------------------------------------------------------------------------#
# SOFTWARE CONFIGURATION for MICAPIPE
#------------------------------------------------------------------------------#
# User defined PATHS
# AFNI
export AFNIDIR="/data/mica1/01_programs/afni-20.2.06"
# ANTS
export ANTSPATH="/data/mica1/01_programs/ants-2.3.4/bin"
# Workbench
export workbench_path="/data/mica1/01_programs/workbench-1.4.2/bin_linux64"
# ICA-FIX
export FIXPATH="/data_/mica1/01_programs/fix"
# FreeSurfer
export FREESURFER_HOME="/data/mica1/01_programs/freesurfer-7.3.2"
# fastsurfer
export FASTSURFER_HOME="/data_/mica1/01_programs/fastsurfer"
export fs_licence="/data_/mica1/01_programs/freesurfer-7.3.2/license.txt"
# FSL 6.0
export FSLDIR="/data_/mica1/01_programs/fsl-6-0-3"
export FSL_DIR="/data_/mica1/01_programs/fsl-6-0-3"
export FSL_BIN="${FSLDIR}/bin"
# MRtrix3 3.0.1
export mrtrixDir="/data_/mica1/01_programs/mrtrix3-3.0.1"
# ITK utils
export itk_dir="/data_/mica1/01_programs/c3d-1.0.0-Linux-x86_64/bin"
# Python 3.7
#export PYTHON_3="/data/mica1/01_programs/micapipe-v0.2.0_conda/micapipe/bin"
# Export fs fs_licence
export fs_licence=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
# Fastsurfer singularity container
export fastsurfer_img=/data_/mica1/01_programs/fastsurfer/fastsurfer-cpu-v2.0.4.sif
unset TMPDIR
# Fastsurfer conda env

#------------------------------------------------------------------------------#
# Remove any other instance from the PATH
# AFNI
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*afni*});IFS=':';echo "${p[*]}";unset IFS)
# ANTS
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*ants*});IFS=':';echo "${p[*]}";unset IFS)
# Workbench binaries
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*workbench*});IFS=':';echo "${p[*]}";unset IFS)
# FSL
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*fsl*});IFS=':';echo "${p[*]}";unset IFS)
# revome any other MRtrix3 version from path
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*mrtrix*});IFS=':';echo "${p[*]}";unset IFS)
# REMOVES any other python configuration from the PATH the conda from the PATH and LD_LIBRARY_PATH variable
PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*conda*});IFS=':';echo "${p[*]}";unset IFS)
LD_LIBRARY_PATH=$(IFS=':';p=($LD_LIBRARY_PATH);unset IFS;p=(${p[@]%%*conda*});IFS=':';echo "${p[*]}";unset IFS)

#------------------------------------------------------------------------------#
# Software configuration
# FreeSurfer 6.0 configuration
source "${FREESURFER_HOME}/FreeSurferEnv.sh"
# FSL 6.0 configuration
source "${FSLDIR}/etc/fslconf/fsl.sh"
# PYTHON 3.7 configuration
unset PYTHONPATH
unset PYTHONHOME
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
conda3_bin=/data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/bin/
source /data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/etc/profile.d/conda.sh

#------------------------------------------------------------------------------#
# Set the libraries paths for mrtrx and fsl
export LD_LIBRARY_PATH="${FSLDIR}/lib:${FSL_BIN}:${mrtrixDir}/lib"

#-----------------------------------------------------------------------------------#
# Export new PATH with al the necessary binaries
#export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FIXPATH}:${FREESURFER_HOME}/bin/:${mrtrixDir}/bin:${mrtrixDir}/lib:${FSLDIR}:${FSL_BIN}:${PYTHON_3}:${FASTSURFER_HOME}:${itk_dir}:${PATH}"
export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FIXPATH}:${FREESURFER_HOME}/bin/:${mrtrixDir}/bin:${mrtrixDir}/lib:${FSLDIR}:${FSL_BIN}:${FASTSURFER_HOME}:${itk_dir}:${conda3_bin}:${PATH}"
conda activate /data/mica1/01_programs/micapipe-v0.2.0_conda/micapipe

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Add the number of threads to use here. Note that this is overwritten by
# $NSLOTS if it exists (i.e. when running on SGE).
local_threads="$1"
if [[ -z $local_threads ]]; then export local_threads=10; fi

# Set basic global variables.
# export MICAPIPE="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Note: As this file is sourced by mica-pipe, this will return the mica-pipe path NOT the path of this script.
if [[ ! -z "$NSLOTS" ]]; then
    export threads="$NSLOTS"
else
    export threads="$local_threads"
fi
export OMP_NUM_THREADS="$threads"

# Where processing will run
if [[ -z "$PROC" ]]; then export PROC="LOCAL-MICA"; fi

# Set tmpDir depending on the node
host=$(echo "$HOSTNAME" | awk -F '.' '{print $1}')
case $host in
    fladgate*|yeatman*|oncilla*) export tmpDir="/host/$host/local_raid/temporaryLocalProcessing" ;;
    cassio*|varro*) export tmpDir="/host/$host/export02/data/temporaryLocalProcessing" ;;
    *) export tmpDir="/data/mica2/temporaryNetworkProcessing" ;;
esac

export SGE_ROOT=/opt/sge
