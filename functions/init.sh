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
# Check if we're running in Docker environment
if [ -f "/opt/freesurfer-7.4.1/SetUpFreeSurfer.sh" ]; then
    # Docker environment paths
    echo "🐳 Detected Docker environment - using Docker paths"
    export AFNIDIR="/opt/afni-23.1.09"
    export ANTSPATH="/opt/ants-2.3.4/bin"
    export workbench_path="/opt/workbench-1.4.2/bin"
    export FIXPATH="/opt/fix"
    export FREESURFER_HOME="/opt/freesurfer-7.4.1"
    export FASTSURFER_HOME="/opt/fastsurfer"
    export fs_licence="/opt/freesurfer-7.4.1/license.txt"
    export FSLDIR="/opt/fsl-6.0.5.1"
    export FSL_DIR="/opt/fsl-6.0.5.1"
    export FSL_BIN="${FSLDIR}/bin"
    export mrtrixDir="/opt/mrtrix3-3.0.1"
    export itk_dir="/opt/c3d/bin"
else
    # Local MICA cluster paths
    echo "🏠 Detected local environment - using MICA cluster paths"
    export AFNIDIR="/data/mica1/01_programs/afni-20.2.06"
    export ANTSPATH="/data/mica1/01_programs/ants-2.3.4/bin"
    export workbench_path="/data/mica1/01_programs/workbench-1.4.2/bin_linux64"
    export FIXPATH="/data_/mica1/01_programs/fix"
    export FREESURFER_HOME="/data/mica1/01_programs/freesurfer-7.3.2"
    export FASTSURFER_HOME="/data_/mica1/01_programs/fastsurfer"
    export fs_licence="/data_/mica1/01_programs/freesurfer-7.3.2/license.txt"
    export FSLDIR="/data_/mica1/01_programs/fsl-6-0-3"
    export FSL_DIR="/data_/mica1/01_programs/fsl-6-0-3"
    export FSL_BIN="${FSLDIR}/bin"
    export mrtrixDir="/data_/mica1/01_programs/mrtrix3-3.0.1"
    export itk_dir="/data_/mica1/01_programs/c3d-1.0.0-Linux-x86_64/bin"
fi
# Python 3.7
#export PYTHON_3="/data/mica1/01_programs/micapipe-v0.2.0_conda/micapipe/bin"
# Export fs fs_licence
if [ -f "/opt/freesurfer-7.4.1/SetUpFreeSurfer.sh" ]; then
    # Docker environment - fs_licence already set above
    export fastsurfer_img="/opt/fastsurfer/fastsurfer.sif"  # Docker path if available
else
    # Local environment
    export fs_licence=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
    export fastsurfer_img=/data_/mica1/01_programs/fastsurfer/fastsurfer-cpu-v2.0.4.sif
fi
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
# FreeSurfer configuration
if [ -f "${FREESURFER_HOME}/SetUpFreeSurfer.sh" ]; then
    source "${FREESURFER_HOME}/SetUpFreeSurfer.sh"
elif [ -f "${FREESURFER_HOME}/FreeSurferEnv.sh" ]; then
    source "${FREESURFER_HOME}/FreeSurferEnv.sh"
else
    echo "⚠️  Warning: FreeSurfer setup script not found at ${FREESURFER_HOME}"
fi

# FSL configuration
if [ -f "${FSLDIR}/etc/fslconf/fsl.sh" ]; then
    source "${FSLDIR}/etc/fslconf/fsl.sh"
else
    echo "⚠️  Warning: FSL setup script not found at ${FSLDIR}/etc/fslconf/fsl.sh"
fi

# PYTHON configuration
unset PYTHONPATH
unset PYTHONHOME
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Conda configuration - only for local environment
if [ ! -f "/opt/freesurfer-7.4.1/SetUpFreeSurfer.sh" ]; then
    conda3_bin=/data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/bin/
    if [ -f "/data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/etc/profile.d/conda.sh" ]; then
        source /data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/etc/profile.d/conda.sh
    fi
fi

#------------------------------------------------------------------------------#
# Set the libraries paths for mrtrx and fsl
export LD_LIBRARY_PATH="${FSLDIR}/lib:${FSL_BIN}:${mrtrixDir}/lib"

#-----------------------------------------------------------------------------------#
# Export new PATH with all the necessary binaries
#export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FIXPATH}:${FREESURFER_HOME}/bin/:${mrtrixDir}/bin:${mrtrixDir}/lib:${FSLDIR}:${FSL_BIN}:${PYTHON_3}:${FASTSURFER_HOME}:${itk_dir}:${PATH}"
if [ -f "/opt/freesurfer-7.4.1/SetUpFreeSurfer.sh" ]; then
    # Docker environment
    export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FIXPATH}:${FREESURFER_HOME}/bin/:${mrtrixDir}/bin:${mrtrixDir}/lib:${FSLDIR}:${FSL_BIN}:${FASTSURFER_HOME}:${itk_dir}:${PATH}"
else
    # Local environment
    export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FIXPATH}:${FREESURFER_HOME}/bin/:${mrtrixDir}/bin:${mrtrixDir}/lib:${FSLDIR}:${FSL_BIN}:${FASTSURFER_HOME}:${itk_dir}:${conda3_bin}:${PATH}"
    if [ -d "/data/mica1/01_programs/micapipe-v0.2.0_conda/micapipe" ]; then
        conda activate /data/mica1/01_programs/micapipe-v0.2.0_conda/micapipe
    fi
fi

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
