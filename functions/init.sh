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
# Detect Docker environment by the presence of /opt/micapipe (which the
# container always has) plus a containerized FreeSurfer install. Tested
# against the comprehensive base image (Dockerfile.mamba-base): the actual
# install paths there are unversioned (e.g. /opt/freesurfer, not
# /opt/freesurfer-7.4.1) so we probe for the directory itself.
if [ -d "/opt/micapipe" ] && [ -d "/opt/freesurfer" ]; then
    # Docker environment paths — match the layout produced by
    # Dockerfile.mamba-base.
    echo "Detected Docker environment - using container paths"
    export AFNIDIR="/opt/afni"
    export ANTSPATH="/opt/ants-2.3.4/bin"
    # MRtrix3 ships in the micapipe conda env, not as a standalone /opt dir.
    export mrtrixDir="/opt/miniconda-latest/envs/micapipe"
    # connectome-workbench: present in the conda env if installed there;
    # otherwise expected on PATH via apt's connectome-workbench package.
    if [ -x "/opt/miniconda-latest/envs/micapipe/bin/wb_command" ]; then
        export workbench_path="/opt/miniconda-latest/envs/micapipe/bin"
    else
        export workbench_path="/usr/bin"
    fi
    export FIXPATH="/opt/fix"
    export FREESURFER_HOME="/opt/freesurfer"
    export FASTSURFER_HOME="/opt/FastSurfer"
    export fs_licence="/opt/licence.txt"
    export FSLDIR="/opt/fsl-6.0.2"
    export FSL_DIR="/opt/fsl-6.0.2"
    export FSL_BIN="${FSLDIR}/bin"
    export itk_dir="/opt/c3d-1.0.0-Linux-x86_64/bin"
    # Activate the micapipe conda env so mrconvert/mrinfo/etc. resolve.
    if [ -f "/opt/miniconda-latest/etc/profile.d/conda.sh" ]; then
        source /opt/miniconda-latest/etc/profile.d/conda.sh
        conda activate micapipe 2>/dev/null || true
    fi
    export fastsurfer_img=""  # No bundled SIF in Docker; use FastSurfer dir directly
else
    # Local MICA cluster paths
    echo "Detected local environment - using MICA cluster paths"
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
    export fastsurfer_img=/data_/mica1/01_programs/fastsurfer/fastsurfer-cpu-v2.0.4.sif
fi
unset TMPDIR

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
# revome any other MRtrix3 version from path — but only outside Docker,
# where conda's mrtrix3 lives in the active env's bin dir.
if [ ! -d "/opt/micapipe" ]; then
    PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*mrtrix*});IFS=':';echo "${p[*]}";unset IFS)
    # REMOVES any other python configuration from the PATH
    PATH=$(IFS=':';p=($PATH);unset IFS;p=(${p[@]%%*conda*});IFS=':';echo "${p[*]}";unset IFS)
    LD_LIBRARY_PATH=$(IFS=':';p=($LD_LIBRARY_PATH);unset IFS;p=(${p[@]%%*conda*});IFS=':';echo "${p[*]}";unset IFS)
fi

#------------------------------------------------------------------------------#
# Software configuration
# FreeSurfer configuration
if [ -f "${FREESURFER_HOME}/SetUpFreeSurfer.sh" ]; then
    source "${FREESURFER_HOME}/SetUpFreeSurfer.sh"
elif [ -f "${FREESURFER_HOME}/FreeSurferEnv.sh" ]; then
    source "${FREESURFER_HOME}/FreeSurferEnv.sh"
else
    echo "Warning: FreeSurfer setup script not found at ${FREESURFER_HOME}"
fi

# FSL configuration
if [ -f "${FSLDIR}/etc/fslconf/fsl.sh" ]; then
    source "${FSLDIR}/etc/fslconf/fsl.sh"
else
    echo "Warning: FSL setup script not found at ${FSLDIR}/etc/fslconf/fsl.sh"
fi

# PYTHON configuration
unset PYTHONPATH
unset PYTHONHOME
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Conda configuration - only for local environment (Docker activated above)
if [ ! -d "/opt/micapipe" ]; then
    conda3_bin=/data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/bin/
    if [ -f "/data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/etc/profile.d/conda.sh" ]; then
        source /data/mica1/01_programs/micapipe-v0.2.0_conda/conda3/etc/profile.d/conda.sh
    fi
fi

#------------------------------------------------------------------------------#
# Set the libraries paths for mrtrx and fsl
export LD_LIBRARY_PATH="${FSLDIR}/lib:${FSL_BIN}:${mrtrixDir}/lib:${LD_LIBRARY_PATH:-}"

#-----------------------------------------------------------------------------------#
# Export new PATH with all the necessary binaries
if [ -d "/opt/micapipe" ]; then
    # Docker environment — mrtrixDir is the conda env root, so its bin is
    # ${mrtrixDir}/bin (already on PATH after `conda activate micapipe`).
    export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FIXPATH}:${FREESURFER_HOME}/bin/:${mrtrixDir}/bin:${FSLDIR}:${FSL_BIN}:${FASTSURFER_HOME}:${itk_dir}:${PATH}"
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
if [[ ! -z "$NSLOTS" ]]; then
    export threads="$NSLOTS"
else
    export threads="$local_threads"
fi
export OMP_NUM_THREADS="$threads"

# Where processing will run
if [[ -z "$PROC" ]]; then export PROC="LOCAL-MICA"; fi

# Set tmpDir depending on the node — but in Docker, /tmp is the bind mount.
if [ -d "/opt/micapipe" ]; then
    export tmpDir="${tmpDir:-/tmp}"
else
    host=$(echo "$HOSTNAME" | awk -F '.' '{print $1}')
    case $host in
        fladgate*|yeatman*|oncilla*) export tmpDir="/host/$host/local_raid/temporaryLocalProcessing" ;;
        cassio*|varro*) export tmpDir="/host/$host/export02/data/temporaryLocalProcessing" ;;
        *) export tmpDir="/data/mica2/temporaryNetworkProcessing" ;;
    esac
fi

export SGE_ROOT=/opt/sge
