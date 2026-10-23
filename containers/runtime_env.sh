#!/bin/bash
if [ "${MICAPIPE_RUNTIME_INITIALIZED:-}" = 1 ]; then return 0; fi
export MICAPIPE_RUNTIME_INITIALIZED=1
export MICAPIPE=/opt/micapipe MICAPIPE_ENV_ROOT=/opt/environments
export MICAPIPE_FIX_BACKEND="${MICAPIPE_FIX_BACKEND:-legacy}"
export PROC=container_micapipe-v1.0.0-beta.2
export FREESURFER_HOME=/opt/freesurfer FASTSURFER_HOME=/opt/FastSurfer
export FSLDIR=/opt/fsl-6.0.7.8 FSL_DIR=/opt/fsl-6.0.7.8
export FSLOUTPUTTYPE=NIFTI_GZ FSLMULTIFILEQUIT=TRUE
export ANTSPATH=/opt/ants-2.6.5/bin
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
export XDG_CACHE_HOME="/tmp/micapipe-cache"
export MPLCONFIGDIR="${XDG_CACHE_HOME}/matplotlib"
export NUMBA_CACHE_DIR="${XDG_CACHE_HOME}/numba"
export PATH="/opt/micapipe/bin:/opt/micapipe:/opt/micapipe/functions:/opt/environments/micapipe/bin:/opt/ants-2.6.5/bin:/opt/afni:/opt/mrtrix3/bin:/opt/environments/mrtrix3/bin:/opt/fsl-6.0.7.8/share/fsl/bin:/opt/fsl-6.0.7.8/bin:/opt/freesurfer/bin:/opt/FastSurfer:/opt/c3d-1.0.0-Linux-x86_64/bin:${PATH}"
export BASH_ENV=/opt/micapipe/containers/runtime_env.sh
