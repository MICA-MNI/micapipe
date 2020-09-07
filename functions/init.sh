#/bin/bash

##################################
### Set variables beneath this ###
##################################

# Add your bin paths here.
afni_path=/data/mica1/01_programs/afni-20.2.06/bin
ants_path=/data/mica1/01_programs/ants-2.3.4/bin
mrtrix_path=/data/mica1/01_programs/mrtrix3-micapipe/bin
workbench_path=/data/mica1/01_programs/workbench/bin_linux64
freesurfer_path=/data/mica1/01_programs/Freesurfer-6.0/bin
PYTHON_PATH=/data_/mica1/01_programs/anaconda/anaconda3/envs/mica_py3.7
FIXPATH=/data_/mica1/01_programs/fix

# FSL donfiguration file
. /data_/mica1/01_programs/fsl_mica/etc/fslconf/fsl.sh

# Add the number of threads to use here. Note that this is overwritten by
# $NSLOTS if it exists (i.e. when running on SGE).
local_threads=10

# Uncomment this and fill in a temporary directory for a custom temporary directory.
# This takes priority over the default.
global_temp_directory=/data/mica2/temporaryNetworkProcessing/

# Uncomment this and fill in host host/temporary directories for
# custom temporary directories for every host. Note that this takes
# priority over global_temp_directory. hostnames and temporary directories should
# be separated by spaces (not commas!). When running on an unspecified host,
# the program will default to the global_temp_directory or default temp directory.
host_temp_dirs=(fladgate.bic.mni.mcgill.ca /host/fladgate/local_raid/temporaryLocalProcessing/ \
                yeatman.bic.mni.mcgill.ca /host/yeatman/local_raid/temporaryLocalProcessing/ \
                cassio.bic.mni.mcgill.ca /host/cassio/export02/data/temporaryLocalProcessing/ \
                oncilla.bic.mni.mcgill.ca /host/oncilla/local_raid/temporaryLocalProcessing)

# Default temporary directory
tmp_file=$(mktemp)
default_temp=$(dirname $tmp_file)
rm -f $tmp_file

# Set basic global variables.
export MICAPIPE="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Note: As this file is sourced by mica-pipe, this will return the mica-pipe path NOT the path of this script.
export OLD_PATH=$PATH
export PATH=${MICAPIPE}:${script_path}:${afni_path}:${ants_path}:${FIXPATH}:${mrtrix_path}:${workbench_path}:${freesurfer_path}:${PYTHON_PATH}:${PATH}
if [[ ! -z $NSLOTS ]]; then
    export CORES=$NSLOTS
else
    export CORES=$local_threads
fi
export OMP_NUM_THREADS=$CORES

# Set the temporary directory
hostname=$(uname -n)
# First try setting from the host specific directories.
for idx in $(seq 0 2 ${#host_temp_dirs[@]}); do
    idx2=$(echo "$idx + 1" | bc)
    if [[ $hostname == ${host_temp_dirs[$idx]} ]]; then
        export tmp=${host_temp_dirs[$idx2]}
        break
    fi
done
# If that didn't work, try setting from the global/default instead.
if [[ -z $tmp ]]; then
    if [[ ! -z $global_temp_directory ]]; then
        export tmp=$global_temp_directory
    else
        export tmp=$default_temp
    fi
fi
