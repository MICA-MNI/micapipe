#!/bin/bash

# Define the variables
SUBJECT_ID="sub-mri3T"
OUTPUT_DIR="/home/bic/eyang/processed_data"
BIDS_DIR="/data_/mica3/BIDS_CI/rawdata"
FS_LIC="$FREESURFER_HOME/license.txt"
SES="01"
T1W_STR="run-2_T1w"  # Example: processing only 'run-02_T1w'
MULTI_FACTOR=45  # Example factor for 3T MP2RAGE, adjust if different

# Run micapipe
micapipe -sub $SUBJECT_ID -out $OUTPUT_DIR -bids $BIDS_DIR \
         -proc_structural -T1wStr $T1W_STR -mf $MULTI_FACTOR \
	 -fs_licence $FS_LIC -ses $SES

