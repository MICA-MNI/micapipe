#!/bin/bash

# test micapipe functional processing
# with changed micapipe script as to save tedana outputs elsewhere
# MAKE SURE PROC SURF WAS RUN BEFORE
# MAKE SURE TO MAKE YOURSELF A MICAPIPE DEV -

source ~/.bashrc
micapipe_dev # this function exists in ~/.bashrc

testpath=/mica/mica1/03_projects/yigu/7t_funcQC/release # a micapipe directory will be created here which is great
tmpDir="${testpath}"/tmp_remask
bidsDir=/data/mica3/BIDS_PNI/data_release/rawdata
#bidsDir="${testpath}"/data/rawdata


#------set the subject and session MANUALLY-----#
subjid="PNC002"
sesid="02"

#-------------call proc func----------#
MICAPIPE=/mica/mica1/03_projects/yigu/software/micapipe

cd $MICAPIPE

./micapipe -sub $subjid \
    -ses $sesid \
    -out $testpath \
    -bids $bidsDir \
    -tmpDir $tmpDir \
    -proc_func \
    -mainScanStr "task-rest_echo-1_bold,task-rest_echo-2_bold,task-rest_echo-3_bold" \
    -nocleanup \
    -mica # this is required as environment.

echo "---------------> DONE PROC-FUNC FOR RUN"
echo "               sub-${subjid}_ses-${sesid}"
