#!/bin/bash

# this re-processes functional scans for the data release
# without bandpass

MICAPIPE=/mica/mica1/03_projects/yigu/software/micapipe
qc_dir=/mica/mica1/03_projects/yigu/7t_funcQC/release
data_name="BIDS_PNI_release"

IFS=$'\n' read -d '' -a subjlist < "${qc_dir}"/data/"${data_name}"_subjects.txt

for subjid in ${subjlist[@]}; do

    for sesid in 01 02 03; do

        idBIDS=sub-${subjid}_ses-${sesid}
        echo "========================reprocessing without filter==================="
        echo $idBIDS

        sh ${MICAPIPE}/7tqc_functest.sh $subjid $sesid

    done
done
echo ":::::::::::::::::::::::done reprocessing, check if ready for tedana QC::::::::::::::::::::"