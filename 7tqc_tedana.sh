#!/bin/bash

# this script resubmits proc_func temp files through tedana
# for testing

source /home/bic/yiguz/.bashrc
micapipe_dev
source $MICAPIPE/functions/init.sh

#-------figure out echo time------#
mainScanStr="task-rest_echo-1_bold,task-rest_echo-2_bold,task-rest_echo-3_bold"
subject_bids="/mica/mica1/03_projects/yigu/7t_funcQC/data/rawdata/sub-PNC003/ses-01"
idBIDS="sub-PNC003_ses-01"
tmp="/mica/mica1/03_projects/yigu/7t_funcQC/tmp/3431_micapipe_proc-func_sub-PNC003_ses-01"

IFS=',' read -ra func_json <<< $mainScanStr
for i in "${!func_json[@]}"; do
    func_json[i]=$(ls "${subject_bids}/func/${idBIDS}_${func_json[$i]}".json 2>/dev/null)
done # Full path
mainScanJson=(${func_json[*]})


unset EchoTime
for json in ${mainScanJson[*]}; do
	EchoTime+=($(grep EchoTime "${json}" | awk -F ' ' '{print $2}' | awk -F ',' '{print $1}'))
done

# retrieve scans for tedana
tedana_dir="${tmp}"/tedana
tmp_func_mask="${tmp}/mainScan01_mean_mask.nii.gz"
scans4tedana=("${tmp}"/mainScan_mc.nii.gz $(ls "${tmp}"/mainScan?_mc.nii.gz))
echo "Files      : ${scans4tedana[*]/${tmp}/}"
echo "EchoNumber : ${EchoNumber[*]}"
#echo "EchoTime   :" "${EchoTime[*]}"

mkdir -p "${tedana_dir}"

tedana -d $(printf "%s " "${scans4tedana[@]}") -e $(printf "%s " "${EchoTime[@]}") --out-dir "${tedana_dir}" --mask "${tmp_func_mask}"
