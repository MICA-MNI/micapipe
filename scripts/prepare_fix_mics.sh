#!/bin/bash
# Standalone paired FIX validation inputs; writes only under the output root.
set -euo pipefail
root="${1:?workspace}" subject="${2:?subject without sub-}" bids="${3:?BIDS root}" derivatives="${4:?derivative root}"
export FSLDIR="$root/environments/fsl-6.0.7.8" FSLOUTPUTTYPE=NIFTI_GZ
export PATH="$FSLDIR/share/fsl/bin:$PATH"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
output="$root/validation/fix/$subject"
mkdir -p "$output/common.ica"/{reg,mc}
raw="$bids/sub-$subject/ses-01/func/sub-${subject}_ses-01_task-rest_acq-AP_bold.nii.gz"
meta="${raw%.nii.gz}.json"
highres="$derivatives/sub-$subject/ses-01/anat/sub-${subject}_ses-01_space-nativepro_T1w_brain.nii.gz"
[[ -f "$raw" && -f "$highres" ]] || { echo "Missing validation inputs for $subject" >&2; exit 1; }
tr=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["RepetitionTime"])' "$meta")
sigma=$(python3 -c 'import sys;print(50/float(sys.argv[1]))' "$tr")
printf '%s\n' "$raw" "$highres" > "$output/input_paths.txt"
sha256sum "$raw" "$meta" "$highres" > "$output/inputs.sha256"
ica="$output/common.ica"
fslroi "$raw" "$output/drop5" 5 -1
mcflirt -in "$output/drop5" -out "$output/motion" -plots
cp "$output/motion.par" "$ica/mc/prefiltered_func_data_mcf.par"
fslmaths "$output/motion" -Tmean "$output/mean"
bet "$output/mean" "$ica/reg/example_func" -f 0.3 -m
cp "$ica/reg/example_func_mask.nii.gz" "$ica/mask.nii.gz"
fslmaths "$output/motion" -bptf "$sigma" -1 -add "$output/mean" "$ica/filtered_func_data"
fslmaths "$ica/filtered_func_data" -Tmean "$ica/mean_func"
cp "$highres" "$ica/reg/highres.nii.gz"
flirt -in "$ica/reg/highres" -ref "$ica/reg/example_func" -omat "$ica/reg/highres2example_func.mat" -dof 6
flirt -in "$ica/reg/highres" -ref "$FSLDIR/data/standard/MNI152_T1_2mm_brain" -omat "$ica/reg/highres2standard.mat" -dof 12
melodic --in="$ica/filtered_func_data" --tr="$tr" --nobet --mask="$ica/mask" \
  --bgthreshold=3 --mmthresh=0.5 --report --Oall --seed=20261022 \
  --outdir="$ica/filtered_func_data.ica"
echo 0 > "$output/prepare.exit"
