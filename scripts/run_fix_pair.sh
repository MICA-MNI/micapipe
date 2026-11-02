#!/bin/bash
set -euo pipefail
root="${1:?workspace}" subject="${2:?subject}"
export FSLDIR="$root/environments/fsl-6.0.7.8"
export PATH="$FSLDIR/share/fsl/bin:$PATH" FSLOUTPUTTYPE=NIFTI_GZ
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
out="$root/validation/fix/$subject"
test -f "$out/prepare.exit"
cd "$out"
# Share one registration for both feature extractors.
common="$out/common.ica"
flirt -in "$common/reg/highres" -ref "$FSLDIR/data/standard/MNI152_T1_1mm_brain"     -omat "$common/reg/highres2standard_1mm.mat" -dof 12
convert_xfm -omat "$common/reg/standard2highres.mat" -inverse "$common/reg/highres2standard_1mm.mat"
convert_xfm -omat "$common/reg/standard2example_func.mat" -concat     "$common/reg/highres2example_func.mat" "$common/reg/standard2highres.mat"
for backend in legacy python; do
    test ! -e "$out/$backend.ica" || { echo "Refusing to reuse comparison outputs" >&2; exit 1; }
    cp -a "$out/common.ica" "$out/$backend.ica"
done
"$root/downloads/fix-legacy/fix" "$out/legacy.ica" \
    "$root/source/functions/MICAMTL_training_15HC_15PX.RData" 20 -m -h 100
"$FSLDIR/bin/fix" "$out/python.ica"     "$root/source/functions/MICAMTL_training_15HC_15PX.pyfix_model" 20 -m -h 100 -lf "$out/python.ica/pyfix.log" -nj 2
"$root/environments/micapipe/bin/python" "$root/source/scripts/compare_fix.py" \
    "$out/legacy.ica" "$out/python.ica" --output "$out/comparison.json"
echo 0 > "$out/pair.exit"
