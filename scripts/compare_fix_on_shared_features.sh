#!/bin/bash
# Isolate model-conversion differences using identical legacy feature tables.
set -euo pipefail
root="$1" subject="$2"
export FSLDIR="$root/environments/fsl-6.0.7.8" FSLOUTPUTTYPE=NIFTI_GZ
export PATH="$FSLDIR/bin:$FSLDIR/share/fsl/bin:$PATH" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
base="$root/validation/fix/$subject"
out="$base/classifier-only.ica"
test ! -e "$out"
cp -a "$base/common.ica" "$out"
mkdir -p "$out/fix"
cp "$base/legacy.ica/fix/features.csv" "$out/fix/features.csv"
cd "$base"
"$FSLDIR/bin/fix" -c "$out" "$root/source/functions/MICAMTL_training_15HC_15PX.pyfix_model" 20 -lf "$out/pyfix.log" -nj 2
labels="$out/fix4melview_MICAMTL_training_15HC_15PX_thr20.txt"
"$FSLDIR/bin/fix" -a "$labels" "$out" -m -h 100 -lf "$out/pyfix.log" -nj 2
"$root/environments/micapipe/bin/python" "$root/source/scripts/compare_fix.py" "$base/legacy.ica" "$out" --output "$base/classifier-only-comparison.json"
