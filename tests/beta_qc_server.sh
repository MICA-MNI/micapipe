#!/bin/bash
set -euo pipefail
[[ $(uname -s) == Linux ]] || { echo "Run on the Linux server." >&2; exit 1; }
root="${1:?Usage: beta_qc_server.sh WORKSPACE [IMAGE]}"
cd "$root"
image="${2:-$PWD/artifacts/micapipe-v1.0.0-beta.2.sif}"
bids="${MICAPIPE_TEST_BIDS:-/data/mica3/BIDS_MICs/rawdata}"
license="${MICAPIPE_TEST_FS_LICENSE:?Set MICAPIPE_TEST_FS_LICENSE to a readable external license}"
test -f "$license"
mkdir -p evidence validation/qc-beta2-HC063 validation/qc-beta2-HC031 tmp/qc-beta2-HC063 tmp/qc-beta2-HC031 tmp/qc-beta2-home
printf 'raise SystemExit(42)\n' > evidence/intentional_capture_failure.py
cp -a --reflink=auto validation/pipeline-HC063-r2/. validation/qc-beta2-HC063/
cp -a --reflink=auto validation/sif-structural-HC031/. validation/qc-beta2-HC031/
find validation/qc-beta2-HC031 -name '*qc-report.pdf' -delete
run_qc() {
 local subject="$1"; shift
 singularity run --cleanenv --containall --env OMP_NUM_THREADS=4,OPENBLAS_NUM_THREADS=4,ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4 \
  --home "$PWD/tmp/qc-beta2-home:$HOME" \
  --bind "$bids:/bids:ro,$PWD/validation/qc-beta2-$subject:/out,$PWD/tmp/qc-beta2-$subject:/tmp,$license:/license.txt:ro" \
  "$@" "$image" -bids /bids -out /out -sub "$subject" -ses 01 -fs_licence /license.txt -threads 4 -QC_subj
}
set +e
run_qc HC063 --bind "$PWD/evidence/intentional_capture_failure.py:/opt/micapipe/functions/nifti_capture.py:ro" > evidence/qc-beta2-intentional-failure.log 2>&1
failure_status=$?
set -e
echo "$failure_status" > evidence/qc-beta2-intentional-failure.status
[[ "$failure_status" != 0 ]]
grep -q 'Subject QC failed' evidence/qc-beta2-intentional-failure.log
run_qc HC063 > evidence/qc-beta2-HC063.log 2>&1
run_qc HC031 > evidence/qc-beta2-HC031.log 2>&1
environments/micapipe/bin/python - <<'PY'
from pathlib import Path
from pypdf import PdfReader
import json
results={}
for subject in ['HC063','HC031']:
 p=next(Path('validation/qc-beta2-'+subject).rglob('*proc_structural_qc-report.pdf'))
 reader=PdfReader(p)
 imgs=[i for page in reader.pages for i in page.images]
 actual=[i for i in imgs if i.image.size!=(2344,234)]
 assert len(actual)>=8,(subject,len(actual))
 assert p.stat().st_size>100000
 results[subject]={'pages':len(reader.pages),'embedded_images':len(imgs),'anatomical_panels':len(actual),'bytes':p.stat().st_size}
results['capture_failure_cli_exit']=int(Path('evidence/qc-beta2-intentional-failure.status').read_text().strip())
results['explicit_license_path_accepted']=True
Path('evidence/qc-beta2-results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results))
PY
