#!/bin/bash
set -euo pipefail
source /opt/micapipe/containers/runtime_env.sh
source /opt/micapipe/functions/micapipe_env.sh
micapipe -h >/tmp/micapipe-help.txt || [ "$?" -eq 1 ]
grep -q proc_structural /tmp/micapipe-help.txt
mrconvert --version
antsRegistration --version
flirt -version
command -v topup eddy_cpu
mri_convert --version
afni -ver
wb_command -version
dcm2niix --version || [ "$?" -eq 3 ]
/opt/fsl-6.0.7.8/bin/python -c 'from pyfix import fix; print("Python FIX imports passed")'
micapipe_env fastsurfer_cpu python -c 'import torch, scipy, nibabel; print("FastSurfer Python:", torch.__version__)'
micapipe_env lamareg lamareg --help >/tmp/lamareg-help.txt
micapipe_env designer designer -help >/tmp/designer-help.txt
micapipe_env tedana tedana --help >/tmp/tedana-help.txt
pyonlm --help >/tmp/pyonlm-help.txt
micapipe-niivue --help >/tmp/niivue-help.txt
python - <<'PY'
import nibabel,numpy,scipy,skfmm,astropy,pygeodesic
from brainspace.plotting import plot_hemispheres
from xhtml2pdf import pisa
print("Processing and QC imports passed")
PY

# Initialization must preserve the pipeline interpreter after native tools are added.
(
 # FreeSurfer's setup script expects ordinary shell options.
 set +eu
 source /opt/micapipe/functions/init.sh 1 >/tmp/micapipe-init.txt 2>&1
 init_status=$?
 set -eu
 [[ "$init_status" == 0 ]]
 python -c 'import sys,nibabel,brainspace; assert sys.prefix == "/opt/environments/micapipe", sys.prefix; print("Pipeline Python isolation passed")'
 before=$(command -v python)
 micapipe_env fastsurfer_cpu python -c 'import sys; assert sys.version_info[:2] == (3,14)'
 [[ "$before" == "$(command -v python)" ]]
)

micapipe_env workflow python -c 'from snakebids import bidsapp,plugins; import snakemake; print("Workflow imports passed")'

# Exercise real PNG generation: imports alone do not catch missing QC panels.
python - <<'PYQC'
import numpy as np,nibabel as nb
from pathlib import Path
d=np.zeros((32,32,32),dtype=np.float32);d[8:24,8:24,8:24]=1
nb.save(nb.Nifti1Image(d,np.eye(4)),'/tmp/qc-smoke.nii.gz')
PYQC
/opt/micapipe/functions/nifti_capture.py -img /tmp/qc-smoke.nii.gz -out /tmp/qc-smoke.png >/tmp/qc-smoke.log
python - <<'PYQC'
from PIL import Image
from pathlib import Path
p=Path('/tmp/qc-smoke.png')
assert p.stat().st_size>1000
assert min(Image.open(p).size)>100
print('QC image generation passed')
PYQC
