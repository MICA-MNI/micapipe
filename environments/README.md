# V1 beta environments

Build and test on Linux servers. The container comes in two variants named
v1-beta-nocuda-DATE and v1-beta-cuda-DATE (date in containers/versions.json);
they differ only in the PyTorch build of FastSurfer and FastScrubber. The
previous candidate was v1.0.0-beta.2.
Pipeline source starts at upstream v1 commit 6df741b7a2f6ddb03b302350ca6539422c9669ad.

| Environment | Purpose |
| --- | --- |
| micapipe | Python 3.11.14, processing, QC, SWM, pyonlm |
| lamareg | Python 3.11.14, TensorFlow 2.15.1/Keras 2.15 and ANTsPy |
| designer | Python 3.11.14, DESIGNER's independently pinned dependencies |
| tedana | Python 3.11.14, existing tedana 0.0.12 algorithm |
| fastsurfer_cpu | Python 3.14.0, FastSurfer 2.6.0 / PyTorch 2.14 CPU (nocuda variant) |
| fastsurfer_cuda | Python 3.14.0, FastSurfer 2.6.0 / PyTorch 2.14 CUDA 12.6 (cuda variant) |
| fastscrubber | Python 3.9.20, FastScrubber (nnU-Net 2.5.1) / PyTorch 2.4.1 CPU or CUDA 12.4, for -alt_fs |
| mrtrix3 | Workbench and dcm2niix; MRtrix3 3.0.8 built from its tag |
| workflow | Python 3.11.14, SnakeBIDS 0.13.1 and Snakemake 8.30.0 |
| FSL | Official FSL 6.0.7.8 packages, Python FIX 0.6.0 |

locks/ contains resolved Python dependencies or explicit conda package URLs
with checksums. Git dependencies use full commit IDs. FSL needs
setuptools=80.9.0: its file-tree=1.3.0 dependency imports pkg_resources.
Do not update FSL's Python dependencies through the micapipe environment.

micapipe_env.sh activates an environment in a subshell. Command wrappers
select LAMAReg, DESIGNER, tedana and pyonlm without altering the caller's Python
environment. FastSurfer uses the same scoped activation.

To add a dependency, edit only the owning .in file, resolve it with the
pinned uv version (0.12.22) and Python interpreter, update its lock, and test
the consuming command and affected pipeline module on the server. Native
dependencies and versions belong in containers/os/Dockerfile.

## Build

Stage pinned downloads described by containers/versions.json under a dedicated
workspace containing source/, downloads/, and tools/. FastScrubber is staged as
downloads/fast_scrubber (a `git archive` of the pinned commit) and, optionally,
downloads/fast_scrubber-weights/fast_scrubber/nnUNet_results; the build checks
the weights' SHA-256 and downloads any that are missing from OSF. Verify
containers/SOURCES.sha256 from that workspace. Then run, once per variant and
in separate workspaces:

    bash source/containers/build_beta_server.sh /path/to/workspace nocuda
    bash source/containers/build_beta_server.sh /path/to/workspace cuda

Each build writes artifacts/micapipe-VERSION.sif and
artifacts/micapipe-VERSION.docker.tar.gz, a Docker image archive packed without
a Docker daemon (`docker load -i` it). The cuda variant needs an NVIDIA driver
on the host (`singularity run --nv` or `docker run --gpus all`) and falls back to
the CPU without one.

The CUDA FastSurfer lock keeps every version of the CPU lock and changes only
the PyTorch builds: it was resolved with the CPU lock (minus torch and
torchvision) as a constraint file and `--torch-backend cu126`. PyTorch 2.14.1
has no CUDA 12.8 build.

The script uses ordinary Docker containers with a workspace mount and a
separate /opt mount on the data disk. It exports the small OS filesystem and
packs it with the installed tools into a Singularity image. It does not replace
a production container or publish an image.

## FIX migration and validation

The .pyfix_model is produced by the official command:

    convert_fix_model -d extracted-training MICAMTL_training_15HC_15PX.RData MICAMTL_training_15HC_15PX.pyfix_model

R is needed for conversion. The beta retains R and Octave and the original FIX
backend for comparison. The candidate defaults to legacy because all three MICs comparisons failed the
combined compatibility gate. See [the validation report](../docs/validation/v1-beta.md).
MICAPIPE_FIX_BACKEND=legacy selects the original
.RData model; python selects the converted model. Custom models must match
the selected backend.

prepare_fix_mics.sh regenerates common MELODIC inputs for paired model validation
from MICs raw resting-state data and existing structural derivatives. Both
backends receive copies of the same input. This standalone test uses FSL motion
correction, masking and 100-second high-pass filtering; it does not constitute
an end-to-end reproduction of micapipe's AFNI preprocessing.

compare_fix.py reports component-label agreement, noise-label Jaccard, and
cleaned data agreement after voxelwise temporal demeaning. The provisional
engineering thresholds are 95% label agreement, Jaccard 0.90, correlation 0.999
and normalized RMSE 0.01. These thresholds are declared before the comparison
and do not establish scientific equivalence. Review component maps and cleaned
time series before promoting the model or removing legacy software.

## NiiVue

    micapipe-niivue --volume scan.nii.gz --output viewer/

This creates a local viewer with bundled NiiVue 0.69.0 and a copy of the volume.
Serve the directory over HTTP on localhost. The viewer does not need a CDN or
upload the scan. The module is decoded from the official package's URI-encoded
esm export; package integrity and the bundle checksum are recorded.
