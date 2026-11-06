#!/bin/bash
# Run as the calling user in containers/os/Dockerfile with /opt on a data disk.
# MICAPIPE_VARIANT picks the PyTorch build of FastSurfer and FastScrubber:
# nocuda (CPU only) or cuda (NVIDIA GPUs; also runs on the CPU).
set -euo pipefail
cd /build
variant="${MICAPIPE_VARIANT:?Set MICAPIPE_VARIANT to nocuda or cuda}"
case "$variant" in
    nocuda) fastsurfer_env=fastsurfer_cpu ;;
    cuda) fastsurfer_env=fastsurfer_cuda ;;
    *) echo "Unknown MICAPIPE_VARIANT: $variant" >&2; exit 1 ;;
esac
export UV_CACHE_DIR=/build/cache/uv-runtime UV_LINK_MODE=copy
export UV_PYTHON_INSTALL_DIR=/opt/python UV_PYTHON_DOWNLOADS=automatic
install -m 755 tools/uv /opt/uv
/opt/uv python install 3.9.20 3.11.14 3.14.0
mkdir -p /opt/environments
for env in micapipe tedana lamareg designer "$fastsurfer_env" fastscrubber workflow; do
    py=3.11.14; lock="$env"; args=()
    case "$env" in
        fastsurfer_cpu) py=3.14.0; lock=fastsurfer-cpu; args+=(--torch-backend cpu) ;;
        fastsurfer_cuda) py=3.14.0; lock=fastsurfer-cuda; args+=(--torch-backend cu126) ;;
        fastscrubber)
            py=3.9.20
            if [[ "$variant" == cuda ]]; then lock=fastscrubber-cuda; args+=(--torch-backend cu124)
            else lock=fastscrubber-cpu; args+=(--torch-backend cpu); fi ;;
    esac
    if [[ ! -x "/opt/environments/$env/bin/python" ]]; then /opt/uv venv --python "$py" "/opt/environments/$env"; fi
    /opt/uv pip sync --python "/opt/environments/$env/bin/python" "${args[@]}" \
        "source/environments/locks/$lock-linux-64.txt"
    /opt/uv pip check --python "/opt/environments/$env/bin/python"
done
MAMBA=/build/tools/micromamba
export MAMBA_ROOT_PREFIX=/build/cache/mamba-runtime
export CONDA_PKGS_DIRS=/build/cache/conda-runtime
if [[ ! -f /opt/fsl-6.0.7.23/.micapipe-ready ]]; then
    "$MAMBA" --no-rc create -y -p /opt/fsl-6.0.7.23 -f source/environments/locks/fsl-6.0.7.23-linux-64.txt
    touch /opt/fsl-6.0.7.23/.micapipe-ready
fi
printf '6.0.7.23' > /opt/fsl-6.0.7.23/etc/fslversion
if [[ ! -f /opt/environments/mrtrix3/.micapipe-ready ]]; then
    "$MAMBA" --no-rc create -y -p /opt/environments/mrtrix3 -f source/environments/locks/mrtrix3-linux-64.txt
    touch /opt/environments/mrtrix3/.micapipe-ready
fi
[[ -d /opt/FastSurfer ]] || cp -a downloads/FastSurfer /opt/FastSurfer
[[ -d /opt/fast_scrubber ]] || cp -a downloads/fast_scrubber /opt/fast_scrubber
[[ -d /opt/mrtrix3 ]] || cp -a downloads/mrtrix3 /opt/mrtrix3
export PATH="/opt/environments/mrtrix3/bin:$PATH"
if [[ ! -f /opt/mrtrix3/.micapipe-ready ]]; then
    (cd /opt/mrtrix3 && ./configure -nogui && NUMBER_OF_PROCESSORS=8 ./build)
    touch /opt/mrtrix3/.micapipe-ready
fi
mkdir -p /opt/ants-2.6.5
unzip -oq downloads/ants-2.6.5.zip -d /opt/ants-2.6.5
# The release zip has a single top-level directory.
if [[ ! -d /opt/ants-2.6.5/bin ]]; then
    nested=(/opt/ants-2.6.5/*)
    mv "${nested[0]}"/* /opt/ants-2.6.5/
    rmdir "${nested[0]}"
fi
mkdir -p /opt/afni /opt/freesurfer-package
tar -xzf downloads/afni-linux_openmp_64.tgz -C /opt/afni --strip-components=1
dpkg-deb -x downloads/freesurfer_ubuntu22-8.2.0_amd64.deb /opt/freesurfer-package
fsroot=$(find /opt/freesurfer-package -name SetUpFreeSurfer.sh -print -quit)
test -n "$fsroot"
ln -sfn "${fsroot%/SetUpFreeSurfer.sh}" /opt/freesurfer
[[ -d /opt/c3d-1.0.0-Linux-x86_64 ]] || cp -a /build/downloads/c3d-1.0.0-Linux-x86_64 /opt/
mkdir -p /opt/micapipe
cp -a source/. /opt/micapipe/
rm -rf /opt/micapipe/.git
mkdir -p /opt/micapipe/bin /opt/fix-legacy
tar -xzf downloads/fix-1.068.tar.gz -C /opt/fix-legacy --strip-components=1
# Preserve the older feature extraction/classification backend for comparison.
cat > /opt/fix-legacy/settings.sh <<'SETTINGS'
FIXVERSION=1.06
export FSL_FIXDIR=/opt/fix-legacy
export FSL_FIX_MATLAB_MODE=2
export FSL_FIX_OCTAVE=/opt/fix-legacy/octave-wrapper
export FSL_FIX_OCEVAL="--eval"
export FSL_FIX_OCFILE=""
export FSL_FIX_FSLMATLAB="${FSLDIR}/etc/matlab"
export FSL_FIX_OCOPTS="--traditional -q --no-window-system"
export FSL_FIX_R_CMD=/usr/bin/R
export FSL_FIX_R_LIBS=/usr/lib/R/site-library
SETTINGS
cp -a source/fsl_conf/. /opt/fsl-6.0.7.23/etc/flirtsch/
cp -a source/surfaces/fsaverage5 /opt/freesurfer/subjects/
cp source/containers/legacy_octave.sh /opt/fix-legacy/octave-wrapper
chmod +x /opt/fix-legacy/octave-wrapper
# FSL bin/share wrappers isolate its Python interpreter from pipeline Python.
for spec in 'lamareg lamareg' 'tedana tedana' 'designer designer' 'designer tmi' 'micapipe pyonlm'; do
    read -r env command <<< "$spec"
    printf '#!/bin/bash\nsource /opt/micapipe/functions/micapipe_env.sh\nmicapipe_env %q %q "$@"\n' \
        "$env" "$command" > "/opt/micapipe/bin/$command"
done
cat > /opt/micapipe/bin/fastscrubber <<'WRAPPER'
#!/bin/bash
source /opt/micapipe/functions/micapipe_env.sh
export PYTHONNOUSERSITE=1
micapipe_env fastscrubber python /opt/fast_scrubber/fastscrubber.py "$@"
WRAPPER
cp source/containers/fix_wrapper.sh /opt/micapipe/bin/fix
cat > /opt/micapipe/bin/micapipe-niivue <<'WRAPPER'
#!/bin/bash
exec /opt/environments/micapipe/bin/python /opt/micapipe/functions/niivue_report.py "$@"
WRAPPER
chmod +x /opt/micapipe/bin/* /opt/micapipe/micapipe /opt/micapipe/functions/*.sh /opt/micapipe/functions/*.py
# Make DESIGNER use this exact MRtrix Python library without a global PYTHONPATH.
printf '/opt/mrtrix3/lib\n' > /opt/environments/designer/lib/python3.11/site-packages/mrtrix3.pth
mkdir -p /opt/environments/lamareg/lib/python3.11/site-packages/lamareg/scripts/models
cp -a /build/downloads/lamareg-models/. /opt/environments/lamareg/lib/python3.11/site-packages/lamareg/scripts/models/
# FastScrubber weights: staged copy, checked (and completed if needed) by its own SHA-256 downloader.
mkdir -p /opt/fast_scrubber/weights
if [[ -d downloads/fast_scrubber-weights ]]; then cp -a downloads/fast_scrubber-weights/. /opt/fast_scrubber/weights/; fi
/opt/environments/fastscrubber/bin/python /opt/fast_scrubber/download_weights.py --dest /opt/fast_scrubber/weights
/opt/micapipe/bin/fastscrubber --help > /dev/null
# Fetch the FastSurfer checkpoints at build time, using its released downloader.
cd /opt/FastSurfer
"/opt/environments/$fastsurfer_env/bin/python" -m FastSurferCNN.download_checkpoints --all
# run_fastsurfer.sh reads its version from BUILD.info when it exists, and falls back
# to git otherwise. /opt/FastSurfer belongs to the build user, not the caller, so at
# run time git refuses it: every FastSurfer call printed four "dubious ownership"
# errors and logged an empty version. Write the cache here, where git can read it.
git config --global --add safe.directory /opt/FastSurfer
"/opt/environments/$fastsurfer_env/bin/python" -m FastSurferCNN.version \
    --sections +git+checkpoints -o /opt/FastSurfer/BUILD.info
grep -q . /opt/FastSurfer/BUILD.info
