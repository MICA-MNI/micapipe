# Containers: where the images are, and how to rebuild them

The scripts in this directory build the `v1` beta containers. They do not run in
CI — the images are built by hand on a Linux server and published to the lab's
shared Singularity directory on BIC.

## Where the newest container is

```
/data_/mica1/01_programs/singularity/micapipe-v1.0.0-beta.3-nocuda-20261103.sif   (24.5 GiB, CPU only)
/data_/mica1/01_programs/singularity/micapipe-v1.0.0-beta.3-cuda-20261103.sif     (30.5 GiB, needs --nv)
```

Also reachable as `/data/mica1/01_programs/singularity/`. Take `nocuda` unless
the host has an NVIDIA GPU. Beside the images:

- `micapipe-v1.0.0-beta.3.README.txt` — checksums, run examples, what is new,
  what was tested and what the limits are
- `micapipe-v1.0.0-beta.3.validation.md` — the full build and acceptance record
- `*.sif.sha256` — one per image

```bash
cd /data_/mica1/01_programs/singularity && sha256sum -c micapipe-v1.0.0-beta.3-nocuda-20261103.sif.sha256
```

beta.3 was built from commit `a81231d`, which reached `v1` in
[MICA-MNI/micapipe#189](https://github.com/MICA-MNI/micapipe/pull/189), so the
published images match `v1`. Every image records its own provenance:

```bash
singularity exec IMAGE.sif cat /opt/micapipe/BUILDINFO.json
singularity exec IMAGE.sif cat /opt/micapipe/CONTAINER_VERSION
```

[`versions.json`](versions.json) in this directory pins what goes into the
image: release name and date, the two variants' torch backends, and the version
or commit of FreeSurfer, FastSurfer, FastScrubber, MRtrix3, ANTs, FSL, pyonlm
and NiiVue.

## Rebuilding

### What the host needs

Docker, `singularity`, `mksquashfs` and `pigz`, on Linux — the build script
refuses to run anywhere else. beta.3 was built on `bb-compxg-01`
(Singularity-CE 4.1.2, squashfs-tools 4.5, pigz 2.6, Docker 29.1.3, 128 cores,
503 GB RAM). Nothing needs root, but the build user needs Docker access.

Budget about 200 GB of disk per variant: 18 GB of downloads, ~55 GB of unpacked
tool tree, ~40 GB of resolver cache and ~60 GB for the two artifacts. The
downloads can be shared between variants with hard links; everything else is
per-variant.

### The workspace

`build_beta_server.sh` takes a workspace directory and a variant, and expects
three things to already be in place:

```
WORKSPACE/
  source/      a micapipe git checkout, at the commit you want to build
  downloads/   the upstream payloads (see below)
  tools/       uv and micromamba, extracted from their archives in downloads/
```

It creates `runtime/`, `cache/`, `tmp/`, `evidence/` and `artifacts/` itself.
The beta workspaces on BIC are `/export02/enning/micapipe-v1-beta-20261028/nocuda`
and `.../cuda`, with their build driver scripts in the directory above them.

**The fast path is to reuse one of those.** Everything in `downloads/`, and the
resolved environments under `runtime/`, are already there:

```bash
cd /export02/enning/micapipe-v1-beta-20261028/nocuda
git -C source fetch https://github.com/MICA-MNI/micapipe.git BRANCH
git -C source checkout --force --detach FETCH_HEAD
rm -rf runtime/opt/micapipe          # the code tree is copied in, not marker-guarded
export MICAPIPE_BUILD_CPUS=32 MICAPIPE_BUILD_MEMORY=96g
bash source/containers/build_beta_server.sh "$PWD" nocuda
```

That is what `build-beta3c.sh` in the parent directory does for both variants at
once, including pre-building the shared OS base image. A rebuild that only
changes micapipe's own code took 15–16 minutes per variant. A first build on an
empty workspace takes hours: every Python environment is resolved from its lock,
FSL and MRtrix3 are installed from conda locks, MRtrix3 is compiled, and the
FastSurfer and FastScrubber checkpoints are downloaded.

### Filling `downloads/` from scratch

There is no fetch script; the tree is assembled once and reused. Items marked
*checksummed* are verified by `sha256sum -c source/containers/SOURCES.sha256`,
which the build script runs from the workspace root before anything else.

| Path under `downloads/` | What | Source |
| --- | --- | --- |
| `freesurfer_ubuntu22-8.2.0_amd64.deb` | FreeSurfer 8.2.0 *(checksummed)* | `https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/8.2.0/freesurfer_ubuntu22-8.2.0_amd64.deb` |
| `ants-2.6.5.zip` | ANTs 2.6.5 *(checksummed)* | `ants-2.6.5-ubuntu-22.04-X64-gcc.zip` from the [v2.6.5 release](https://github.com/ANTsX/ANTs/releases/tag/v2.6.5), renamed |
| `afni-linux_openmp_64.tgz` | AFNI *(checksummed)* | `https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz`, renamed |
| `fix-1.068.tar.gz` | legacy FIX 1.068 *(checksummed)* | FSL's FIX distribution. The old `fmrib.ox.ac.uk/~steve/ftp` link is gone, so copy the one already in the workspace |
| `uv.tar.gz` | uv 0.12.22 *(checksummed)* | `https://github.com/astral-sh/uv/releases/download/0.12.22/uv-x86_64-unknown-linux-gnu.tar.gz`, renamed |
| `micromamba-2.5.0-2.tar.bz2` | micromamba 2.5.0-2 *(checksummed)* | `https://github.com/mamba-org/micromamba-releases/releases/download/2.5.0-2/micromamba-linux-64.tar.bz2`, renamed |
| `fast_scrubber/` | FastScrubber at `cec34a1` | https://github.com/zihuaihuai/fast_scrubber — `fastscrubber.py` and `download_weights.py` are *checksummed* |
| `FastSurfer/` | FastSurfer 2.6.0 at `1dc0247` | https://github.com/Deep-MI/FastSurfer.git |
| `mrtrix3/` | MRtrix3 3.0.8 at `464173c` | https://github.com/MRtrix3/mrtrix3.git |
| `c3d-1.0.0-Linux-x86_64/` | Convert3D 1.0.0, unpacked | https://sourceforge.net/projects/c3d/files/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz |
| `vendor-src/lamareg/` | LAMAReg checkout | https://github.com/MICA-MNI/LAMAReg.git — only `lamareg/SynthSeg/model_handler.py` is read, for the model URLs |
| `lamareg-models/` | SynthSeg weights | written by `stage_lamareg_models.py`, pinned to LAMAR-Models `5f43003`; the build script runs it |
| `fast_scrubber-weights/` | optional nnU-Net weight cache | optional: `download_weights.py` fetches or completes the weights at build time from https://osf.io/x95g7 |

Then extract the two tool archives so `tools/uv` and `tools/micromamba` exist —
`install_runtime.sh` runs those binaries directly rather than the archives.

Everything else comes from the repository: the Python environments from
`environments/locks/*.txt`, FSL 6.0.7.23 and MRtrix3's dependencies from their
conda locks, the two vendored `.deb`s in [`os/vendor`](os/vendor), and
micapipe's own code from `source/`.

### What the scripts do

| Script | Role |
| --- | --- |
| `build_beta_server.sh` | the entry point: `build_beta_server.sh WORKSPACE nocuda\|cuda` |
| `os/Dockerfile` | the OS base image, tagged `micapipe-v1-beta-os:20261022` |
| `install_runtime.sh` | runs inside that image with `/opt` bind-mounted to the data disk, and installs every tool |
| `prepare_singularity_metadata.sh`, `runtime_env.sh` | the image's startup hooks and environment |
| `package_beta.sh` | squashes `runtime/` into the `.sif` and smoke-checks it |
| `package_docker.sh` | writes the same filesystem as a Docker archive, with no daemon |
| `docker_entrypoint.sh` | the Docker image's entrypoint |
| `fix_wrapper.sh`, `legacy_octave.sh` | the FIX wrappers; legacy FIX with R and Octave is still the default |
| `export_conda_lock.py` | regenerates a conda lock from an installed prefix |
| `stage_lamareg_models.py` | downloads the pinned SynthSeg weights |

The variant decides only the torch backends: `nocuda` installs CPU FastSurfer
and FastScrubber, `cuda` installs cu126 and cu124. The `cuda` image runs on a
CPU-only host too; it is just larger.

### Things that will bite you

- **The version stamp comes from `versions.json`**, as
  `release-variant-date` — beta.3 produced
  `micapipe-v1.0.0-beta.3-nocuda-20261103.sif`. `package_beta.sh` and
  `package_docker.sh` **refuse to overwrite an existing artifact**, so either
  bump `release`/`date` for a new candidate or move the previous pass's
  artifacts aside first. Rebuilding the same stamp also means anything already
  running from that SIF should be finished before you start.
- **Reruns are incremental**, guarded by `.micapipe-ready` marker files under
  `runtime/opt`. That is what makes a rebuild 15 minutes. micapipe's own code is
  the exception: it is copied in unguarded, so `rm -rf runtime/opt/micapipe`
  before rebuilding, or you ship the old code.
- `MICAPIPE_BUILD_CPUS` (default 12) and `MICAPIPE_BUILD_MEMORY` (default 48g)
  cap the install container and `mksquashfs`.
- The build writes `BUILDINFO.json` and `CONTAINER_VERSION` into
  `/opt/micapipe`, and the source commit into the SIF's labels, from the
  `source/` checkout's `HEAD` — so build from a real commit, not a dirty tree.

### Outputs

```
WORKSPACE/artifacts/micapipe-VERSION.sif
WORKSPACE/artifacts/micapipe-VERSION.docker.tar.gz      # docker load -i ...
WORKSPACE/evidence/container-VERSION.sha256
WORKSPACE/evidence/container-inspect.json
```

## Checking a fresh build

The in-image smoke test is the first thing to run. It checks the pipeline's
help; MRtrix, ANTs, FSL, FreeSurfer, AFNI, Connectome Workbench and dcm2niix;
that `topup`, `eddy_cpu` and `eddy_cuda` are all present and that the FSL
version is the pinned one; that the GPUs are hidden by default and that
`restore_gpus` hands them back; that Python FIX imports; that FastSurfer and
FastScrubber have torch with the expected CUDA support and that the FastScrubber
weights are in place; that LAMAReg, DESIGNER, tedana, pyonlm and the NiiVue
report all answer `--help`; that `init.sh` leaves the pipeline interpreter
isolated; and that a QC PNG is really generated rather than merely imported. It
takes about a minute.

```bash
singularity exec --cleanenv IMAGE.sif bash /opt/micapipe/tests/beta_smoke.sh
```

On a GPU host, run it again with `--nv` on the `cuda` image: that is the mode
where the GPU eddy check actually has something to check.

Static and smoke checks are not enough on their own — both beta.3 defects
(FastSurfer's multi-GPU `DataParallel` hang, and `proc_dwi` never reaching
`eddy_cuda`) passed every static check and only showed up in a real run. Run a
pilot subject through the modules you care about, and run the GPU one with every
GPU on the host visible. The beta.3 acceptance scripts are on BIC under
`/export02/enning/micapipe-acceptance-20261103/scripts`, and
[`docs/validation/v1-beta-3.md`](../docs/validation/v1-beta-3.md) records what
they covered and what they did not.

## Publishing to the shared directory

The shared directory sits at about 98% full, so check `df -h /data_/mica1`
first; a beta pair is roughly 55 GB. Publishing means copying each `.sif` in,
writing its `.sha256`, and putting a README and the validation record beside
them — `publish_beta3.sh` in the beta workspace's parent directory does that for
beta.3 and is the easiest thing to adapt. Verify at the destination afterwards,
not only at the source:

```bash
cd /data_/mica1/01_programs/singularity && sha256sum -c micapipe-v1.0.0-beta.3-*.sif.sha256
```
