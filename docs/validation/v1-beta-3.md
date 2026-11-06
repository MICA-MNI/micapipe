# V1 beta.3 rebuild validation

Built and checked 2026-11-03 on a BIC Linux server (128 cores, 8x RTX A6000,
driver 580.95.05). This is still an engineering beta for lab review; the GPU
path is now exercised end to end, but nothing here is scientific validation.

beta.3 exists for one reason: beta.2 (20261023) predates both the experimental
`-alt_fs` surfaces and the FSL 6.0.7.23 bump that gives us a usable
`eddy_cuda`, so nothing on `v1` had a matching image.

## What was built

| | |
| --- | --- |
| Source commit | `a81231d` (branch `enning/v1-beta3-rebuild`) |
| Variants | `micapipe-v1.0.0-beta.3-nocuda-20261103.sif` (24.5 GiB), `micapipe-v1.0.0-beta.3-cuda-20261103.sif` (30.5 GiB) |
| SHA-256 (nocuda) | `99c296a7c92e548598be08006a8af66863321137f040ec34cd39e853b9c7a2cc` |
| SHA-256 (cuda) | `60f355a0473fbcd9e27cc6b5848033919e108ea52e801ddd49acb8b6636907a3` |

The source commit is `v1` plus the rebuild branch, so the images run slightly
ahead of `v1` until that branch is merged. Nothing else is on top of `v1`: the
MRIQC run-selection work is not in these images.

Versions are as recorded in [`containers/versions.json`](../../containers/versions.json):
FastSurfer 2.6.0 (FreeSurfer 8.2.0), FreeSurfer 8.2.0, MRtrix3 3.0.8,
ANTs 2.6.5, **FSL 6.0.7.23** (`fsl-eddy-cuda 2602.0`, Python FIX 0.10.0),
pyonlm, NiiVue 0.69.0.

Where the published images live, and how to rebuild them, is in
[`containers/README.md`](../../containers/README.md).

### This was an incremental rebuild

The `/opt` trees from the 20261028 build were reused. `/opt/fsl-6.0.7.8` and
`/opt/micapipe` were deleted first, so the FSL prefix was solved fresh from
`environments/locks/fsl-6.0.7.23-linux-64.txt` and the pipeline code tree was
re-copied from the new commit; every other tool is marker-guarded in
`containers/install_runtime.sh` and was re-verified rather than re-downloaded.
`containers/SOURCES.sha256` is unchanged since the 20261028 build and was
checked at the start of each build, and the pinned FastScrubber weight hashes
were re-checked inside every finished image.

One consequence worth stating: the OS layer came from the cached
`micapipe-v1-beta-os:20261022` Docker build, whose apt packages are not
snapshot-pinned, so bit-for-bit reproduction of these images is not claimed.

### Three builds were made; only the third is the candidate

Each build passed every static, FSL and smoke check. The pilots are what found
the two real defects, and both were in the pipeline rather than the packaging,
so the sequence is worth recording.

**First build (`2b65d3b`) -- FastSurfer could not use a multi-GPU host.** Its
CPU pilot passed in full. Its GPU pilot died inside FastSurfer: with all eight
GPUs visible, FastSurfer wraps its model in `torch.nn.DataParallel`
(`FastSurferCNN/inference.py`: `torch.cuda.device_count() > 1 and
device.index is None`), which sat in the first batch of coronal inference for
77 minutes and then failed with `NCCL Error 6: remote process exited or there
was a network error`. The GPUs were otherwise idle, so this was not
contention; re-running the same image with `CUDA_VISIBLE_DEVICES=0` got all
three predictions and the surface reconstruction, which isolated the cause.
`proc_surf` then made it worse by carrying on: it ran FastScrubber,
`clean_brain_mask.py` crashed on the aseg FastSurfer never wrote, and the
module only reported `INCOMPLETE` at the end.

**Second build (`be0fe0a`) -- the CUDA eddy had never run.** Both pilots
passed, including the GPU surfaces, so the FastSurfer fix held. But the eddy
command the GPU pilot recorded was `eddy_cpu`, and `proc_dwi` had taken 183
minutes. `eddy_cuda` is sound in the image -- it is on `PATH`, all its
libraries including `libcuda.so.1` resolve under `--nv`, and
`mrtrix3.fsl.eddy_binary(True)` returns it -- and `dwifslpreproc` does attempt
it first. It failed at `cudaGetDevice` with
`cudaErrorNoDevice: no CUDA-capable device is detected`, and MRtrix quietly
re-ran `eddy_cpu`.

The cause was micapipe's own GPU gating. `functions/utilities.sh` sets
`CUDA_VISIBLE_DEVICES=''` for the whole pipeline (added in `5cfc7d2` for a
BIC-SGE failure in `proc_dwi`), and only `restore_gpus` undoes it;
`01_proc-surf.sh` calls it for FastSurfer and FastScrubber, `02_proc-dwi.sh`
never did. So the CUDA eddy that PR #187 raised FSL to 6.0.7.23 for could not
run on any host, and nothing in the logs said why: `dwifslpreproc` writes its
reason to `eddy_cuda_failure_output.txt` inside a scratch directory the run
then deletes. Reproducing it needed `-nocleanup`.

**Third build (the candidate).** It carries the fixes for both:

- `proc_surf` keeps only the first GPU the caller offered, so FastSurfer never
  takes the DataParallel path. One GPU segments a scan in about a minute.
- Both FastSurfer passes and the mask cleanup stop the module where they fail,
  instead of running on with missing inputs.
- `proc_dwi` hands the GPUs back around the `dwifslpreproc` call, and only
  where the host really has one, so a node without a GPU keeps the quiet CPU
  path `5cfc7d2` was protecting.
- FastSurfer's version cache (`BUILD.info`) is written at build time, so its
  logs name the FastSurfer that built the image instead of printing four
  "detected dubious ownership" git errors and an empty version.
- `tests/beta_smoke.sh` now checks the GPU-eddy mechanism, not just that the
  binary exists: that the pipeline still hides the GPUs, that `restore_gpus`
  hands them back, and that `eddy_cuda` then stops reporting
  `cudaErrorNoDevice`. An image whose headline feature is unreachable should
  fail its own smoke test.

## Static and smoke checks

Nine runs on the candidate, with no `[FAIL]` lines
(`/export02/enning/micapipe-acceptance-20261103/runs` on the build host):

| check | nocuda | cuda (CPU) | cuda (`--nv`) |
| --- | --- | --- | --- |
| image acceptance (`acc_static.sh`) | PASS | PASS | PASS |
| FSL/eddy (`acc_fsl3.sh`) | PASS | PASS | PASS |
| in-image smoke (`tests/beta_smoke.sh`) | PASS | PASS | PASS |

`acc_static.sh` confirmed the `CONTAINER_VERSION`/`BUILDINFO.json` stamps and
source commit, the `-alt_fs` / `-alt_fs_tta` / `-surf_dir` help entries, the
pinned FastScrubber script and all five fold checkpoints, the per-variant
torch backends (`cpu` / `cu126` + `cu124`, with `torch.cuda.is_available()`
True under `--nv` and False without), and that `-alt_fs -freesurfer` and
`-alt_fs_tta` without `-alt_fs` are refused before any processing starts.

`acc_fsl3.sh` added the beta.3-specific checks:

- `/opt/fsl-6.0.7.23` is the only FSL prefix in the image, `$FSLDIR` points at
  it and `etc/fslversion` reads `6.0.7.23`.
- `topup`, `applytopup`, `eddy_cpu`, `eddy_cuda`, `eddy_quad`, `flirt` and
  `fslmaths` all resolve inside that prefix.
- `mrtrix3.fsl.eddy_binary(True)` returns `eddy_cuda` inside the image, so
  `dwifslpreproc` selects the CUDA eddy on its own -- MRtrix 3.0.8 has no
  `-use_cuda` flag, it tries CUDA first and falls back to the CPU binary.
- Python FIX is 0.10.0 and the legacy FIX wrapper points at the new prefix.
- Under `--nv`, `eddy_cuda` resolves its CUDA libraries and reaches FSL's own
  argument notice. Without `--nv` it stops at
  `libcuda.so.1: cannot open shared object file`.

The smoke test's own GPU-eddy check reported `GPU eddy is reachable on this
host` under `--nv`, and `No GPU here, so eddy_cpu is the expected path` in both
CPU modes.

Both images carry FastSurfer's version cache: `/opt/FastSurfer/BUILD.info`
reads `2.6.0+1dc0247` with the checkpoint hashes under it, and the pilots
logged `Version: 2.6.0+1dc0247` and zero `detected dubious ownership` errors,
where the first build printed four of them per FastSurfer call and recorded an
empty version.

## Pilots

Run on `/data/mica3/BIDS_CI` `sub-mri3T` `ses-01` with unmodified images
(nothing bound over `/opt`), 16 threads, in two parallel chains:
`-proc_structural -proc_surf -alt_fs -post_structural`, then `acc_verify.py`,
and for the GPU chain `-proc_dwi` with the reverse-phase b0 as `-dwi_rpe`.
The GPU chain runs with every GPU on the host visible, which is the condition
the first build died on.

Both chains finished clean, with no `[FAIL]` line in either:

| | CPU chain (`nocuda`) | GPU chain (`cuda --nv`) |
| --- | --- | --- |
| `acc_verify.py` | PASS | PASS |
| surfaces | 146 min | 109 min |
| `proc_dwi` | not run | 15 min |
| `post_structural` | 28 `.surf.gii`, 20 parcellations | 28 `.surf.gii`, 20 parcellations |
| white surface topology | Euler 2 both hemispheres | Euler 2 both hemispheres |
| mean thickness (lh / rh) | 2.169 / 2.160 mm | 2.166 / 2.136 mm |
| `-alt_fs` mask vs FastSurfer's | 1416.0 mL vs 1498.7 mL, Dice 0.97 | 1415.9 mL vs 1498.8 mL, Dice 0.97 |
| `detected dubious ownership` | 0 | 0 |

`acc_verify.py` checks properties rather than a reimplementation of
`clean_brain_mask.py`, because the first attempt at it reported a false
failure: it compared `mask.mgz` against the largest connected piece of the raw
prediction, which ignores the deep brain labels the script unions in first. It
now asserts that `mask.mgz` is a single connected piece with no interior holes,
that it keeps the whole largest piece of the prediction (0 voxels dropped, from
18 and 20 pieces respectively), that neither `aseg.auto.mgz` nor
`aseg.presurf.mgz` has any label outside it, that no deep brain label lies
outside it, that `norm.mgz` and `brainmask.mgz` are zero outside it, and that
the result differs from FastSurfer's own mask.

The GPU chain's recorded eddy command is the point of the exercise:

```
eddy_cuda --imain=eddy_in.nii --mask=eddy_mask.nii --acqp=eddy_config.txt \
  --index=eddy_indices.txt --bvecs=bvecs --bvals=bvals --topup=field \
  --data_is_shelled --slm=linear --repol --out=dwi_post_eddy --verbose
```

The same subject on the previous build recorded `eddy_cpu` and took 183
minutes for `proc_dwi`; on this one `eddy_cuda` runs and it takes 15. The
`eddy_cpu` fallback is not untested either -- it is what the previous build's
GPU chain ran end to end, and the smoke test asserts the no-GPU branch
directly.

## Not established

Unchanged from [beta.2](v1-beta-2.md): the converted FIX classifier failed the
compatibility gate, so legacy FIX with R and Octave stays the default; manual
visual sign-off of QC reports is still pending; scientific equivalence is not
claimed. Specific to this candidate: bit-for-bit rebuilds are not claimed (OS
apt packages are not pinned and the OS layer was cached), GPU and CPU runs are
not compared against each other for numerical equivalence, `-alt_fs_tta` is
built but not exercised, only one subject was piloted, and `-alt_fs` has still
not been re-evaluated against the 10-session MICs cohort that motivated its
hole-filling and deep-label restoration.

One gap found here and deliberately left alone: `micapipe` runs every module as
`$COMMAND | tee`, so the pipeline's exit status is `tee`'s. A module that
cannot even start -- an unreadable script, for instance -- still leaves
micapipe reporting success. Changing that touches every module and every qsub
path, so it does not belong in a container rebuild.
