# Changelog

All notable changes to micapipe will be documented in this file.

## [Unreleased] — comprehensive cleanup branch

### CI / build infrastructure

- **Two-stage Docker build** kept; the slim `Dockerfile.main` now copies
  only the fast-changing code layer (`micapipe`, `functions/`,
  `micapipe.py`). Static assets (parcellations, surfaces, MNI152,
  fsl_conf, R_config, fix_settings) live in the comprehensive base
  image. Iterative main rebuild on a warm runner: ~1 second.
- **Three workflow split**: `lint.yml` runs on every push/PR on
  ubuntu-latest (<2 min), `build-main.yml` is the fast self-hosted
  code-iteration path (~14 s end-to-end with smoke tests), and
  `build-base.yml` only triggers on changes to `Dockerfile.mamba-base`
  or the conda env files. SIF + `sample_test.sh` is opt-in via the
  `ci:build` PR label.
- **`build_singularity.sh` auto-picks transport** between
  `docker-daemon://` (fast) and `docker-archive://` (works when Docker
  data root is tight). Override with
  `MICAPIPE_SIF_TRANSPORT={auto,docker-daemon,docker-archive}`.
- **Server-side build scripts parameterized** via `MICAPIPE_*` env vars
  and CLI flags — no more server-hardcoded paths.

### Pipeline

- **#156 — LAMAReg replaces regSynth** as the default cross-modality
  registration. CLI flags `-regSynth`, `-regAffine`, `-reg_nonlinear`,
  `-microstructural_reg` removed. Two-stage robust registration emits
  primary + secondary warpfields and DICE-score CSV under `xfm/`.
  `02_proc-{dwi,func,flair}.sh`, `03_MPC{,-SWM,SC}.sh`, snakebids
  rules, and tests all updated.
- **#153 — MNI atlas updated** to MNI NLIN SYM 09a 1mm (already merged
  via #164).
- **#162 — Argument parser hardened**:
  - Wrapper now uses `while [ "$#" -gt 0 ]; do case "$1"` instead of
    the old `for arg in "$@"` (which iterated a frozen snapshot and
    silently swallowed malformed flag/value pairs).
  - New `require_value` helper invoked by critical value-taking flags
    (`-sub`, `-out`, `-bids`, `-ses`, `-tmpDir`, `-fs_licence`, `-T1`,
    `-threads`). Catches `micapipe -sub -out /foo` (missing value)
    with a clear error.
- **#162 — FLAIR session-isolated tmp**: `02_proc-flair.sh` tmp path
  now embeds nanosecond timestamp + PID alongside `RANDOM` and
  `idBIDS`, so two sequential ses-01/ses-02 invocations cannot collide
  on a shared `tmpDir`. Trap quoted to handle paths with spaces.
- **`cleanup()` hardened**: refuses to delete suspicious paths (empty,
  `/`, no separator) so a mis-set `tmpDir` cannot wipe the working
  tree. `bids_variables_unset` now also clears `idBIDS`, `ses`, and
  `bids_flair`.

### Container

- **`functions/init.sh` Docker detection rewritten**. The previous
  detection probed `/opt/freesurfer-7.4.1/SetUpFreeSurfer.sh` (which
  doesn't exist in the actual built image) and fell through to MICA
  cluster paths. Now probes `/opt/micapipe && /opt/freesurfer` and
  points each tool var at the actual install path:
  - FreeSurfer: `/opt/freesurfer` (unversioned)
  - FSL: `/opt/fsl-6.0.2`
  - MRtrix3: `/opt/miniconda-latest/envs/micapipe` (conda env, not a
    standalone /opt directory)
  - AFNI: `/opt/afni`
  - ANTs: `/opt/ants-2.3.4/bin`
  - c3d: `/opt/c3d-1.0.0-Linux-x86_64/bin`
  - FIX: `/opt/fix`
  - FastSurfer: `/opt/FastSurfer`
- **Wrapper auto-sources `init.sh` inside the container** (detected
  by `/opt/micapipe + /opt/freesurfer`). Previously only `-mica`,
  `-qsub`, and `-qall` modes sourced init.sh, so default invocations
  inside the container had no `FREESURFER_HOME`/`FSLDIR`/conda env.
- **Base image gaps filled** (PR #174 + #175 merged into this branch):
  - `dcm2niix` (1.0.20240202) — required for DICOM → NIfTI conversion.
  - `connectome-workbench` (`wb_command`) — required by every micapipe
    module.

### Cleanup

- 10 root-level debugging-journal markdown files removed.
- Stale `regAffine` default removed from `snakebids.yml`.
- `tests/sample_test.sh` — drop the now-removed `-regSynth` and
  `-microstructural_reg` flags (they would error against the new
  strict parser).
- Obsolete `.github/workflows/blank.yml` template removed.

---

## [v0.2.3] — previous release

See [the GitHub releases page](https://github.com/MICA-MNI/micapipe/releases)
for entries before v1.0.0-beta.
