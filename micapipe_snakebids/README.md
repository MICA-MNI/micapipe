# micapipe_snakebids

Runs micapipe as a [Snakemake](https://snakemake.github.io) workflow through a
[SnakeBIDS](https://github.com/khanlab/snakebids) app: one job per micapipe
module and session, sessions in parallel, and a job counted as done only when
micapipe itself recorded the module as `COMPLETED`.

## Running it

The micapipe image already carries Snakemake 8.30 and SnakeBIDS 0.13 in
`/opt/environments/workflow`, so the simplest way is to run the app inside the
image, where `micapipe` is on `PATH`:

```bash
singularity exec --cleanenv --containall \
  --bind /path/to/bids:/bids:ro,/path/to/derivatives:/out,/path/to/scratch:/tmp,/path/to/license.txt:/opt/licence.txt:ro \
  --pwd /out micapipe-v1.0.0-beta.3-nocuda-20261103.sif \
  /opt/environments/workflow/bin/python /opt/micapipe/micapipe_snakebids/run.py \
  /bids /out participant \
  --participant-label HC001 --modules structural surface post_structural \
  --fs-licence /opt/licence.txt --tmp-dir /tmp --threads 8 --cores 16
```

Bind real disk over `/tmp`: under `--containall` the container's own `/tmp` is a
64 MB tmpfs, and micapipe's tools (MRtrix among them) fill it within minutes.
The app refuses to start inside such a `/tmp` without `--tmp-dir`.

Or run it from a host Python that has Snakemake and SnakeBIDS, and let every
micapipe call go through the image:

```bash
python micapipe_snakebids/run.py /path/to/bids /path/to/derivatives participant \
  --container-image /path/to/micapipe-v1.0.0-beta.3-nocuda-20261103.sif \
  --fs-licence /path/to/license.txt --tmp-dir /path/to/scratch \
  --modules dwi sc --threads 8 --cores 16
```

With `--container-image` the scratch (`--tmp-dir`, default `<derivatives>/tmp`)
is bound as the container's `/tmp`.

The derivatives directory is micapipe's `-out`: micapipe writes
`micapipe_v0.2.0/` there, and the app adds `logs/micapipe_snakebids/` (one log
and one benchmark per job), `config/` and `.snakemake/`. Add `-n -p` for a dry
run that prints every micapipe call it would make. Any other Snakemake option
(`--cores`, `--keep-going`, `--profile`, ...) is passed through.

## Modules

`--modules` takes any of

| module | micapipe | waits for |
| --- | --- | --- |
| `structural` | `-proc_structural` | |
| `surface` | `-proc_surf` | structural |
| `post_structural` | `-post_structural` | surface |
| `flair` | `-proc_flair` | post_structural |
| `dwi` | `-proc_dwi` | post_structural |
| `sc` | `-SC` | dwi |
| `mpc` | `-MPC` | post_structural |
| `mpc_swm` | `-MPC_SWM` | post_structural |
| `gd` | `-GD` | post_structural |
| `asl` | `-proc_asl` | post_structural |
| `func` | `-proc_func`, one job per scan | post_structural |
| `swm` | `-SWM` | post_structural |
| `qc_subj` | `-QC_subj` | structural |

A requested module brings in the modules it waits for, so `--modules sc` alone
runs structural, surface, post_structural, dwi and SC. A module whose record
already says `COMPLETED` is not run again.

Within one session the modules run one at a time, in the order of the table:
micapipe rewrites session-level files as it goes, and two modules writing them
at once can corrupt them. Different sessions run in parallel, up to `--cores`.

## How a job is judged

Each job runs one micapipe call and then requires micapipe's record for that
module, `micapipe_v0.2.0/sub-X/ses-Y/QC/sub-X_ses-Y_module-<name>.json`, to say
`"Status": "COMPLETED"` for that subject, session and module. micapipe runs every
module as `$COMMAND | tee`, so it exits 0 even when the module failed; the
record is what tells. A job whose record is missing or `INCOMPLETE` fails, with
the reason at the end of its log. `qc_subj` is judged by its structural report.

The record names follow micapipe: `proc_surf-fastsurfer` (or `-freesurfer`),
`proc_dwi[_acq-X]`, `SC-<tracts>[_acq-X]`, `MPC-<acq>`, `MPC-SWM-<acq>`,
`proc_func-desc-se_<scan>` (`desc-me_` for multi-echo).

## Options

Every micapipe option has a flag; `run.py --help` lists them all. A flag left
out is not passed, so micapipe keeps its own default.

| flag | micapipe |
| --- | --- |
| `--threads N` | `-threads` (capped by `--cores`) |
| `--fs-licence FILE` | `-fs_licence` (needed whenever structural, surface or post_structural run) |
| `--tmp-dir DIR`, `--nocleanup` | `-tmpDir`, `-nocleanup` |
| `--t1w-str`, `--uni`, `--mf` | `-T1wStr`, `-uni`, `-mf` |
| `--freesurfer`, `--alt-fs`, `--alt-fs-tta`, `--surf-dir`, `--surf-t1` | `-freesurfer`, `-alt_fs`, `-alt_fs_tta`, `-surf_dir`, `-T1` |
| `--atlas A B ...` | `-atlas A,B` |
| `--flair-scan-str` | `-flairScanStr` |
| `--dwi-main`, `--dwi-rpe`, `--dwi-processed`, `--rpe-all`, `--b0thr`, `--no-bvalue-scaling`, `--dwi-upsample`, `--dwi-acq` | the `-proc_dwi` options |
| `--tracts`, `--keep-tck`, `--auto-tract`, `--weighted-sc`, `--tck` | the `-SC` options |
| `--microstructural-img`, `--mpc-acq` | `-MPC` and `-MPC_SWM` options |
| `--asl-scan-str`, `--m0-scan-str` | `-aslScanStr`, `-m0ScanStr` |
| `--func-scans`, `--func-pe`, `--func-rpe`, `--main-scan-run`, `--phase-reversal-run`, `--topup-config`, `--icafix-training`, `--ses-anat`, `--nsr`, `--gsr`, `--no-fix`, `--drop-tr`, `--no-fc` | the `-proc_func` options |
| `--session-label` | only these sessions |
| `--container-image`, `--container-exe`, `--container-binds`, `--container-gpu` | run micapipe through singularity/apptainer (`--nv` with `--container-gpu`) |

Images of a session can be named instead of spelled out: `--dwi-main
acq-b2000-91_dir-AP_dwi`, `--func-pe task-rest_acq-APse_bold`,
`--microstructural-img T1map` and `--surf-t1 T1w` each become that session's
`sub-X_ses-Y_<name>.nii.gz` in the right folder. A value with a `/` is used as a
path.

`--func-scans` names the scans to process, one job each
(`task-rest_acq-AP_bold`; join the echoes of a multi-echo scan with commas).
Without it, each session processes the rest scan micapipe would pick itself
(`task-rest_acq-AP_*bold`), and the app refuses to guess between several.

`--ses-anat 01` with `--modules func` processes a later session's func against
session 01's anatomy, and waits for session 01's structural chain.

Options can also be set in a YAML file passed with `--configfile`, under the
same names as the flags (`fs_licence`, `dwi_acq`, ...). Booleans are `true` /
`false`.

## Sessions

The app processes every `sub-*/ses-*` whose `anat/` holds a `*_T1w` or `*_UNIT1`
image, whatever other entities it has (`run-`, `acq-`); `--participant-label`,
`--exclude-participant-label` and `--session-label` narrow that. micapipe picks
the image itself (see `--t1w-str`). Datasets without `ses-` directories are not
supported yet.

## Not handled here

- micapipe's group QC (`-QC`) has no rule.
- `micapipe_v0.2.0/micapipe_processed_sub.csv` is rewritten by every module of
  every session, without a lock, so sessions running at once can drop each
  other's lines from it. The per-module QC records are not affected, and they
  are what this workflow reads.
- A module re-run because something before it changed resumes from whatever
  micapipe left on disk; micapipe has no notion of upstream changes.

## Tests

```bash
python micapipe_snakebids/tests/test_micapipe_workflow.py
```

covers module ordering, record names, flag rendering (including that nothing is
ever rendered that micapipe's parser would loop on), session discovery and the
record check.
