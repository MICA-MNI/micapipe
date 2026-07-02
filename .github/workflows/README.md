# Container images & CI/CD

How micapipe's Docker/Singularity images are built and how the GitHub Actions
pipelines behave. If you opened a PR and a check was **skipped** (a ⊘ icon
rather than a ✓ or ✗), that is almost always **by design** — see
[Fast path vs. full path](#fast-path-vs-full-path) below.

## Two-stage image architecture

The image is split in two so day-to-day code changes rebuild in seconds
instead of hours.

| Stage | Dockerfile | Image | Contents | Build time |
|-------|------------|-------|----------|------------|
| **Base** (stable) | `Dockerfile.mamba-base` | `ghcr.io/mica-mni/micapipe-comprehensive-base:latest` | Every neuroimaging dependency: conda/mamba, FSL, FreeSurfer, AFNI, ANTs, MRtrix3, workbench, dcm2niix, FastSurfer (+ `fastsurfer_cpu` env), R, and the legacy X/Mesa libs. ~52 GB. | slow (45–240 min) |
| **Main** (fast) | `Dockerfile.main` | `ghcr.io/mica-mni/micapipe:latest` | `FROM` the base, then `COPY`s only the micapipe code (`micapipe`, `micapipe.py`, `functions/`). | seconds–minutes |

The base is treated as a **stable artifact**: it is built rarely, pushed to
`ghcr.io`, and pulled (never rebuilt) by the fast workflow. The main image is
just a thin code layer on top.

Supporting scripts:
- `prepare_build_context.sh` — stages the build context (`--target base|main`).
- `build_singularity.sh` — converts the Docker image to a Singularity `.sif`.
- `tests/sample_test.sh` — end-to-end pipeline run against a sample subject.

## Workflows

| Workflow | File | Runner | When it runs |
|----------|------|--------|--------------|
| **Lint** | `lint.yml` | GitHub-hosted (`ubuntu-latest`) | Every push to `master`/`v1`/`enning/**` and every PR. Fast. |
| **Build base image** | `build-base.yml` | self-hosted | Push to `master` that changes `Dockerfile.mamba-base` / `micapipe_environment*.yml`; PRs touching those paths **only with the `ci:build-base` label**; or manual dispatch. |
| **Build main image** | `build-main.yml` | self-hosted | Push to `master`/`v1`/`enning/**` and PRs (ignoring `*.md`, `docs/**`, and the base-image inputs). |

### Lint (`lint.yml`)
Pure static checks, no containers, on a GitHub-hosted runner. Five jobs:
Bash syntax (`bash -n`), Python syntax (`py_compile`), YAML lint
(`yaml.safe_load`), Dockerfile lint (`hadolint`), and Snakemake parse
(Snakefile parse-only, no execution). This is the gate that protects `v1`.

### Build main image (`build-main.yml`) — the one you'll see most

#### Fast path (default)
On every qualifying push/PR:
1. Ensure the base image is present (pull from `ghcr.io`, **never build**).
2. Build the slim main image from `Dockerfile.main`.
3. Smoke-test: `docker run … micapipe -h`, plus argument-parser regression
   checks.

Total wall time is ~10 s on a warm runner (base already cached).

#### Full path (opt-in)
Three additional steps — **Build Singularity SIF**, **Point sample_test at the
SIF**, and **Run sample_test** — are gated and **skipped by default** (they
show as ⊘). `sample_test` is a multi-hour end-to-end run, so it only fires
when explicitly requested:

- add the **`ci:build`** label to a PR, **or**
- trigger **Run workflow** manually (`workflow_dispatch`) with
  `run_sample_test = true`.

The final *Cleanup CI artifacts* step always runs (`if: always()`), dropping
the run-specific image tag and SIF while keeping `:latest`.

> **A skipped (⊘) SIF/sample_test step is not a failure.** It means the full
> path wasn't requested for that run. The run is green as long as lint, the
> slim build, and the smoke test pass.

### Build base image (`build-base.yml`)
Rebuilds the comprehensive base (the multi-hour tool install). On PRs it is
guarded behind the **`ci:build-base`** label so a 45–60 min rebuild doesn't
fire on every push to a base-image PR. On a `push` event (or manual dispatch
with `push = true`) it logs in to `ghcr.io` and pushes the new base.

## CI labels

| Label (on a PR) | Effect |
|-----------------|--------|
| `ci:build` | `build-main` runs the **full** path (SIF + `sample_test`). |
| `ci:build-base` | `build-base` actually **rebuilds** the base image for that PR. |

## Self-hosted runner

`build-base` and `build-main` run on a **self-hosted** runner because they need
a Docker daemon, Singularity, the ~52 GB base image cached locally, and large
scratch space — none of which a GitHub-hosted runner provides.

- **Host:** a BIC server (`bb-compxg-01`), reached via the `bic-login` jump
  host. It has Docker, Singularity, and the `/export03` scratch volume the
  workflows reference.
- **Install dir:** `/export03/data/enning/actions-runner`
- **Label:** default `self-hosted` (matches `runs-on: self-hosted`).

### Operating the runner

```bash
ssh bb-compxg-01                       # via the bic-login ProxyJump
cd /export03/data/enning/actions-runner

# Is it alive?
pgrep -af Runner.Listener              # a running PID = online
tail -f runner.log                     # "Listening for Jobs" / job results

# Start it (no sudo): survives logout, NOT reboot
nohup ./run.sh > runner.log 2>&1 &

# Start it as a service (needs sudo): survives reboot
sudo ./svc.sh install <user> && sudo ./svc.sh start && sudo ./svc.sh status
```

If jobs sit **queued** for a long time, the runner is offline — restart it as
above. Every self-hosted job will queue (not fail) until a runner picks it up.

### CI input/output paths

All `sample_test` I/O lives under one **shared, non-personal** root on the
runner's scratch volume (no more `.../enning/...` home paths):

```
/export03/data/action-runner/           # MICAPIPE_CI_ROOT
├── micapipe_v1_beta.sif                 # default SIF (CI overrides with the fresh build)
├── rawdata/                             # input BIDS dataset
├── license_fc.txt                       # FreeSurfer licence
├── tmp/                                 # scratch
├── sif/micapipe_ci_<run_id>.sif         # run-specific SIF (deleted on cleanup)
└── output/                              # subject outputs, one timestamped dir per run
    ├── singularity_0.2.3_<timestamp>/
    └── singularity_0.2.3_<timestamp>_freesurfer/
```

Every path is overridable via env var — `MICAPIPE_CI_ROOT` moves the whole
tree, or set `MICAPIPE_CI_BIDS` / `MICAPIPE_CI_LICENSE` / `MICAPIPE_CI_SIF` /
`MICAPIPE_CI_TMP` / `MICAPIPE_CI_OUTDIR` individually (see the top of
`tests/sample_test.sh`). Outputs are timestamped, so runs never overwrite each
other.

`rawdata`, `license_fc.txt` and `micapipe_v1_beta.sif` are **symlinks** to the
existing shared copies (`/data/mica3/BIDS_CI/...`, `/data_/mica1/...`), so no
data is duplicated. For convenience there is also a shortcut next to the input
dataset:

```
/data/mica3/BIDS_CI/ci_output -> /export03/data/action-runner/output
```

**Access / permissions.** The tree is group `mica` with setgid **and** default
ACLs (`g:mica:rwx`), so every artifact the runner writes is readable, editable
and deletable by anyone in the group — not just whoever the runner ran as. To
grant a new person access, add them to the `mica` group; no per-run `chmod` is
needed. If you ever recreate the tree, re-apply the ACLs:

```bash
setfacl -R -m g:mica:rwx    /export03/data/action-runner
setfacl -R -d -m g:mica:rwx /export03/data/action-runner   # default for new files
```

### Registering a new runner
Requires repo **Admin** (Maintainer is not enough). Generate a token at
*Settings → Actions → Runners → New self-hosted runner*, then on the host:

```bash
./config.sh --url https://github.com/MICA-MNI/micapipe --token <TOKEN> \
            --name bb-compxg-01 --work _work --unattended --replace
```
