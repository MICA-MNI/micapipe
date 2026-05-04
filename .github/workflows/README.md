# CI workflows

Three workflows, each on the right tier:

| Workflow | Runner | Triggers | Wall time |
|---|---|---|---|
| [`lint.yml`](lint.yml) | `ubuntu-latest` | every push, every PR | <2 min |
| [`build-main.yml`](build-main.yml) | `self-hosted` (BIC) | code changes (paths-ignore docs/base) | **~14 s** |
| [`build-base.yml`](build-base.yml) | `self-hosted` (BIC) | only when `Dockerfile.mamba-base` or `micapipe_environment*.yml` change | 45–60 min |

## Design intent

The comprehensive base image (`ghcr.io/mica-mni/micapipe-comprehensive-base`) is a stable artifact: it bundles every neuroimaging tool (FreeSurfer, FSL, ANTs, AFNI, MRtrix3, Workbench, dcm2niix, …) and only changes when one of those tool versions or `micapipe_environment.yml` changes. Rebuilding it is a 45–60 minute operation.

Day-to-day code iteration goes through the slim `Dockerfile.main` layer, which copies only `micapipe`, `functions/`, and `micapipe.py`. Rebuilding that layer is a 1-second operation against the cached base.

`build-main.yml` therefore **pulls** the base from the registry and **never rebuilds** it. The fast path is:

```
ensure base cached         0.04 s
build slim main image      1.0 s
docker run -h              6.0 s   (cold container start)
parser regression #1       1.7 s   (micapipe -sub → must error)
parser regression #2       3.0 s   (micapipe -sub -out → must error)
cleanup                    2.0 s
─────────────────────────
total                      ~14 s
```

Singularity SIF conversion (10–15 min for the 50GB image) is **not** in the fast path. SIF + `sample_test.sh` runs only when a PR carries the `ci:build` label, or via `workflow_dispatch` with `run_sample_test=true`.

## When does each workflow fire?

| You change… | `lint.yml` | `build-main.yml` | `build-base.yml` |
|---|---|---|---|
| `functions/*.sh`, `micapipe`, `*.py` | ✓ | ✓ | — |
| `Dockerfile.main` | ✓ | ✓ | — |
| `Dockerfile.mamba-base` | ✓ | — | ✓ |
| `micapipe_environment.yml` | ✓ | — | ✓ |
| Docs (`*.md`, `docs/`) | ✓ | — | — |
| Workflow files | ✓ | — (paths-ignore) | — (paths-ignore unless its own .yml) |

## Labels

- **`ci:build`** — on a PR, run `build-main.yml`'s full path (build SIF, run `sample_test.sh`). Multi-hour. Remove when done.
- **`ci:build-base`** — on a PR that touches `Dockerfile.mamba-base`, run the 45-60 min base rebuild. Otherwise, base PRs only run lint until the label is set.

## Server setup

The self-hosted runner lives on `bb-compxg-01` under `/data_/mica1/03_projects/actions-runner/`. See [`ci-testing/README.md`](../../ci-testing/README.md) for runner build/setup. If `build-main` or `build-base` jobs sit in `pending` indefinitely, the runner agent is probably offline — start it with `./run.sh` or `./svc.sh start` from that directory.

## Local equivalent of CI's lint

Run the same checks locally before pushing:

```bash
# Bash syntax across every .sh + the wrapper
find . -path ./.git -prune -o -path './docs/build*' -prune -o -name '*.sh' -print \
  | xargs -I{} bash -n {}
bash -n micapipe functions/micapipe_anonymize

# Python syntax
find . -path ./.git -prune -o -path '*__pycache__*' -prune -o -name '*.py' -print \
  | xargs -I{} python3 -m py_compile {}

# YAML
find . -name '*.yml' -o -name '*.yaml' \
  | xargs -I{} python3 -c "import yaml; yaml.safe_load(open('{}'))"

# Snakefile parse
python3 -c "import compileall; compileall.compile_file('micapipe_snakebids/workflow/Snakefile')"
```
