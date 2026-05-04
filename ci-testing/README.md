# CI/Testing Scripts and Tools

Scripts that support our CI/testing infrastructure. These are separate
from the main micapipe Docker build process; they exist to bootstrap
self-hosted GitHub Actions runners on the BIC servers.

## Scripts

- `build_actions_runner.sh` — Build Docker image for the self-hosted GitHub
  Actions runner with the micapipe SIF embedded.
- `build_actions_runner_ci.sh` — CI-specific runner build (used in the
  workflow).
- `prepare_action_runner.sh` — Prepare the host environment before
  building the runner image.
- `build_test_runner.sh` — Build Docker image for ad-hoc testing.

## Usage

For self-hosted GitHub Actions runners:

```bash
cd ci-testing
./prepare_action_runner.sh
./build_actions_runner.sh
```

For an ad-hoc test environment:

```bash
cd ci-testing
./build_test_runner.sh
```

These scripts assume MICA workstation paths
(`/data_/mica1/03_projects/actions-runner`,
`/host/cassio/export03/data/enning/singularity/`); adjust per host.
