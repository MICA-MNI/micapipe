# CI/Testing Scripts and Tools

This folder contains scripts and documentation for Continuous Integration (CI) and testing workflows.

## Scripts

### Action Runner Setup
- `build_actions_runner.sh` - Build Docker image for GitHub Actions self-hosted runner
- `build_actions_runner_ci.sh` - CI-specific action runner build script  
- `prepare_action_runner.sh` - Prepare environment for action runner setup
- `build_test_runner.sh` - Build Docker image for testing environment

### Documentation
- `ACTIONS_RUNNER_README.md` - Detailed guide for setting up GitHub Actions runners
- `ACTION_RUNNER_BASE_SETUP.md` - Base setup instructions for action runners

## Usage

These scripts are for setting up CI infrastructure and testing environments. They are separate from the main micapipe Docker build process but support automated testing and deployment workflows.

### For GitHub Actions Self-Hosted Runners:
```bash
cd ci-testing
./prepare_action_runner.sh
./build_actions_runner.sh
```

### For Test Environment:
```bash
cd ci-testing  
./build_test_runner.sh
```

## Note

These scripts are organized separately from the main build process to avoid confusion with the core micapipe Docker build workflow.