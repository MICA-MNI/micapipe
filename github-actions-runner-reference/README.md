# GitHub Actions Runner Reference Files

This directory contains reference files from the GitHub Actions runner setup for micapipe.

## Files

- `Dockerfile`: The Dockerfile for building the GitHub Actions runner container with Singularity support
- `run_docker.sh`: Script to run the Docker container with appropriate volume mounts

## Original Location

These files were originally located at:
`/data_/mica1/03_projects/actions-runner/`

## Notes

- The Dockerfile sets up an Ubuntu 20.04 base with Singularity and GitHub Actions runner
- Includes Go 1.17.13 for Singularity compilation
- Pre-configures the runner for the MICA-MNI/micapipe repository
- Sets up volume mounts for CI data and Singularity cache directories

These files are for reference only and should not be modified in this location.