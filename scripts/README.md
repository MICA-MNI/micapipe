# Micapipe Build Scripts

This directory contains scripts to help build and deploy micapipe containers.

## build_container.sh

A comprehensive script to build Docker and Singularity containers for micapipe, extracted from the CI workflow.

### Features
- Builds Docker images with customizable tags
- Optional Singularity image building
- Build progress tracking and timing
- Colored output for better readability
- Comprehensive error checking and validation
- Support for clean builds (no cache)

### Prerequisites
- Docker installed and running
- Singularity (optional, only if building .sif files)
- Sufficient disk space for builds

### Usage

#### Basic build:
```bash
./scripts/build_container.sh
```

#### Build with custom tag:
```bash
./scripts/build_container.sh --tag micapipe:v0.2.4
```

#### Build both Docker and Singularity images:
```bash
./scripts/build_container.sh --tag micapipe:v0.2.4 --singularity
```

#### Clean build (no cache):
```bash
./scripts/build_container.sh --no-cache
```

#### Build with custom Singularity directory:
```bash
./scripts/build_container.sh --singularity --singularity-dir /custom/path
```

### Options
- `-t, --tag TAG`: Docker tag to use (default: micapipe:latest)
- `-s, --singularity`: Also build Singularity image
- `-d, --singularity-dir DIR`: Directory for Singularity build (default: /tmp/singularity_build)
- `-n, --no-cache`: Build without using Docker cache
- `-h, --help`: Show help message

### Output
The script provides:
- Real-time build progress
- Build timing information
- Image size information
- Usage examples after successful build
- Detailed error messages if build fails

### Example Output
```
[INFO] Checking prerequisites...
[SUCCESS] Prerequisites check passed
[INFO] Building Docker image with tag: micapipe:v0.2.4
[SUCCESS] Docker image built successfully in 1234 seconds
[SUCCESS] Singularity image built successfully in 567 seconds
[SUCCESS] Build completed successfully!
```

## Based on CI Workflow
This script is based on the steps from `.github/workflows/ci_test.yml` and includes:
1. Docker image building
2. Singularity conversion (optional)
3. Proper error handling and validation
4. Progress tracking and reporting