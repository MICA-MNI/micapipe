# Micapipe Build Scripts

This directory contains scripts to help build and deploy micapipe containers.

## build_container.sh

A comprehensive script to build Docker and Singularity containers for micapipe, extracted from the CI workflow.

### Features
- Builds Docker images with customizable tags
- Optional Singularity image building
- **CUDA support toggle** - Enable/disable GPU acceleration
- Build progress tracking and timing
- Colored output for better readability
- Comprehensive error checking and validation
- Support for clean builds (no cache)

### Prerequisites
- Docker installed and running
- Singularity (optional, only if building .sif files)
- NVIDIA drivers and Docker GPU runtime (optional, only for CUDA builds)
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
./scripts/build_container.sh --no-cache --tag micapipe:latest
```

#### Enable CUDA support:
```bash
./scripts/build_container.sh --cuda --tag micapipe:v0.2.4-cuda
```

#### Build CUDA-enabled container with Singularity:
```bash
./scripts/build_container.sh --cuda --singularity --tag micapipe:v0.2.4-cuda
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
- `-c, --cuda`: Enable CUDA support in the container (default: disabled)
- `-h, --help`: Show help message

### Output
The script provides:
- Real-time build progress
- Build timing information
- Image size information
- CUDA support status
- GPU runtime detection (for CUDA builds)
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

## CUDA Support

The build script supports optional CUDA acceleration for GPU-enabled processing. When CUDA is enabled:

### What gets installed:
- CUDA runtime libraries (CUDA 11.8)
- cuDNN for deep learning
- GPU-enabled PyTorch and TensorFlow
- GPU-enabled FastSurfer environment
- GPU-enabled onnxruntime

### Requirements for CUDA builds:
- NVIDIA GPU with CUDA 11.8+ support
- NVIDIA drivers installed on host system
- Docker with NVIDIA runtime (nvidia-docker2)

### Usage:
```bash
# Build CUDA-enabled container
./scripts/build_container.sh --cuda --tag micapipe:cuda

# Run with GPU support
docker run --gpus all -it micapipe:cuda

# Run Singularity with GPU support  
singularity exec --nv micapipe.sif micapipe --help
```

### Default behavior:
- CUDA is **disabled by default** to preserve current CPU-only behavior
- Use `--cuda` flag to explicitly enable GPU support
- CPU-only builds work on any system without GPU requirements