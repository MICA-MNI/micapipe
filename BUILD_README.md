# MICApipe Container Build Scripts

## Quick Start

### Comprehensive Build (Recommended)
```bash
./server_build_test.sh
```
- Complete build, test, and Singularity conversion
- Creates: `micapipe:latest` (Docker) + `/opt/micapipe/micapipe_v1-beta.sif` (Singularity)
- Comprehensive logging and validation

### Simple Singularity Build (Faster)
```bash
./build_container.sh
```
- Direct Singularity .sif build (Docker as intermediate step)
- Creates: `/opt/micapipe/micapipe_v1-beta.sif` (Singularity only)
- Automatically removes intermediate Docker container

## Output Files

**Primary Output:** `/opt/micapipe/micapipe_v1-beta.sif`
- Singularity container with all updated neuroimaging tools
- MRtrix 3.0.7, FreeSurfer 7.4.1, FastSurfer v2.4.2, etc.
- Custom location: `MICAPIPE_SIF_PATH=/custom/path ./build_container.sh`

**Docker Container:** `micapipe:latest` (only from server_build_test.sh)
- Available only when using comprehensive build script
- Intermediate containers are automatically cleaned up

## Build Options

### Simple Build Script
```bash
./build_container.sh --help      # Show all options
./build_container.sh --cuda      # Build with CUDA support  
./build_container.sh --no-cache  # Build without cache
./build_container.sh --output /custom/path  # Custom output directory
```

### Comprehensive Build Script
```bash
./server_build_test.sh help            # Show all options
./server_build_test.sh fsl-only        # Test FSL section only
./server_build_test.sh singularity-only # Convert existing Docker to .sif
```

## Requirements

- Docker installed and running
- At least 10GB free disk space
- At least 4GB RAM
- Internet connection
- Singularity (optional, for .sif creation)

## Troubleshooting

Check build logs in `./build_logs/` directory for detailed error information.

For comprehensive deployment guide, see: [SERVER_DEPLOYMENT.md](SERVER_DEPLOYMENT.md)