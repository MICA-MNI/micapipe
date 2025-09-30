# MICApipe Container Build Scripts

## Quick Start

### Comprehensive Build (Recommended)
```bash
./server_build_test.sh
```
- Complete build, test, and Singularity conversion
- Creates: `micapipe:latest` (Docker) + `/opt/micapipe/micapipe_v1-beta.sif` (Singularity)
- Comprehensive logging and validation

### Simple Build (Faster)
```bash
./build_container.sh
```
- Quick Docker container build
- Creates: `micapipe:latest` (Docker only)
- Minimal logging

## Output Files

**Docker Container:** `micapipe:latest`
- Contains all updated neuroimaging tools
- MRtrix 3.0.7, FreeSurfer 7.4.1, FastSurfer v2.4.2, etc.

**Singularity Container:** `/opt/micapipe/micapipe_v1-beta.sif`
- Automatically created by `server_build_test.sh`
- Custom location: `MICAPIPE_SIF_PATH=/custom/path ./server_build_test.sh`

## Build Options

```bash
./server_build_test.sh help    # Show all options
./server_build_test.sh fsl-only # Test FSL section only
./fix_fsl_build.sh            # Fix FSL installation issues
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