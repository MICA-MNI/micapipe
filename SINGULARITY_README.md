# MICApipe Docker to Singularity Conversion

This repository contains the tools to convert MICApipe Docker images to Singularity SIF format for use on HPC systems without Docker support.

## Prerequisites

- Docker with the MICApipe image built
- Singularity/Apptainer installed
- 12TB+ storage space available
- No sudo rights required

## Quick Start

### Step 1: Build Docker Images

Build the base and main Docker images:

```bash
# Build base image
docker build -f Dockerfile.base -t ghcr.io/mica-mni/micapipe-base:latest .

# Build main image  
docker build -f Dockerfile.main -t ghcr.io/mica-mni/micapipe:latest .
```

### Step 2: Convert to Singularity SIF

Use the automated conversion script:

```bash
./build_singularity_via_tar.sh
```

This will:
1. Export Docker image to tar archive (~10-15 minutes)
2. Convert tar to SIF format (~15-20 minutes) 
3. Output: `/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif`

## Monitoring Progress

The conversion process takes 30-40 minutes total. Monitor progress:

```bash
# Watch tar export (Step 1)
watch -n 10 'ls -lh /host/cassio/export03/data/enning/micapipe_docker_*.tar 2>/dev/null || echo "Exporting..."'

# Watch SIF creation (Step 2)
watch -n 10 'ls -lh /host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif 2>/dev/null || echo "Building SIF..."'
```

## Files

### Essential Files
- `build_singularity_auto.sh` - Main Singularity build script
- `build_singularity_via_tar.sh` - Wrapper for tar-based conversion
- `Dockerfile.base` - Base image Dockerfile
- `Dockerfile.main` - Main image Dockerfile
- `Dockerfile` - Current working Dockerfile
- `micapipe_environment.yml` - Conda environment specification

### Core MICApipe Components
- `micapipe/` - Main Python package
- `functions/` - Processing scripts and utilities
- `parcellations/` - Brain atlases and parcellation data
- `surfaces/` - Surface processing templates
- `MICs60_T1-atlas/` - MNI template atlas
- `MNI152Volumes/` - Standard space templates

## Storage Requirements

- **Docker images**: ~200GB total
- **Tar export**: ~100GB (temporary)
- **Final SIF**: ~100GB
- **Cache/tmp**: ~50GB
- **Total needed**: ~350GB during build, ~150GB final

## Troubleshooting

### Build Hangs
If the build appears stuck:
```bash
# Check if processes are running
ps aux | grep docker
ps aux | grep singularity

# Check disk space
df -h /host/cassio/export03/data/enning/
```

### Memory Issues
The build uses significant memory:
- Docker build: Up to 32GB RAM
- Singularity conversion: Up to 16GB RAM

### Permission Issues
All operations use your user permissions:
- No sudo required
- Files created with your ownership
- Cache stored in your data directory

## Output

Success produces:
- SIF file: `/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif`
- Size: ~100GB
- Ready for HPC deployment

## Testing

Test the SIF file:
```bash
singularity run /host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif --help
```

## Environment Variables

The script uses these locations by default:
- `CACHE_ROOT=/host/cassio/export03/data/enning`
- `SINGULARITY_DIR=/host/cassio/export03/data/enning/singularity`
- `SINGULARITY_TMPDIR=/host/cassio/export03/data/enning/.singularity_tmp`
- `SINGULARITY_CACHEDIR=/host/cassio/export03/data/enning/.singularity_cache`

Override by setting environment variables before running the script.