# MICApipe Container Build Guide

Simple container building for MICApipe with no sudo privileges required.

## Quick Start

### Basic Build (CPU-only)
```bash
./build_container.sh
```

### Build with CUDA Support
```bash
./build_container.sh --cuda
```

### Build without Cache (clean build)
```bash
./build_container.sh --no-cache
```

## Requirements

- **Docker**: Access without sudo (user in docker group)
- **Singularity**: For .sif container creation
- **Memory**: 8GB+ available RAM recommended
- **Disk**: 20GB+ free space
- **Network**: Stable connection for downloads (3-5GB total)

## Build Process

1. **System Checks**: Verifies Docker/Singularity access and resources
2. **Cleanup**: Removes old containers and images
3. **Docker Build**: Creates intermediate Docker image (30-60 minutes)
4. **Container Test**: Validates the built container
5. **Singularity Convert**: Creates .sif file from Docker image
6. **Cleanup**: Removes intermediate Docker image

## Output

- **Default Location**: `/data_/mica1/01_programs/singularity/micapipe_v1-beta.sif`
- **Fallback Location**: `$HOME/micapipe_v1-beta.sif`
- **Build Logs**: `build_logs/container_build_YYYYMMDD_HHMMSS.log`

## Usage

### Test Container
```bash
singularity exec /path/to/micapipe_v1-beta.sif micapipe --help
```

### Run with Data
```bash
singularity exec --bind /path/to/data:/data \
  /path/to/micapipe_v1-beta.sif \
  micapipe -sub 001 -out /data/derivatives -bids /data
```

## Troubleshooting

### Docker Permission Issues
```bash
# Ask admin to add you to docker group
sudo usermod -aG docker $USER
newgrp docker  # or logout/login
```

### Memory Issues (exit code 137)
- Close other applications
- Ask admin about increasing swap space
- Build during off-peak hours

### Network Issues
- Check internet connectivity
- Try building at different times
- Use `--no-cache` for clean retry

### Low Disk Space
```bash
# Clean Docker system
docker system prune -f

# Use custom output path
./build_container.sh --singularity-path /path/with/space
```

## Options

| Option | Description |
|--------|-------------|
| `--cuda` | Enable CUDA support for GPU acceleration |
| `--no-cache` | Build without Docker cache (clean build) |
| `--singularity-path PATH` | Custom output directory for .sif file |
| `--help` | Show usage information |

## Examples

```bash
# Standard CPU build
./build_container.sh

# GPU-enabled build
./build_container.sh --cuda

# Clean build without cache
./build_container.sh --no-cache

# Custom output location
./build_container.sh --singularity-path /custom/path

# Combined options
./build_container.sh --cuda --no-cache --singularity-path /opt/containers
```

## Support

If build issues persist:
1. Check build logs in `build_logs/`
2. Contact system administrator for Docker/Singularity setup
3. Consider building on a system with more resources