# MICApipe Dependencies Pre-Download Guide

This guide explains how to use the separate download script to pre-download FSL and FreeSurfer before building the Docker container.

## Quick Start

### 1. Download Dependencies First

```bash
# Download both FSL and FreeSurfer to ./downloads directory
./download_dependencies.sh

# Or specify custom directory
./download_dependencies.sh --download-dir /path/to/your/downloads
```

### 2. Build Container Using Pre-Downloaded Files

```bash
# Build using the downloaded files
./build_container.sh --downloads-dir ./downloads

# Or with custom path
./build_container.sh --downloads-dir /path/to/your/downloads
```

## Detailed Usage

### Download Script Options

```bash
./download_dependencies.sh [OPTIONS]

Options:
  --download-dir DIR    Directory to download files (default: ./downloads)
  --fsl-only           Download only FSL
  --freesurfer-only    Download only FreeSurfer
  --force              Force re-download even if files exist
  --verify             Only verify existing files, don't download
  --help, -h           Show help message
```

### Examples

```bash
# Download only FSL (faster, ~600MB)
./download_dependencies.sh --fsl-only

# Download only FreeSurfer (~3GB)
./download_dependencies.sh --freesurfer-only

# Force re-download all files
./download_dependencies.sh --force

# Verify existing downloads
./download_dependencies.sh --verify

# Use server cache directory
./download_dependencies.sh --download-dir /data_/mica1/01_programs/micapipe_downloads
```

### Build Script Integration

```bash
# Build with pre-downloaded dependencies
./build_container.sh --downloads-dir ./downloads

# Combine with other options
./build_container.sh --cuda --downloads-dir ./downloads

# Use both cache and downloads (cache takes priority)
./build_container.sh --cache-dir /cache --downloads-dir ./downloads

# Use custom temporary directory for build operations
./build_container.sh --custom-tmpdir /path/to/custom/temp --downloads-dir ./downloads
```

## Benefits of Pre-Downloading

1. **Reliability**: Download files once, build multiple times without network issues
2. **Speed**: No download time during Docker build process
3. **Bandwidth**: Reduces repeated downloads on servers
4. **Offline**: Build containers without internet access after initial download
5. **Debugging**: Separate download failures from build failures

## File Sizes and Requirements

- **FSL 6.0.2**: ~600MB compressed
- **FreeSurfer 7.4.1**: ~3GB compressed
- **Total**: ~3.6GB

Ensure sufficient disk space before downloading.

## Priority Order

The Docker build process checks for dependencies in this order:

1. **Cache directory** (if `--cache-dir` specified)
2. **Downloads directory** (if `--downloads-dir` specified)  
3. **Network download** (fallback)

## Troubleshooting

### Download Failures

```bash
# Check network connectivity
curl -I https://fsl.fmrib.ox.ac.uk
ping surfer.nmr.mgh.harvard.edu

# Retry with force download
./download_dependencies.sh --force

# Download individually
./download_dependencies.sh --fsl-only
./download_dependencies.sh --freesurfer-only
```

### Verification Issues

```bash
# Check file sizes
ls -lh ./downloads/

# Verify files manually
./download_dependencies.sh --verify

# Expected files:
# - fsl-6.0.2-centos6_64.tar.gz (≥500MB)
# - freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz (≥2GB)
```

### Build Issues

```bash
# Verify downloads directory is accessible
ls -la /path/to/downloads

# Check Docker can mount the directory
docker run --rm -v /path/to/downloads:/test ubuntu ls -la /test

# Build without pre-downloads as fallback
./build_container.sh  # Will download during build
```

## Integration with Existing Workflows

### For Development

```bash
# One-time setup
./download_dependencies.sh --download-dir ~/micapipe_deps

# Regular builds
./build_container.sh --downloads-dir ~/micapipe_deps
```

### For CI/CD

```bash
# Cache dependencies in CI
./download_dependencies.sh --download-dir $CI_CACHE_DIR/micapipe_deps

# Build in CI
./build_container.sh --downloads-dir $CI_CACHE_DIR/micapipe_deps
```

### For Server Deployment

```bash
# Setup shared downloads location
sudo mkdir -p /data_/mica1/01_programs/micapipe_downloads
sudo ./download_dependencies.sh --download-dir /data_/mica1/01_programs/micapipe_downloads

# All users can build using shared downloads
./build_container.sh --downloads-dir /data_/mica1/01_programs/micapipe_downloads
```