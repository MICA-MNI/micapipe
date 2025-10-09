# GitHub Actions Runner with Embedded Micapipe

This workflow creates a self-contained GitHub Actions runner with micapipe SIF file embedded inside.

## 🎯 **Purpose**
- **CI Mode**: Run micapipe tests in GitHub Actions workflows
- **Standalone Mode**: Use the runner independently with embedded micapipe
- **No Downloads**: SIF file is baked into the Docker image (no runtime downloads)

## 📋 **Complete Workflow**

### Step 1: Build Micapipe SIF (Simple v1 Method)
```bash
# Kill any stuck builds first
./kill_build.sh

# Build SIF using simple v1 method (reliable)
./build_singularity_v1.sh
```

**Expected output**: `/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif`

### Step 2: Build Actions Runner with Embedded SIF
```bash
# This will copy SIF and build the runner Docker image
./build_actions_runner.sh
```

**What this does**:
- Copies SIF to actions-runner directory
- Updates Dockerfile to embed the SIF
- Builds `micapipe-runner` Docker image
- Embeds ~100GB SIF inside the image

### Step 3: Deploy the Runner
```bash
cd /data_/mica1/03_projects/actions-runner
./run_docker.sh
```

## 🔍 **Verification**

Test that SIF is embedded:
```bash
docker run --rm micapipe-runner ls -la /data_/mica1/01_programs/micapipe-v0.2.0/
```

Test micapipe functionality:
```bash
docker run --rm micapipe-runner singularity run /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v1_beta.sif --help
```

## 📁 **File Locations**

**In the runner container**:
- SIF file: `/data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v1_beta.sif`
- Singularity cache: `/data/mica1/03_projects/enning/BIDS_CI/sing_cache`
- Singularity temp: `/data/mica1/03_projects/enning/BIDS_CI/sing_tmp`

**On the host system**:
- Actions runner: `/data_/mica1/03_projects/actions-runner/`
- SIF source: `/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif`

## 🚀 **Usage**

### CI Mode (GitHub Actions)
The runner will automatically pick up CI jobs from your GitHub repository.

### Standalone Mode
```bash
# Run micapipe directly
docker exec -it micapipe-runner singularity run /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v1_beta.sif [args]

# Interactive shell
docker exec -it micapipe-runner bash
```

## 🛠️ **Benefits**

✅ **Self-contained**: No external dependencies  
✅ **Fast startup**: No SIF download time  
✅ **Reliable**: Uses simple v1 build method  
✅ **Flexible**: Works for CI and standalone use  
✅ **Consistent**: Same environment every time  

## 📊 **Resource Requirements**

- **Docker image**: ~150GB (includes Ubuntu + SIF + tools)
- **Build time**: 45-90 minutes (SIF build + Docker build)
- **Runtime**: Minimal overhead, direct SIF execution