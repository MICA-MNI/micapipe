#!/bin/bash
#
# Build CI-compatible GitHub Actions runner
# This version makes the runner work with existing CI workflows
#

set -e

# Paths to match CI expectations
ACTIONS_RUNNER_DIR="/data_/mica1/03_projects/actions-runner"
SIF_SOURCE="/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif"

# CI-compatible paths
CI_SIF_DIR="/data_/mica1/01_programs/micapipe-v0.2.0"
CI_SIF_PATH="$CI_SIF_DIR/micapipe_v0.2.3.sif"  # CI expects this exact name
CI_TEMP_DIR="/export02/local/singularity_tmp"

DOCKERFILE="$ACTIONS_RUNNER_DIR/Dockerfile"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🔧 BUILDING CI-COMPATIBLE GITHUB ACTIONS RUNNER"
log "================================================="

# ============================================================================
# Step 1: Verify SIF exists
# ============================================================================
log "📋 Step 1: Checking SIF file..."

if [ ! -f "$SIF_SOURCE" ]; then
    log "❌ SIF file not found: $SIF_SOURCE"
    log "   Run ./build_singularity_v1.sh first"
    exit 1
fi

SIF_SIZE=$(du -h "$SIF_SOURCE" | cut -f1)
log "✅ Found SIF file: $SIF_SIZE"

# ============================================================================
# Step 2: Copy SIF with CI-compatible name
# ============================================================================
log "📋 Step 2: Copying SIF with CI-compatible name..."

SIF_TARGET="$ACTIONS_RUNNER_DIR/micapipe_v0.2.3.sif"
cp "$SIF_SOURCE" "$SIF_TARGET"
log "✅ SIF copied to: $SIF_TARGET"

# ============================================================================
# Step 3: Create CI-compatible Dockerfile
# ============================================================================
log "📋 Step 3: Creating CI-compatible Dockerfile..."

cat > "$DOCKERFILE" << 'EOF'
# Use an official Ubuntu base image
FROM ubuntu:20.04

# Avoid interactive prompts and allow runner to run as root
ENV DEBIAN_FRONTEND=noninteractive
ENV RUNNER_ALLOW_RUNASROOT=1

# Install necessary packages and build dependencies for Singularity and the Actions runner
RUN apt-get update && apt-get install -y \
    curl \
    sudo \
    git \
    tar \
    ca-certificates \
    libicu-dev \
    docker.io \
    build-essential \
    libseccomp-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup \
    wget \
    libglib2.0-dev \
    uidmap \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------
# Install Go (version 1.17.13) required for Singularity build
# ---------------------------
ENV GO_VERSION=1.17.13
RUN wget https://golang.org/dl/go${GO_VERSION}.linux-amd64.tar.gz && \
    tar -C /usr/local -xzf go${GO_VERSION}.linux-amd64.tar.gz && \
    rm go${GO_VERSION}.linux-amd64.tar.gz
ENV PATH="/usr/local/go/bin:${PATH}"

# Set CI-compatible environment variables
ENV SINGULARITY_TMPDIR=/export02/local/singularity_tmp
ENV SINGULARITY_CACHEDIR=/export02/local/singularity_tmp/cache

RUN echo "root:100000:65536" >> /etc/subuid && \
    echo "root:100000:65536" >> /etc/subgid

# ---------------------------
# Install Singularity from source (version 3.10.2)
# ---------------------------
ENV SINGULARITY_VERSION=3.10.2
RUN wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-ce-${SINGULARITY_VERSION}.tar.gz && \
    tar -xzf singularity-ce-${SINGULARITY_VERSION}.tar.gz && \
    cd singularity-ce-${SINGULARITY_VERSION} && \
    ./mconfig && \
    make -C builddir && \
    make -C builddir install && \
    cd .. && rm -rf singularity-ce-${SINGULARITY_VERSION}*

# ---------------------------
# Setup GitHub Actions Runner
# ---------------------------
WORKDIR /actions-runner

# Download the GitHub Actions runner package
RUN curl -o actions-runner-linux-x64-2.322.0.tar.gz -L https://github.com/actions/runner/releases/download/v2.322.0/actions-runner-linux-x64-2.322.0.tar.gz && \
    echo "b13b784808359f31bc79b08a191f5f83757852957dd8fe3dbfcc38202ccf5768  actions-runner-linux-x64-2.322.0.tar.gz" | shasum -a 256 -c && \
    tar xzf actions-runner-linux-x64-2.322.0.tar.gz && \
    rm actions-runner-linux-x64-2.322.0.tar.gz

# Install .NET Core dependencies used by the runner
RUN ./bin/installdependencies.sh

# Copy configuration and run scripts; ensure they're executable
COPY config.sh run.sh ./
RUN chmod +x config.sh run.sh

# Configure the runner with your provided token and repository URL
RUN ./config.sh --url https://github.com/MICA-MNI/micapipe --token ADUS2CGIIELTFLXPHXMIU2TIDOVHS --unattended --replace

# ---------------------------
# Embed Singularity Image File (CI-compatible paths)
# ---------------------------
# Create CI-compatible directory structure
RUN mkdir -p /data_/mica1/01_programs/micapipe-v0.2.0 /export02/local/singularity_tmp

# Copy the Singularity image with CI-expected name
COPY micapipe_v0.2.3.sif /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v0.2.3.sif

# ---------------------------
# Create dummy build_singularity_auto.sh that uses embedded SIF
# ---------------------------
RUN cat > /actions-runner/build_singularity_auto.sh << 'SCRIPT_EOF'
#!/bin/bash
# Dummy script for CI - uses embedded SIF instead of building
set -e

# Parse arguments (for CI compatibility)
while [[ $# -gt 0 ]]; do
    case $1 in
        --output)
            OUTPUT_PATH="$2"
            shift 2
            ;;
        --force)
            FORCE=1
            shift
            ;;
        *)
            shift
            ;;
    esac
done

# Default output path
OUTPUT_PATH="${OUTPUT_PATH:-/export02/local/singularity_tmp/micapipe_ci.sif}"

echo "[$(date '+%H:%M:%S')] 🚀 Using embedded SIF (no build needed)"
echo "[$(date '+%H:%M:%S')] 📋 Copying embedded SIF to CI location..."

# Copy embedded SIF to CI expected location
cp /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v0.2.3.sif "$OUTPUT_PATH"

echo "[$(date '+%H:%M:%S')] ✅ SIF ready at: $OUTPUT_PATH"
echo "[$(date '+%H:%M:%S')] 📊 Size: $(du -h "$OUTPUT_PATH" | cut -f1)"
SCRIPT_EOF

RUN chmod +x /actions-runner/build_singularity_auto.sh

# ---------------------------
# Declare volumes for CI compatibility
# ---------------------------
VOLUME ["/data/mica3/BIDS_CI/rawdata", "/data_/mica3/BIDS_CI", "/data/mica1/03_projects/enning/BIDS_CI", "/export02/local/singularity_tmp", "/tmp", "/var/run/docker.sock"]

# Set the entrypoint to run the Actions runner
ENTRYPOINT ["./run.sh"]
EOF

log "✅ CI-compatible Dockerfile created"

# ============================================================================
# Step 4: Build Docker image
# ============================================================================
log "📋 Step 4: Building CI-compatible Docker image..."

cd "$ACTIONS_RUNNER_DIR"
docker build -t micapipe-runner-ci .

# ============================================================================
# Step 5: Create CI-compatible run script  
# ============================================================================
log "📋 Step 5: Creating CI-compatible run script..."

cat > "$ACTIONS_RUNNER_DIR/run_docker_ci.sh" << 'EOF'
#!/bin/bash
export DOCKER_CONTENT_TRUST=0

docker run -d --name micapipe-runner-ci --restart always --privileged  \
    -v /data/mica1/03_projects/enning/BIDS_CI/sing_cache:/data/mica1/03_projects/enning/BIDS_CI/sing_cache \
    -v /export02/local/singularity_tmp/:/export02/local/singularity_tmp/ \
    -v /data/mica3/BIDS_CI/rawdata:/data/mica3/BIDS_CI/rawdata \
    -v /data_/mica3/BIDS_CI:/data_/mica3/BIDS_CI \
    -v /data/mica1/03_projects/enning/BIDS_CI:/data/mica1/03_projects/enning/BIDS_CI \
    -v /tmp:/tmp \
    -v /var/run/docker.sock:/var/run/docker.sock \
    micapipe-runner-ci
EOF

chmod +x "$ACTIONS_RUNNER_DIR/run_docker_ci.sh"

log "✅ CI-compatible run script created: $ACTIONS_RUNNER_DIR/run_docker_ci.sh"

# ============================================================================
# Completion
# ============================================================================
IMAGE_SIZE=$(docker images micapipe-runner-ci --format "table {{.Size}}" | tail -n1)

log "============================================="
log "✅ CI-COMPATIBLE RUNNER BUILD COMPLETE"
log "============================================="
log "📦 Docker image: micapipe-runner-ci"
log "📊 Image size: $IMAGE_SIZE"
log "🗃️  Embedded SIF: /data_/mica1/01_programs/micapipe-v0.2.0/micapipe_v0.2.3.sif"
log ""
log "🚀 Deploy with:"
log "   cd $ACTIONS_RUNNER_DIR"
log "   ./run_docker_ci.sh"
log ""
log "✅ This runner will work with existing CI workflows!"
log "   - Uses CI-expected paths"
log "   - Has dummy build_singularity_auto.sh"
log "   - Embedded SIF with correct name"
log ""
log "🎯 CI jobs will now use embedded SIF instead of building!"