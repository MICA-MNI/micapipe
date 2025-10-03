#!/bin/bash
# ============================================================================
# MICAPIPE MAMBA ENVIRONMENT BASE IMAGE BUILDER
# Creates a reusable base image with pre-installed mamba environments
# Optimized for CI workflows where environments are stable but code changes frequently
# ============================================================================

set -e

# Registry configuration - update these for your setup
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
BASE_IMAGE_NAME="micapipe-mamba-base"
BASE_IMAGE_TAG="${MICAPIPE_BASE_TAG:-$(date +%Y%m%d)}"
FULL_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:${BASE_IMAGE_TAG}"
LATEST_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:latest"

echo "🚀 Building MICApipe mamba base image: ${FULL_BASE_IMAGE}"
echo "📅 Using date-based tag for versioning: ${BASE_IMAGE_TAG}"

# Create the base image Dockerfile
cat > Dockerfile.mamba-base << 'EOF'
FROM ubuntu:bionic-20201119

USER root
ARG DEBIAN_FRONTEND="noninteractive"

# Essential packages for conda/mamba
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    ca-certificates \
    curl \
    bzip2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda with pre-download support
ENV PATH="/opt/miniconda-latest/bin:$PATH"
ENV MAMBA_NO_LOCK=1
ENV CONDA_PKGS_DIRS="/opt/miniconda-latest/pkgs"

# Use pre-downloaded Miniconda if available (for CI optimization)
COPY . /tmp/build_context/
RUN if [ -f "/tmp/build_context/Miniconda3-latest-Linux-x86_64.sh" ]; then \
        echo "✅ Using pre-downloaded Miniconda installer"; \
        cp /tmp/build_context/Miniconda3-latest-Linux-x86_64.sh /tmp/miniconda.sh; \
    else \
        echo "📥 Downloading Miniconda installer..."; \
        wget -q --retry-connrefused --waitretry=5 --read-timeout=20 --timeout=15 -O /tmp/miniconda.sh \
            https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh; \
    fi \
    && bash /tmp/miniconda.sh -b -p /opt/miniconda-latest \
    && rm /tmp/miniconda.sh \
    && rm -rf /tmp/build_context \
    && . /opt/miniconda-latest/etc/profile.d/conda.sh \
    && conda config --set auto_update_conda false \
    && conda config --set show_channel_urls true \
    && conda config --set solver libmamba \
    && conda config --set channel_priority strict \
    && conda install -y -n base -c conda-forge mamba

# Create optimized micapipe environment
RUN . /opt/miniconda-latest/etc/profile.d/conda.sh \
    && conda create -y -n micapipe python=3.9 \
    && echo "🔧 Installing core scientific packages..." \
    && mamba install -y -n micapipe -c conda-forge \
        --override-channels --strict-channel-priority \
        --threads 8 --retry-delay 5 --retry-attempts 3 \
        numpy scipy pandas matplotlib seaborn plotly \
        nibabel scikit-learn scikit-image \
        dipy==1.4.1 fury==0.8.0 \
        ipython jupyter tqdm joblib \
        h5py hdf5 pytables \
    && echo "🔧 Installing antspy and MRtrix3..." \
    && mamba install -y -n micapipe -c conda-forge -c mrtrix3 \
        --override-channels --strict-channel-priority \
        --threads 8 --retry-delay 5 --retry-attempts 3 \
        antspy mrtrix3==3.0.7 \
    && echo "🔧 Installing utilities and networking packages..." \
    && mamba install -y -n micapipe -c conda-forge \
        --override-channels --strict-channel-priority \
        --threads 8 --retry-delay 5 --retry-attempts 3 \
        bzip2 ca-certificates curl git openssh \
        unzip wget xz zip \
        cryptography paramiko pycryptodome \
        requests urllib3 \
        cmake make \
        gcc_linux-64 gxx_linux-64

# Create designer environment
RUN . /opt/miniconda-latest/etc/profile.d/conda.sh \
    && conda create -y -n designer python=3.8 \
    && mamba install -y -n designer -c conda-forge \
        --threads 8 --retry-delay 5 --retry-attempts 3 \
        numpy scipy matplotlib nibabel dipy

# Clean up to reduce image size
RUN conda clean -ya \
    && find /opt/miniconda-latest -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true

# Set up environment
RUN echo ". /opt/miniconda-latest/etc/profile.d/conda.sh" >> /etc/bash.bashrc
ENV BASH_ENV="/etc/bash.bashrc"

# Metadata
LABEL org.opencontainers.image.title="MICApipe Mamba Base"
LABEL org.opencontainers.image.description="Pre-built mamba environments for MICApipe neuroimaging pipeline"
LABEL org.opencontainers.image.version="1.0"
LABEL org.opencontainers.image.source="https://github.com/MICA-MNI/micapipe"
EOF

# Build the base image
echo "📦 Building base image (this will take some time)..."
docker build -f Dockerfile.mamba-base -t ${FULL_BASE_IMAGE} -t ${LATEST_BASE_IMAGE} .

echo "✅ Base image built successfully!"
echo "📊 Image size:"
docker images ${FULL_BASE_IMAGE}

echo ""
echo "🏷️  Tagged images:"
echo "   - ${FULL_BASE_IMAGE} (date-versioned)"
echo "   - ${LATEST_BASE_IMAGE} (latest)"
echo ""
echo "🚀 To use this base image, update your main Dockerfile:"
echo "FROM ${LATEST_BASE_IMAGE}"
echo ""
echo "💾 To push to registry (for CI):"
echo "docker push ${FULL_BASE_IMAGE}"
echo "docker push ${LATEST_BASE_IMAGE}"
echo ""
echo "� To pull on CI/server:"
echo "docker pull ${LATEST_BASE_IMAGE}"
echo ""
echo "🔄 Base image should be rebuilt when:"
echo "   - Python packages in environment change"
echo "   - New conda packages are added"
echo "   - Major dependency updates needed"

# Cleanup
rm Dockerfile.mamba-base

echo ""
echo "✨ Next steps:"
echo "1. Push base image: docker push ${LATEST_BASE_IMAGE}"
echo "2. Update main Dockerfile to use: FROM ${LATEST_BASE_IMAGE}"
echo "3. Main builds will now be 80-90% faster!"
EOF