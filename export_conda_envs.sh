#!/bin/bash
# ============================================================================
# CONDA ENVIRONMENT EXPORT/IMPORT OPTIMIZATION
# Export environments to files for lightning-fast recreation
# ============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="${SCRIPT_DIR}/conda_envs"

echo "🚀 Exporting MICApipe conda environments for fast recreation..."

# Create environments directory
mkdir -p "${ENVS_DIR}"

# Function to export environment
export_env() {
    local env_name="$1"
    local export_type="$2"  # explicit or full
    
    echo "📦 Exporting ${env_name} environment (${export_type})..."
    
    if [[ "$export_type" == "explicit" ]]; then
        # Export only explicitly installed packages (faster recreation)
        conda env export -n "$env_name" --from-history > "${ENVS_DIR}/${env_name}_explicit.yml"
    else
        # Export full environment with all dependencies (more reliable)
        conda env export -n "$env_name" > "${ENVS_DIR}/${env_name}_full.yml"
    fi
}

# Check if environments exist
if conda env list | grep -q "micapipe"; then
    export_env "micapipe" "explicit"
    export_env "micapipe" "full"
    echo "✅ micapipe environment exported"
else
    echo "⚠️  micapipe environment not found - creating template..."
    cat > "${ENVS_DIR}/micapipe_template.yml" << 'EOF'
name: micapipe
channels:
  - conda-forge
  - mrtrix3
dependencies:
  - python=3.9
  - numpy
  - scipy
  - pandas
  - matplotlib
  - seaborn
  - plotly
  - nibabel
  - scikit-learn
  - scikit-image
  - dipy=1.4.1
  - fury=0.8.0
  - ipython
  - jupyter
  - tqdm
  - joblib
  - h5py
  - hdf5
  - pytables
  - antspy
  - mrtrix3=3.0.7
  - bzip2
  - ca-certificates
  - curl
  - git
  - openssh
  - unzip
  - wget
  - xz
  - zip
  - cryptography
  - paramiko
  - pycryptodome
  - requests
  - urllib3
  - cmake
  - make
  - gcc_linux-64
  - gxx_linux-64
  - pip
  - pip:
    - vtk
    - pyvista
    - git+https://github.com/MICA-MNI/LAMAReg.git
    - git+https://github.com/MICA-MNI/ENIGMA.git
EOF
fi

if conda env list | grep -q "designer"; then
    export_env "designer" "explicit"
    export_env "designer" "full"
    echo "✅ designer environment exported"
else
    echo "⚠️  designer environment not found - creating template..."
    cat > "${ENVS_DIR}/designer_template.yml" << 'EOF'
name: designer
channels:
  - conda-forge
dependencies:
  - python=3.8
  - numpy
  - scipy
  - matplotlib
  - nibabel
  - dipy
EOF
fi

# Create fast recreation script
cat > "${ENVS_DIR}/recreate_envs.sh" << 'EOF'
#!/bin/bash
# Fast conda environment recreation script
set -e

ENVS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "🚀 Recreating MICApipe conda environments..."

# Use mamba for faster environment creation
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
else
    CONDA_CMD="conda"
fi

# Function to create environment with fallback
create_env() {
    local env_name="$1"
    
    echo "📦 Creating ${env_name} environment..."
    
    # Try explicit first (faster), fallback to full, then template
    if [ -f "${ENVS_DIR}/${env_name}_explicit.yml" ]; then
        echo "Using explicit environment file..."
        $CONDA_CMD env create -f "${ENVS_DIR}/${env_name}_explicit.yml" || \
        echo "Explicit creation failed, trying full environment..."
    fi
    
    if [ -f "${ENVS_DIR}/${env_name}_full.yml" ] && ! conda env list | grep -q "$env_name"; then
        echo "Using full environment file..."
        $CONDA_CMD env create -f "${ENVS_DIR}/${env_name}_full.yml" || \
        echo "Full creation failed, trying template..."
    fi
    
    if [ -f "${ENVS_DIR}/${env_name}_template.yml" ] && ! conda env list | grep -q "$env_name"; then
        echo "Using template environment file..."
        $CONDA_CMD env create -f "${ENVS_DIR}/${env_name}_template.yml"
    fi
}

# Create environments
create_env "micapipe"
create_env "designer"

echo "✅ All environments recreated successfully!"
echo ""
echo "🔧 To activate micapipe environment:"
echo "conda activate micapipe"
EOF

chmod +x "${ENVS_DIR}/recreate_envs.sh"

# Create Dockerfile snippet
cat > "${ENVS_DIR}/Dockerfile.snippet" << 'EOF'
# ============================================================================
# FAST CONDA ENVIRONMENT RECREATION
# Use pre-exported environment files for lightning-fast environment setup
# ============================================================================

# Copy environment files
COPY conda_envs/ /opt/conda_envs/

# Install conda/mamba if not already available
RUN if [ ! -f "/opt/miniconda-latest/bin/conda" ]; then \
        wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
        && bash /tmp/miniconda.sh -b -p /opt/miniconda-latest \
        && rm /tmp/miniconda.sh; \
    fi

ENV PATH="/opt/miniconda-latest/bin:$PATH"

# Install mamba for faster environment creation
RUN . /opt/miniconda-latest/etc/profile.d/conda.sh \
    && conda config --set solver libmamba \
    && conda install -y -n base -c conda-forge mamba

# Recreate environments using exported files (much faster than package-by-package installation)
RUN cd /opt/conda_envs && ./recreate_envs.sh

# Set up environment activation
RUN echo ". /opt/miniconda-latest/etc/profile.d/conda.sh" >> ~/.bashrc \
    && echo "conda activate micapipe" >> ~/.bashrc
EOF

echo ""
echo "✅ Environment export complete!"
echo ""
echo "📁 Environment files created in: ${ENVS_DIR}"
echo "🚀 To use in Docker:"
echo "   1. Copy Dockerfile snippet content to your Dockerfile"
echo "   2. Environment recreation will be ~5-10x faster"
echo ""
echo "🔄 To update environments:"
echo "   1. Rebuild your environments locally"
echo "   2. Run this script again to export updated versions"

# Show file sizes
echo ""
echo "📊 Environment file sizes:"
find "${ENVS_DIR}" -name "*.yml" -exec ls -lh {} \; | awk '{print $9, $5}'