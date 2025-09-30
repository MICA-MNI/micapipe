#!/bin/bash

# Quick Setup Script for Server
# Copy and run this on your server to get ready for building

echo "🚀 Micapipe Container Build Setup"
echo "================================="
echo

# Server paths
SIF_DIR="/data_/mica1/01_programs/singularity"
TMP_DIR="/host/cassio/export03/data/enning/tmp"

echo "Setting up directories..."

# Check if directories exist and are accessible
echo "Checking directory accessibility..."

# Check SIF directory
if [ -d "$SIF_DIR" ]; then
    if [ -w "$SIF_DIR" ]; then
        echo "✅ SIF directory exists and is writable: $SIF_DIR"
    else
        echo "❌ SIF directory exists but is not writable: $SIF_DIR"
        echo "   Please ask your admin to run: chmod 755 $SIF_DIR && chown $USER:$USER $SIF_DIR"
        SIF_ACCESS=false
    fi
else
    echo "❌ SIF directory does not exist: $SIF_DIR"
    echo "   Please ask your admin to run: mkdir -p $SIF_DIR && chown $USER:$USER $SIF_DIR"
    SIF_ACCESS=false
fi

# Check TMP directory
if [ -d "$TMP_DIR" ]; then
    if [ -w "$TMP_DIR" ]; then
        echo "✅ TMP directory exists and is writable: $TMP_DIR"
    else
        echo "❌ TMP directory exists but is not writable: $TMP_DIR"
        echo "   Please ask your admin to run: chmod 755 $TMP_DIR && chown $USER:$USER $TMP_DIR"
        TMP_ACCESS=false
    fi
else
    # Try to create it without sudo
    if mkdir -p "$TMP_DIR" 2>/dev/null; then
        echo "✅ Created TMP directory: $TMP_DIR"
        chmod 755 "$TMP_DIR" 2>/dev/null || true
    else
        echo "❌ Cannot create TMP directory: $TMP_DIR"
        echo "   Please ask your admin to run: mkdir -p $TMP_DIR && chown $USER:$USER $TMP_DIR"
        TMP_ACCESS=false
    fi
fi
echo

# Check space
echo "📊 Disk space check:"
echo "SIF directory:"
df -h "$SIF_DIR"
echo "TMP directory:"
df -h "$TMP_DIR"
echo

# Check prerequisites
echo "🔍 Prerequisites check:"

if command -v docker &> /dev/null; then
    echo "✅ Docker: $(docker --version)"
else
    echo "❌ Docker: Not found"
fi

if command -v singularity &> /dev/null; then
    echo "✅ Singularity: $(singularity --version)"
else
    echo "❌ Singularity: Not found"
fi

if command -v git &> /dev/null; then
    echo "✅ Git: $(git --version)"
else
    echo "❌ Git: Not found"
fi

echo
echo "🔧 Next steps:"
echo "1. Clone the repository:"
echo "   git clone https://github.com/MICA-MNI/micapipe.git"
echo "   cd micapipe"
echo "   git checkout docker-container-updates"
echo
echo "2. Run the build:"
echo "   # CPU-only build"
echo "   ./build_server.sh"
echo
echo "   # CUDA-enabled build"
echo "   ./build_server.sh --cuda"
echo
echo "   # Clean build without cache"
echo "   ./build_server.sh --no-cache"
echo

# Set environment variable persistently
echo "📝 Setting environment variables..."
echo "export SINGULARITY_TMPDIR=$TMP_DIR" >> ~/.bashrc
echo "export DOCKER_CONTENT_TRUST=0" >> ~/.bashrc
echo "✅ Added SINGULARITY_TMPDIR and DOCKER_CONTENT_TRUST to ~/.bashrc"

echo
echo "🎯 Ready to build! Make sure to source your bashrc or start a new session:"
echo "   source ~/.bashrc"