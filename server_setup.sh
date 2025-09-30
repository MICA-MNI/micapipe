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

# Create directories with proper permissions
sudo mkdir -p "$SIF_DIR" 2>/dev/null || true
sudo mkdir -p "$TMP_DIR" 2>/dev/null || true

# Set ownership
sudo chown $USER:$USER "$SIF_DIR"
sudo chown $USER:$USER "$TMP_DIR"

# Set permissions
chmod 755 "$SIF_DIR"
chmod 755 "$TMP_DIR"

echo "✅ Directories created:"
echo "   SIF: $SIF_DIR"
echo "   TMP: $TMP_DIR"
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
echo "📝 Setting SINGULARITY_TMPDIR environment variable..."
echo "export SINGULARITY_TMPDIR=$TMP_DIR" >> ~/.bashrc
echo "✅ Added to ~/.bashrc"

echo
echo "🎯 Ready to build! Make sure to source your bashrc or start a new session:"
echo "   source ~/.bashrc"