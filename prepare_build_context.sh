#!/bin/bash
# ============================================================================ 
# PREPARE BUILD CONTEXT WITH PRE-DOWNLOADED FILES
# Copy pre-downloaded files to build context so Docker can find them
# ============================================================================

set -euo pipefail

echo "📦 Preparing build context with pre-downloaded files"
echo "=================================================="

# Check if we're in the right location
if [[ ! -f "Dockerfile.mamba-base" ]]; then
    echo "❌ Dockerfile.mamba-base not found"
    echo "   Please run this from your micapipe source directory"
    exit 1
fi

# The files that exist on your server
DOWNLOADS_DIR="/host/cassio/export03/data/enning/downloads"
SOURCE_FILES=(
    "freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
    "fsl-6.0.2-centos6_64.tar.gz"  
    "afni-linux_openmp_64.tgz"
)

# Copy files to current directory so Docker build context can find them
echo "📁 Copying pre-downloaded files to build context..."

for file in "${SOURCE_FILES[@]}"; do
    if [[ -f "$DOWNLOADS_DIR/$file" ]]; then
        echo "✅ Copying $file..."
        cp "$DOWNLOADS_DIR/$file" ./
    else
        echo "⚠️  $file not found at $DOWNLOADS_DIR/$file"
    fi
done

# Also look for FIX and Miniconda (different naming)
echo ""
echo "🔍 Looking for additional files..."

# Check for FIX with various names
if [[ -f "$DOWNLOADS_DIR/fix-1.068.tar.gz" ]]; then
    echo "✅ Copying fix-1.068.tar.gz..."
    cp "$DOWNLOADS_DIR/fix-1.068.tar.gz" ./
elif [[ -f "$DOWNLOADS_DIR/fix.tar.gz" ]]; then
    echo "✅ Copying fix.tar.gz..."
    cp "$DOWNLOADS_DIR/fix.tar.gz" ./fix-1.068.tar.gz
else
    echo "⚠️  FSL FIX not found"
fi

# Check for Miniconda with various names  
if [[ -f "$DOWNLOADS_DIR/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh" ]]; then
    echo "✅ Copying Miniconda3-py39_22.11.1-1-Linux-x86_64.sh..."
    cp "$DOWNLOADS_DIR/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh" ./
elif [[ -f "$DOWNLOADS_DIR/Miniconda3-latest-Linux-x86_64.sh" ]]; then
    echo "✅ Copying Miniconda3-latest-Linux-x86_64.sh..."
    cp "$DOWNLOADS_DIR/Miniconda3-latest-Linux-x86_64.sh" ./Miniconda3-py39_22.11.1-1-Linux-x86_64.sh
else
    echo "⚠️  Miniconda not found"
fi

echo ""
echo "📋 Build context files:"
ls -la *.tar.gz *.tgz *.sh 2>/dev/null || echo "No pre-downloaded files in build context"

echo ""
echo "✅ Build context prepared!"
echo "💡 Now you can run: docker build -f Dockerfile.mamba-base -t micapipe-base ."