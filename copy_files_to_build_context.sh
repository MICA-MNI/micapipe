#!/bin/bash
# ============================================================================
# QUICK FIX: Copy pre-downloaded files to build context
# Run this on your server before building Docker image
# ============================================================================

echo "🔧 QUICK FIX: Copy pre-downloaded files to build context"
echo "======================================================"

# Your server setup
DOWNLOADS_DIR="/host/cassio/export03/data/enning/downloads"

echo "📍 Current directory: $PWD"
echo "📁 Downloads directory: $DOWNLOADS_DIR"
echo ""

# Check if files exist and copy them
echo "📦 Copying files to current directory (build context)..."

if [[ -f "$DOWNLOADS_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]]; then
    echo "✅ Copying FreeSurfer..."
    cp "$DOWNLOADS_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ./
else
    echo "❌ FreeSurfer not found"
fi

if [[ -f "$DOWNLOADS_DIR/fsl-6.0.2-centos6_64.tar.gz" ]]; then
    echo "✅ Copying FSL..."
    cp "$DOWNLOADS_DIR/fsl-6.0.2-centos6_64.tar.gz" ./
else
    echo "❌ FSL not found"
fi

if [[ -f "$DOWNLOADS_DIR/afni-linux_openmp_64.tgz" ]]; then
    echo "✅ Copying AFNI..."
    cp "$DOWNLOADS_DIR/afni-linux_openmp_64.tgz" ./
else
    echo "❌ AFNI not found"
fi

# Look for FIX with different names
if [[ -f "$DOWNLOADS_DIR/fix-1.068.tar.gz" ]]; then
    echo "✅ Copying FIX..."
    cp "$DOWNLOADS_DIR/fix-1.068.tar.gz" ./
elif [[ -f "$DOWNLOADS_DIR/fix.tar.gz" ]]; then
    echo "✅ Copying FIX (renaming)..."
    cp "$DOWNLOADS_DIR/fix.tar.gz" ./fix-1.068.tar.gz
else
    echo "❌ FIX not found"
fi

# Look for Miniconda
if [[ -f "$DOWNLOADS_DIR/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh" ]]; then
    echo "✅ Copying Miniconda..."
    cp "$DOWNLOADS_DIR/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh" ./
else
    echo "❌ Miniconda not found"
fi

echo ""
echo "📋 Files in build context:"
ls -la *.tar.gz *.tgz *.sh 2>/dev/null || echo "No files found"

echo ""
echo "✅ Ready! Now run your Docker build command"