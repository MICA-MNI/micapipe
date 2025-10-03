#!/bin/bash
set -euo pipefail

# MICApipe Server Migration Script
# ================================
# This script copies source code to the server and sets up the build environment

# Configuration
SERVER_BASE_DIR="/host/cassio/export03/data/enning"
DOWNLOADS_DIR="$SERVER_BASE_DIR/downloads"
BUILD_DIR="$DOWNLOADS_DIR"  # Build directly in downloads directory!
BACKUP_DIR="$SERVER_BASE_DIR/downloads_backup"

echo "🚚 MICApipe Server Migration"
echo "============================"

# Check if server directory is accessible
if [[ ! -d "$SERVER_BASE_DIR" ]]; then
    echo "❌ Server directory not accessible: $SERVER_BASE_DIR"
    echo "   Please ensure the server is mounted or the path is correct"
    exit 1
fi

echo "✅ Server directory accessible: $SERVER_BASE_DIR"

# Check and backup downloads directory first
echo "🔍 Checking downloads directory..."
if [[ -d "$DOWNLOADS_DIR" ]]; then
    echo "✅ Downloads directory found: $DOWNLOADS_DIR"
    
    # Check for key files
    FSL_FILE="$DOWNLOADS_DIR/fsl-6.0.2-centos6_64.tar.gz"
    FS_FILE="$DOWNLOADS_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
    AFNI_FILE="$DOWNLOADS_DIR/afni-linux_openmp_64.tgz"
    FSL_FIX_FILE="$DOWNLOADS_DIR/fix-1.068.tar.gz"
    MINICONDA_FILE="$DOWNLOADS_DIR/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh"
    
    if [[ -f "$FSL_FILE" ]]; then
        FSL_SIZE=$(du -h "$FSL_FILE" | cut -f1)
        echo "   ✅ FSL found: $FSL_SIZE"
    else
        echo "   ❌ FSL not found: $FSL_FILE"
        echo "   Please run ./download_dependencies.sh first"
        exit 1
    fi
    
    if [[ -f "$FS_FILE" ]]; then
        FS_SIZE=$(du -h "$FS_FILE" | cut -f1)
        echo "   ✅ FreeSurfer found: $FS_SIZE"
    else
        echo "   ❌ FreeSurfer not found: $FS_FILE"
        echo "   Please run ./download_dependencies.sh first"
        exit 1
    fi
    
    if [[ -f "$AFNI_FILE" ]]; then
        AFNI_SIZE=$(du -h "$AFNI_FILE" | cut -f1)
        echo "   ✅ AFNI found: $AFNI_SIZE"
    else
        echo "   ⚠️  AFNI not found: $AFNI_FILE (will download during build)"
    fi
    
    if [[ -f "$FSL_FIX_FILE" ]]; then
        FSL_FIX_SIZE=$(du -h "$FSL_FIX_FILE" | cut -f1)
        echo "   ✅ FSL FIX found: $FSL_FIX_SIZE"
    else
        echo "   ⚠️  FSL FIX not found: $FSL_FIX_FILE (will download during build)"
    fi
    
    if [[ -f "$MINICONDA_FILE" ]]; then
        MINICONDA_SIZE=$(du -h "$MINICONDA_FILE" | cut -f1)
        echo "   ✅ Miniconda found: $MINICONDA_SIZE"
    else
        echo "   ⚠️  Miniconda not found: $MINICONDA_FILE (will download during build)"
    fi
    
    # Create backup if it doesn't exist
    if [[ ! -d "$BACKUP_DIR" ]]; then
        echo "💾 Creating backup of downloads directory..."
        cp -r "$DOWNLOADS_DIR" "$BACKUP_DIR"
        echo "✅ Backup created: $BACKUP_DIR"
    else
        echo "✅ Backup already exists: $BACKUP_DIR"
    fi
    
else
    echo "❌ Downloads directory not found: $DOWNLOADS_DIR"
    echo "   Please run ./download_dependencies.sh first"
    exit 1
fi

# Create build directory if it doesn't exist (same as downloads)
echo "📁 Using downloads directory as build directory: $BUILD_DIR"

# Check for existing source and decide if we need to copy
SOURCE_CHANGED=false
if [[ ! -f "$BUILD_DIR/.last_sync" ]]; then
    SOURCE_CHANGED=true
    echo "🔄 First time setup - copying all source files"
elif [[ -n "$(find . -name '*.sh' -o -name 'Dockerfile' -o -name '*.py' -o -name '*.md' -newer "$BUILD_DIR/.last_sync" 2>/dev/null)" ]]; then
    SOURCE_CHANGED=true
    echo "🔄 Source files changed - updating build directory"
else
    echo "✅ Source files up to date"
fi

# Copy source files to downloads directory (excluding downloads to avoid recursion)
echo "📋 Copying source files to downloads directory..."

# Use rsync for efficient copying, but exclude downloads directory itself
if command -v rsync >/dev/null 2>&1; then
    rsync -av \
        --exclude='.git' \
        --exclude='downloads' \
        --exclude='__pycache__' \
        --exclude='*.pyc' \
        --exclude='.DS_Store' \
        --exclude='build_logs' \
        --exclude='test_data' \
        ./ "$BUILD_DIR/"
else
    # Fallback to cp if rsync not available
    find . -maxdepth 1 -type f -exec cp {} "$BUILD_DIR/" \;
    find . -maxdepth 1 -type d ! -name '.' ! -name 'downloads' ! -name '.git' -exec cp -r {} "$BUILD_DIR/" \;
fi

# Mark sync time
touch "$BUILD_DIR/.last_sync"
echo "✅ Source files copied to downloads directory"

# Downloads are already here - no need to copy them!
echo "✅ Downloads already in build directory (same location)"

# Verify downloads are accessible
if [[ -f "$BUILD_DIR/fsl-6.0.2-centos6_64.tar.gz" ]]; then
    BUILD_FSL_SIZE=$(du -h "$BUILD_DIR/fsl-6.0.2-centos6_64.tar.gz" | cut -f1)
    echo "   ✅ FSL ready: $BUILD_FSL_SIZE"
else
    echo "   ❌ FSL not found in build directory!"
    exit 1
fi

if [[ -f "$BUILD_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]]; then
    BUILD_FS_SIZE=$(du -h "$BUILD_DIR/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" | cut -f1)
    echo "   ✅ FreeSurfer ready: $BUILD_FS_SIZE"
else
    echo "   ❌ FreeSurfer not found in build directory!"
    exit 1
fi

if [[ -f "$BUILD_DIR/afni-linux_openmp_64.tgz" ]]; then
    BUILD_AFNI_SIZE=$(du -h "$BUILD_DIR/afni-linux_openmp_64.tgz" | cut -f1)
    echo "   ✅ AFNI ready: $BUILD_AFNI_SIZE"
else
    echo "   ⚠️  AFNI not found in build directory (will download during build)"
fi

if [[ -f "$BUILD_DIR/fix-1.068.tar.gz" ]]; then
    BUILD_FSL_FIX_SIZE=$(du -h "$BUILD_DIR/fix-1.068.tar.gz" | cut -f1)
    echo "   ✅ FSL FIX ready: $BUILD_FSL_FIX_SIZE"
else
    echo "   ⚠️  FSL FIX not found in build directory (will download during build)"
fi

# Make scripts executable
chmod +x "$BUILD_DIR"/*.sh 2>/dev/null || true

echo ""
echo "🎯 Migration Complete!"
echo "======================"
echo "Build/Downloads directory: $BUILD_DIR"
echo "Downloads backup: $BACKUP_DIR"
echo ""
echo "🚀 Starting Docker build..."
echo ""

# Automatically run the build using pushd/popd
pushd "$BUILD_DIR"
./build_container_server.sh
BUILD_EXIT_CODE=$?
popd

if [[ $BUILD_EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "✅ Container build completed successfully!"
    echo "🐳 Container: micapipe:latest"
else
    echo ""
    echo "❌ Container build failed (exit code: $BUILD_EXIT_CODE)"
    echo "📋 Check build logs in: $BUILD_DIR/build_logs/"
fi