#!/bin/bash
set -euo pipefail

# MICApipe Server Migration Script
# ================================
# This script copies source code to the server and sets up the build environment

# Configuration
SERVER_BASE_DIR="/host/cassio/export03/data/enning"
DOWNLOADS_DIR="$SERVER_BASE_DIR/downloads"
BUILD_DIR="$SERVER_BASE_DIR/micapipe_build"
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

# Create build directory if it doesn't exist
if [[ ! -d "$BUILD_DIR" ]]; then
    echo "📁 Creating build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
else
    echo "📁 Build directory exists: $BUILD_DIR"
fi

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

# Copy source files if needed
if [[ "$SOURCE_CHANGED" == "true" ]]; then
    echo "📋 Copying source files to server..."
    
    # Use rsync for efficient copying
    if command -v rsync >/dev/null 2>&1; then
        rsync -av --delete \
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
        cp -r . "$BUILD_DIR/"
        # Clean up unwanted files
        rm -rf "$BUILD_DIR/.git" "$BUILD_DIR/__pycache__" "$BUILD_DIR"/*.pyc "$BUILD_DIR/.DS_Store" "$BUILD_DIR/build_logs" 2>/dev/null || true
    fi
    
    # Mark sync time
    touch "$BUILD_DIR/.last_sync"
    echo "✅ Source files copied successfully"
fi

# Set up downloads directory in build location
echo "📦 Setting up downloads in build directory..."

# Copy downloads directly into build directory (not symlink to avoid any path issues)
if [[ ! -d "$BUILD_DIR/downloads" ]] || [[ -n "$(find "$DOWNLOADS_DIR" -newer "$BUILD_DIR/downloads" -type f 2>/dev/null)" ]]; then
    echo "   📋 Copying downloads to build directory..."
    # Remove existing downloads if it exists
    rm -rf "$BUILD_DIR/downloads" 2>/dev/null || true
    
    # Copy downloads directory
    cp -r "$DOWNLOADS_DIR" "$BUILD_DIR/downloads"
    echo "✅ Downloads copied to build directory"
    
    # Verify copy was successful
    if [[ -f "$BUILD_DIR/downloads/fsl-6.0.2-centos6_64.tar.gz" ]]; then
        BUILD_FSL_SIZE=$(du -h "$BUILD_DIR/downloads/fsl-6.0.2-centos6_64.tar.gz" | cut -f1)
        echo "   ✅ FSL in build: $BUILD_FSL_SIZE"
    else
        echo "   ❌ FSL copy failed!"
        exit 1
    fi
    
    if [[ -f "$BUILD_DIR/downloads/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]]; then
        BUILD_FS_SIZE=$(du -h "$BUILD_DIR/downloads/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" | cut -f1)
        echo "   ✅ FreeSurfer in build: $BUILD_FS_SIZE"
    else
        echo "   ❌ FreeSurfer copy failed!"
        exit 1
    fi
else
    echo "✅ Downloads already up to date in build directory"
fi

# Make scripts executable
chmod +x "$BUILD_DIR"/*.sh 2>/dev/null || true

echo ""
echo "🎯 Migration Complete!"
echo "======================"
echo "Build directory: $BUILD_DIR"
echo "Downloads backup: $BACKUP_DIR"
echo ""
echo "Next steps:"
echo "1. cd $BUILD_DIR"
echo "2. ./build_container_server.sh"
echo ""
echo "Or run: pushd $BUILD_DIR && ./build_container_server.sh && popd"