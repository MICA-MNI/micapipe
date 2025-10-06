#!/bin/bash
set -euo pipefail

# MICApipe Two-Stage Build Strategy - Server Migration Script
# ============================================================
# This script migrates micapipe code to server for two-stage builds:
#   STAGE 1: Build comprehensive base image (Dockerfile.base) - rarely
#   STAGE 2: Build fast main image (Dockerfile.main) - frequently
# 
# Copies from: ~/micapipe (local development)
# Copies to: /host/cassio/export03/data/enning/downloads (server with space + pre-downloaded files)

# Configuration
SERVER_BASE_DIR="/host/cassio/export03/data/enning"
DOWNLOADS_DIR="$SERVER_BASE_DIR/downloads"
BUILD_DIR="$DOWNLOADS_DIR"  # Build directly in downloads directory (has unlimited space)
BACKUP_DIR="$SERVER_BASE_DIR/downloads_backup"
HOME_MICAPIPE="$PWD"  # Use current directory as source

echo "🚚 MICApipe Two-Stage Build - Server Migration"
echo "=============================================="
echo "📍 Server base: $SERVER_BASE_DIR"
echo "📁 Downloads/Build: $DOWNLOADS_DIR"
echo "🏠 Source code: $HOME_MICAPIPE"
echo ""
echo "📦 Two-Stage Strategy:"
echo "   Stage 1 (Dockerfile.base): Build base with ALL tools (45-90 min, rarely)"
echo "   Stage 2 (Dockerfile.main): Build main with code only (3-5 min, frequently)"
echo ""

# Verify source micapipe directory exists
if [[ ! -d "$HOME_MICAPIPE" ]]; then
    echo "❌ Source micapipe directory not found: $HOME_MICAPIPE"
    echo "   Please run this script from your micapipe directory"
    exit 1
fi

# Check if we're on the right branch
pushd "$HOME_MICAPIPE"
CURRENT_BRANCH=$(git branch --show-current)
if [[ "$CURRENT_BRANCH" != "comprehensive-base-image" ]]; then
    echo "⚠️  Warning: Not on comprehensive-base-image branch (current: $CURRENT_BRANCH)"
    echo "   Switching to comprehensive-base-image branch..."
    git checkout comprehensive-base-image
    git pull origin comprehensive-base-image
fi
popd

# Check if server directory is accessible
if [[ ! -d "$SERVER_BASE_DIR" ]]; then
    echo "❌ Server directory not accessible: $SERVER_BASE_DIR"
    echo "   Please ensure the server is mounted or the path is correct"
    exit 1
fi

echo "✅ Server directory accessible: $SERVER_BASE_DIR"

# Check and verify downloads directory with pre-downloaded files
echo "🔍 Checking downloads directory and pre-downloaded files..."
if [[ -d "$DOWNLOADS_DIR" ]]; then
    echo "✅ Downloads directory found: $DOWNLOADS_DIR"
    
    # Check for key files
    REQUIRED_FILES=(
        "fsl-6.0.2-centos6_64.tar.gz"
        "freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
        "afni-linux_openmp_64.tgz"
        "fix-1.068.tar.gz"
        "Miniconda3-py39_22.11.1-1-Linux-x86_64.sh"
    )
    
    MISSING_FILES=()
    for file in "${REQUIRED_FILES[@]}"; do
        FILEPATH="$DOWNLOADS_DIR/$file"
        if [[ -f "$FILEPATH" ]]; then
            FILE_SIZE=$(du -h "$FILEPATH" | cut -f1)
            echo "   ✅ $file: $FILE_SIZE"
        else
            echo "   ❌ $file: NOT FOUND"
            MISSING_FILES+=("$file")
        fi
    done
    
    if [[ ${#MISSING_FILES[@]} -gt 0 ]]; then
        echo ""
        echo "⚠️  Missing pre-downloaded files:"
        for file in "${MISSING_FILES[@]}"; do
            echo "     - $file"
        done
        echo ""
        echo "Build will continue but missing files will be downloaded (slower)."
        read -p "Continue anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "❌ Migration cancelled"
            exit 1
        fi
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
    echo "   Please ensure pre-downloaded files are available at the server location"
    exit 1
fi

# Check for source file changes
SOURCE_CHANGED=false
if [[ ! -f "$BUILD_DIR/.last_sync_twostage" ]]; then
    SOURCE_CHANGED=true
    echo "🔄 First time two-stage setup - copying all source files"
elif [[ -n "$(find "$HOME_MICAPIPE" -name '*.sh' -o -name 'Dockerfile*' -o -name '*.py' -o -name '*.md' -newer "$BUILD_DIR/.last_sync_twostage" 2>/dev/null)" ]]; then
    SOURCE_CHANGED=true
    echo "🔄 Source files changed - updating build directory"
else
    echo "✅ Source files up to date"
fi

if $SOURCE_CHANGED; then
    echo "📋 Copying two-stage build files to server..."
    
    # Copy two-stage Dockerfiles
    echo "   Copying Dockerfile.base (Stage 1 - comprehensive base)..."
    cp "$HOME_MICAPIPE/Dockerfile.base" "$BUILD_DIR/"
    
    echo "   Copying Dockerfile.main (Stage 2 - fast micapipe)..."
    cp "$HOME_MICAPIPE/Dockerfile.main" "$BUILD_DIR/"
    
    # Copy build scripts
    echo "   Copying build scripts..."
    cp "$HOME_MICAPIPE/build_base_image_server.sh" "$BUILD_DIR/"
    cp "$HOME_MICAPIPE/build_main_image_server.sh" "$BUILD_DIR/"
    
    # Copy documentation
    if [[ -f "$HOME_MICAPIPE/DOCKERFILE_COMPREHENSIVE_REVIEW.md" ]]; then
        cp "$HOME_MICAPIPE/DOCKERFILE_COMPREHENSIVE_REVIEW.md" "$BUILD_DIR/"
    fi
    
    # Copy other essential files using rsync for efficiency
    if command -v rsync >/dev/null 2>&1; then
        echo "   Syncing all source files..."
        rsync -av \
            --exclude='.git' \
            --exclude='downloads' \
            --exclude='__pycache__' \
            --exclude='*.pyc' \
            --exclude='.DS_Store' \
            --exclude='build_logs' \
            --exclude='test_data' \
            "$HOME_MICAPIPE"/ "$BUILD_DIR/"
    else
        # Fallback to cp if rsync not available
        echo "   Copying files (rsync not available)..."
        find "$HOME_MICAPIPE" -maxdepth 1 -type f \( -name "*.sh" -o -name "Dockerfile*" -o -name "*.md" -o -name "*.py" \) -exec cp {} "$BUILD_DIR/" \;
        find "$HOME_MICAPIPE" -maxdepth 1 -type d ! -name '.' ! -name 'downloads' ! -name '.git' -exec cp -r {} "$BUILD_DIR/" \;
    fi
    
    # Mark sync time
    touch "$BUILD_DIR/.last_sync_twostage"
    echo "✅ Two-stage build files copied to server"
else
    echo "✅ Source files already up to date"
fi

# Make scripts executable
chmod +x "$BUILD_DIR"/*.sh 2>/dev/null || true

# Verify critical files are in place
echo "🔍 Verifying two-stage build setup..."
STRATEGY_FILES=(
    "Dockerfile.base"
    "Dockerfile.main"
    "build_base_image_server.sh"
    "build_main_image_server.sh"
)

for file in "${STRATEGY_FILES[@]}"; do
    if [[ -f "$BUILD_DIR/$file" ]]; then
        echo "   ✅ $file"
    else
        echo "   ❌ $file - MISSING!"
        echo "Error: Critical file missing. Migration failed."
        exit 1
    fi
done

echo ""
echo "🎯 Migration Complete!"
echo "======================"
echo "Build directory: $BUILD_DIR"
echo "Strategy: Two-Stage Build (95% faster CI builds)"
echo ""
echo "📋 Build Workflow:"
echo ""
echo "   STAGE 1 (Rarely - when dependencies change):"
echo "   ------------------------------------------"
echo "   cd $BUILD_DIR"
echo "   ./build_base_image_server.sh"
echo "   ⏱️  Time: 45-90 minutes"
echo "   📦 Result: micapipe-base:latest with ALL tools"
echo ""
echo "   STAGE 2 (Frequently - every code change):"
echo "   ------------------------------------------"
echo "   cd $BUILD_DIR"
echo "   ./build_main_image_server.sh"
echo "   ⏱️  Time: 3-5 minutes (95% faster!)"
echo "   📦 Result: micapipe:latest ready to use"
echo ""
echo "⚡ Performance Benefits:"
echo "   - Traditional build: 45-90 minutes every time"
echo "   - Base image build: 45-90 minutes (one-time only)"
echo "   - Main image build: 3-5 minutes (every code change)"
echo "   - Time saved per CI run: ~60 minutes!"
echo ""

# Ask if user wants to start the base image build immediately
read -p "🚀 Start Stage 1 base image build now? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "🏗️  Starting Stage 1: Base image build..."
    echo "⏱️  Expected time: 45-90 minutes (one-time setup)"
    echo "📋 This will build the base image with all neuroimaging tools:"
    echo "   - FSL 6.0.2, FreeSurfer 7.4.1, AFNI"
    echo "   - MRtrix3 3.0.7, FastSurfer 2.4.2"
    echo "   - DESIGNER, LAMAReg, SWM, Synb0/SynBOLD"
    echo "   - Conda/Mamba environments"
    echo ""
    
    # Change to build directory and run base image build
    pushd "$BUILD_DIR"
    ./build_base_image_server.sh
    BASE_BUILD_EXIT_CODE=$?
    popd
    
    if [[ $BASE_BUILD_EXIT_CODE -eq 0 ]]; then
        echo ""
        echo "✅ Stage 1 complete! Base image built successfully!"
        echo "🎯 Ready for Stage 2 fast builds!"
        echo ""
        echo "🚀 To build Stage 2 main image (3-5 minutes):"
        echo "   cd $BUILD_DIR && ./build_main_image_server.sh"
    else
        echo ""
        echo "❌ Base image build failed (exit code: $BASE_BUILD_EXIT_CODE)"
        echo "📋 Check build logs for details"
        echo "💡 You can retry later with: cd $BUILD_DIR && ./build_comprehensive_base_server.sh"
    fi
else
    echo ""
    echo "✅ Migration complete - ready for manual builds"
    echo ""
    echo "🔧 When ready to build:"
    echo "   cd $BUILD_DIR"
    echo "   ./build_comprehensive_base_server.sh  # One-time base build"
    echo "   ./build_fast_ci_server.sh            # Fast CI builds"
fi