#!/bin/bash
set -euecho "🚚 MICApipe Two-Stage Build - Server Migration"
echo "=============================================="
echo ""
echo "⚠️  IMPORTANT: This script copies files TO the server"
echo "   Docker builds will happen ON THE SERVER at: $BUILD_DIR"
echo "   NOT in your home directory!"
echo ""
echo "📍 Server base: $SERVER_BASE_DIR"
echo "📁 Build directory (ON SERVER): $BUILD_DIR"
echo "🏠 Source code (LOCAL): $HOME_MICAPIPE"
echo ""
echo "📦 Two-Stage Strategy:"
echo "   Stage 1 (Dockerfile.base): Build base with ALL tools (45-90 min, rarely)"
echo "   Stage 2 (Dockerfile.main): Build main with code only (3-5 min, frequently)"
echo ""
echo "💡 Simple approach: Pre-downloaded files + build files in SAME directory"
echo ""# MICApipe Two-Stage Build Strategy - Server Migration Script
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
BUILD_DIR="$DOWNLOADS_DIR"  # Build in same directory as pre-downloaded files (SIMPLER!)
BACKUP_DIR="$SERVER_BASE_DIR/downloads_backup"
HOME_MICAPIPE="$PWD"  # Use current directory as source

echo "🚚 MICApipe Two-Stage Build - Server Migration"
echo "=============================================="
echo "📍 Server base: $SERVER_BASE_DIR"
echo "📁 Build directory: $BUILD_DIR"
echo "🏠 Source code: $HOME_MICAPIPE"
echo ""
echo "📦 Two-Stage Strategy:"
echo "   Stage 1 (Dockerfile.base): Build base with ALL tools (45-90 min, rarely)"
echo "   Stage 2 (Dockerfile.main): Build main with code only (3-5 min, frequently)"
echo ""
echo "� Simple approach: Pre-downloaded files + build files in SAME directory"
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
        echo "⚠️  Missing files detected - will download during build (slower)"
        echo "   Continuing automatically..."
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

# Create build directory if it doesn't exist
if [[ ! -d "$BUILD_DIR" ]]; then
    echo "📁 Creating build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    SOURCE_CHANGED=true
fi

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
    
    # Copy .dockerignore
    echo "   Copying .dockerignore..."
    cp "$HOME_MICAPIPE/.dockerignore" "$BUILD_DIR/"
    
    # Copy build scripts
    echo "   Copying build scripts..."
    cp "$HOME_MICAPIPE/build_comprehensive_base_server.sh" "$BUILD_DIR/"
    cp "$HOME_MICAPIPE/build_main_image_server.sh" "$BUILD_DIR/"
    
    # Copy ONLY necessary directories for Docker build context
    echo "   Copying R_config directory..."
    cp -r "$HOME_MICAPIPE/R_config" "$BUILD_DIR/" 2>/dev/null || true
    
    # Copy other essential config directories
    echo "   Copying essential config directories..."
    for dir in parcellations surfaces MNI152Volumes MICs60_T1-atlas fsl_conf functions; do
        if [[ -d "$HOME_MICAPIPE/$dir" ]]; then
            cp -r "$HOME_MICAPIPE/$dir" "$BUILD_DIR/" 2>/dev/null || true
        fi
    done
    
    # Copy essential scripts
    echo "   Copying essential scripts..."
    find "$HOME_MICAPIPE" -maxdepth 1 -type f \( \
        -name "*.py" -o \
        -name "micapipe" -o \
        -name "*.yaml" -o \
        -name "*.toml" -o \
        -name "*.txt" \
    \) ! -name "*test*" -exec cp {} "$BUILD_DIR/" \; 2>/dev/null || true
    
    # Copy documentation
    if [[ -f "$HOME_MICAPIPE/DOCKERFILE_COMPREHENSIVE_REVIEW.md" ]]; then
        cp "$HOME_MICAPIPE/DOCKERFILE_COMPREHENSIVE_REVIEW.md" "$BUILD_DIR/"
    fi
    if [[ -f "$HOME_MICAPIPE/readme.md" ]]; then
        cp "$HOME_MICAPIPE/readme.md" "$BUILD_DIR/"
    fi
    
    echo ""
    echo "   ✅ Pre-downloaded files already in: $BUILD_DIR"
    echo "   ✅ Build files copied to same directory (SIMPLE!)"
    
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
    "Dockerfile.mamba-base"
    "Dockerfile.main"
    "build_comprehensive_base_server.sh"
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
echo "⚠️  CRITICAL: DO NOT build from ~/micapipe (home directory)!"
echo "   Builds MUST happen from: $BUILD_DIR"
echo ""
echo "📋 Build Workflow:"
echo ""
echo "   STAGE 1 (Rarely - when dependencies change):"
echo "   ------------------------------------------"
echo "   cd $BUILD_DIR"
echo "   ./build_comprehensive_base_server.sh"
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

# Interactive build selection
echo "🤔 What would you like to build next?"
echo ""
echo "Options:"
echo "   1) Build Base Image Only (45-90 min) - Contains all neuroimaging tools"
echo "   2) Build Main Image Only (3-5 min) - Uses existing base + adds micapipe code"
echo "   3) Build Both (Base then Main) - Full rebuild"
echo "   4) Nothing - Just migrate files"
echo ""
echo "💡 Recommendation: If you have recent base image, choose option 2 (Main Only)"
echo ""
read -p "Enter your choice (1/2/3/4): " -n 1 -r BUILD_CHOICE
echo
echo

case $BUILD_CHOICE in
    1)
        echo "🏗️  Building Base Image Only..."
        echo "📍 Build location: $BUILD_DIR"
        echo "⏱️  Expected time: 45-90 minutes"
        echo ""
        pushd "$BUILD_DIR"
        ./build_comprehensive_base_server.sh
        BASE_BUILD_EXIT_CODE=$?
        popd
        
        if [[ $BASE_BUILD_EXIT_CODE -eq 0 ]]; then
            echo "✅ Base image built successfully!"
            echo "� To build main image later: cd $BUILD_DIR && ./build_main_image_server.sh"
        else
            echo "❌ Base image build failed (exit code: $BASE_BUILD_EXIT_CODE)"
        fi
        ;;
    2)
        echo "⚡ Building Main Image Only (Fast!)..."
        echo "📍 Build location: $BUILD_DIR"
        echo "⏱️  Expected time: 3-5 minutes"
        echo ""
        pushd "$BUILD_DIR"
        ./build_main_image_server.sh
        MAIN_BUILD_EXIT_CODE=$?
        popd
        
        if [[ $MAIN_BUILD_EXIT_CODE -eq 0 ]]; then
            echo "✅ Main image built successfully!"
            echo "🎯 micapipe:latest is ready to use!"
        else
            echo "❌ Main image build failed (exit code: $MAIN_BUILD_EXIT_CODE)"
        fi
        ;;
    3)
        echo "🏗️  Building Both Images (Base then Main)..."
        echo "� Build location: $BUILD_DIR"
        echo "⏱️  Expected time: 50-95 minutes total"
        echo ""
        pushd "$BUILD_DIR"
        
        echo "🔧 Step 1: Building base image..."
        ./build_comprehensive_base_server.sh
        BASE_BUILD_EXIT_CODE=$?
        
        if [[ $BASE_BUILD_EXIT_CODE -eq 0 ]]; then
            echo "✅ Base image complete! Starting main image..."
            ./build_main_image_server.sh
            MAIN_BUILD_EXIT_CODE=$?
            
            if [[ $MAIN_BUILD_EXIT_CODE -eq 0 ]]; then
                echo "✅ Both images built successfully!"
                echo "🎯 micapipe:latest is ready to use!"
            else
                echo "❌ Main image build failed (exit code: $MAIN_BUILD_EXIT_CODE)"
            fi
        else
            echo "❌ Base image build failed (exit code: $BASE_BUILD_EXIT_CODE)"
            echo "   Skipping main image build"
        fi
        popd
        ;;
    4)
        echo "📁 Files migrated only - no builds started"
        echo ""
        echo "🚀 To build later:"
        echo "   Base: cd $BUILD_DIR && ./build_comprehensive_base_server.sh"
        echo "   Main: cd $BUILD_DIR && ./build_main_image_server.sh"
        ;;
    *)
        echo "❌ Invalid choice. Files migrated only."
        echo ""
        echo "🚀 To build later:"
        echo "   Base: cd $BUILD_DIR && ./build_comprehensive_base_server.sh"
        echo "   Main: cd $BUILD_DIR && ./build_main_image_server.sh"
        ;;
esac