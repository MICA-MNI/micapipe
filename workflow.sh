#!/bin/bash
# ============================================================================
# MICAPIPE COMPREHENSIVE BASE WORKFLOW
# Simple workflow manager for the comprehensive base image strategy
# Coordinates between ~/micapipe (development) and server deployment
# ============================================================================

set -e

echo "🎯 MICApipe Comprehensive Base Image Workflow"
echo "=============================================="
echo ""

# Configuration
HOME_MICAPIPE="$HOME/micapipe"
SERVER_BASE="/host/cassio/export03/data/enning"
SERVER_DOWNLOADS="$SERVER_BASE/downloads"

# Function to show status
show_status() {
    echo "📍 Current Status:"
    
    # Check home development environment
    if [[ -d "$HOME_MICAPIPE" ]]; then
        pushd "$HOME_MICAPIPE" >/dev/null
        CURRENT_BRANCH=$(git branch --show-current 2>/dev/null || echo "unknown")
        echo "   🏠 Home development: $HOME_MICAPIPE (branch: $CURRENT_BRANCH)"
        popd >/dev/null
    else
        echo "   ❌ Home development not found: $HOME_MICAPIPE"
    fi
    
    # Check server environment
    if [[ -d "$SERVER_BASE" ]]; then
        echo "   🖥️  Server base: $SERVER_BASE ✅"
        if [[ -d "$SERVER_DOWNLOADS" ]]; then
            echo "   📁 Downloads/Build: $SERVER_DOWNLOADS ✅"
            
            # Check for key files
            local files_found=0
            local total_files=7
            [[ -f "$SERVER_DOWNLOADS/fsl-6.0.2-centos6_64.tar.gz" ]] && ((files_found++))
            [[ -f "$SERVER_DOWNLOADS/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]] && ((files_found++))
            [[ -f "$SERVER_DOWNLOADS/Dockerfile.mamba-base" ]] && ((files_found++))
            [[ -f "$SERVER_DOWNLOADS/Dockerfile.minimal" ]] && ((files_found++))
            [[ -f "$SERVER_DOWNLOADS/build_comprehensive_base_server.sh" ]] && ((files_found++))
            [[ -f "$SERVER_DOWNLOADS/build_fast_ci_server.sh" ]] && ((files_found++))
            [[ -f "$SERVER_DOWNLOADS/.last_sync_comprehensive" ]] && ((files_found++))
            
            echo "   📊 Setup completeness: $files_found/$total_files files ready"
        else
            echo "   📁 Downloads/Build: $SERVER_DOWNLOADS ❌"
        fi
    else
        echo "   ❌ Server not accessible: $SERVER_BASE"
    fi
    
    # Check for base image
    if docker images | grep -q "micapipe-comprehensive-base"; then
        echo "   🐳 Base image: Available ✅"
    else
        echo "   🐳 Base image: Not built yet ❌"
    fi
    
    echo ""
}

# Function to show menu
show_menu() {
    echo "🔧 Available Actions:"
    echo "   1) Show detailed status"
    echo "   2) Migrate source to server (from ~/micapipe)"
    echo "   3) Build comprehensive base image (45-90 min, one-time)"
    echo "   4) Build main image (3-5 min, fast CI)"
    echo "   5) Complete setup (migrate + base build)"
    echo "   6) Development workflow (sync + fast build)"
    echo "   7) Open server deployment guide"
    echo "   q) Quit"
    echo ""
}

# Function for action 1: detailed status
detailed_status() {
    echo "📊 Detailed Status Check"
    echo "========================"
    
    show_status
    
    echo "🔍 Checking pre-downloaded files:"
    REQUIRED_FILES=(
        "fsl-6.0.2-centos6_64.tar.gz"
        "freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
        "afni-linux_openmp_64.tgz"
        "fix-1.068.tar.gz"
        "Miniconda3-py39_22.11.1-1-Linux-x86_64.sh"
    )
    
    for file in "${REQUIRED_FILES[@]}"; do
        if [[ -f "$SERVER_DOWNLOADS/$file" ]]; then
            size=$(du -h "$SERVER_DOWNLOADS/$file" 2>/dev/null | cut -f1 || echo "unknown")
            echo "   ✅ $file ($size)"
        else
            echo "   ❌ $file (missing)"
        fi
    done
    
    echo ""
    echo "🔍 Checking Docker images:"
    if docker images | grep -q "micapipe-comprehensive-base"; then
        docker images | grep "micapipe-comprehensive-base" | while read line; do
            echo "   ✅ $line"
        done
    else
        echo "   ❌ No comprehensive base images found"
    fi
    
    if docker images | grep -q "micapipe" | grep -v "comprehensive-base"; then
        docker images | grep "micapipe" | grep -v "comprehensive-base" | head -3 | while read line; do
            echo "   🐳 $line"
        done
    fi
}

# Function for action 2: migrate
migrate_to_server() {
    echo "🚚 Migrating Source to Server"
    echo "============================="
    
    if [[ ! -f "./migrate_comprehensive_base_to_server.sh" ]]; then
        echo "❌ Migration script not found. Please run from ~/micapipe directory."
        return 1
    fi
    
    ./migrate_comprehensive_base_to_server.sh
}

# Function for action 3: build base
build_base_image() {
    echo "🏗️  Building Comprehensive Base Image"
    echo "===================================="
    
    if [[ ! -d "$SERVER_DOWNLOADS" ]]; then
        echo "❌ Server downloads directory not found. Run migration first."
        return 1
    fi
    
    echo "📁 Changing to server build directory: $SERVER_DOWNLOADS"
    pushd "$SERVER_DOWNLOADS"
    
    if [[ -f "./build_comprehensive_base_server.sh" ]]; then
        ./build_comprehensive_base_server.sh
    else
        echo "❌ Base build script not found. Run migration first."
        popd
        return 1
    fi
    
    popd
}

# Function for action 4: build main
build_main_image() {
    echo "⚡ Building Main Image (Fast CI)"
    echo "==============================="
    
    if [[ ! -d "$SERVER_DOWNLOADS" ]]; then
        echo "❌ Server downloads directory not found. Run migration first."
        return 1
    fi
    
    echo "📁 Changing to server build directory: $SERVER_DOWNLOADS"
    pushd "$SERVER_DOWNLOADS"
    
    if [[ -f "./build_fast_ci_server.sh" ]]; then
        ./build_fast_ci_server.sh
    else
        echo "❌ Fast CI build script not found. Run migration first."
        popd
        return 1
    fi
    
    popd
}

# Function for action 5: complete setup
complete_setup() {
    echo "🎯 Complete Setup (Migration + Base Build)"
    echo "=========================================="
    
    echo "Step 1: Migrating source to server..."
    migrate_to_server || return 1
    
    echo ""
    echo "Step 2: Building comprehensive base image..."
    build_base_image || return 1
    
    echo ""
    echo "✅ Complete setup finished!"
    echo "🚀 Ready for fast CI builds with ./workflow.sh -> option 4"
}

# Function for action 6: development workflow
dev_workflow() {
    echo "⚡ Development Workflow (Sync + Fast Build)"
    echo "=========================================="
    
    echo "Step 1: Syncing latest changes to server..."
    migrate_to_server || return 1
    
    echo ""
    echo "Step 2: Fast CI build..."
    build_main_image || return 1
    
    echo ""
    echo "✅ Development workflow complete!"
    echo "🐳 Ready to test: docker run --rm -it micapipe:latest"
}

# Function for action 7: open guide
open_guide() {
    if [[ -f "./SERVER_DEPLOYMENT_GUIDE.md" ]]; then
        echo "📖 Opening SERVER_DEPLOYMENT_GUIDE.md..."
        if command -v code >/dev/null 2>&1; then
            code ./SERVER_DEPLOYMENT_GUIDE.md
        elif command -v less >/dev/null 2>&1; then
            less ./SERVER_DEPLOYMENT_GUIDE.md
        else
            cat ./SERVER_DEPLOYMENT_GUIDE.md
        fi
    else
        echo "❌ Server deployment guide not found in current directory"
    fi
}

# Main script
if [[ $# -eq 0 ]]; then
    # Interactive mode
    while true; do
        show_status
        show_menu
        read -p "Choose an action (1-7, q): " choice
        echo ""
        
        case $choice in
            1) detailed_status ;;
            2) migrate_to_server ;;
            3) build_base_image ;;
            4) build_main_image ;;
            5) complete_setup ;;
            6) dev_workflow ;;
            7) open_guide ;;
            q|Q) echo "👋 Goodbye!"; exit 0 ;;
            *) echo "❌ Invalid choice. Please select 1-7 or q." ;;
        esac
        
        echo ""
        read -p "Press Enter to continue..."
        echo ""
    done
else
    # Command line mode
    case $1 in
        status) show_status ;;
        detailed-status) detailed_status ;;
        migrate) migrate_to_server ;;
        build-base) build_base_image ;;
        build-main) build_main_image ;;
        setup) complete_setup ;;
        dev) dev_workflow ;;
        guide) open_guide ;;
        *) 
            echo "Usage: $0 [status|detailed-status|migrate|build-base|build-main|setup|dev|guide]"
            echo "   Or run without arguments for interactive mode"
            exit 1
            ;;
    esac
fi