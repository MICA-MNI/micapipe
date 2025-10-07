#!/bin/bash
set -euo pipefail

# MICApipe Server Space Cleanup Script
# =====================================
# Checks disk usage and cleans up unnecessary files on server

SERVER_BASE="/host/cassio/export03/data/enning"
DOWNLOADS_DIR="$SERVER_BASE/downloads"

echo "🧹 MICApipe Server Space Cleanup"
echo "================================="
echo ""

# Check if we're on the server
if [[ ! -d "$SERVER_BASE" ]]; then
    echo "❌ Server directory not accessible: $SERVER_BASE"
    echo "   Please run this script from the server"
    exit 1
fi

# Show current disk usage
echo "📊 Current Disk Usage:"
echo "---------------------"
df -h "$SERVER_BASE" | tail -1
echo ""

# Show space used by main directories
echo "📁 Space Used by Directories:"
echo "-----------------------------"
if [[ -d "$DOWNLOADS_DIR" ]]; then
    echo "Downloads directory:"
    du -sh "$DOWNLOADS_DIR"
    echo ""
fi

# Check Docker disk usage
echo "🐳 Docker Disk Usage:"
echo "--------------------"
if command -v docker >/dev/null 2>&1; then
    docker system df || echo "Cannot access Docker"
else
    echo "Docker not found"
fi
echo ""

# Interactive cleanup options
echo "🗑️  Cleanup Options:"
echo "-------------------"
echo ""

read -p "1. Clean Docker cache (images, containers, build cache)? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "   Cleaning Docker system..."
    docker system prune -a -f --volumes || echo "Failed to clean Docker"
    echo "   ✅ Docker cleaned"
    echo ""
fi

read -p "2. Remove old Docker build logs? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [[ -d "$DOWNLOADS_DIR" ]]; then
        echo "   Removing build logs..."
        rm -f "$DOWNLOADS_DIR"/build_*.log 2>/dev/null || true
        rm -f "$DOWNLOADS_DIR"/*.log 2>/dev/null || true
        echo "   ✅ Build logs removed"
    fi
    echo ""
fi

read -p "3. Remove old micapipe_build subdirectory (if exists)? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [[ -d "$DOWNLOADS_DIR/micapipe_build" ]]; then
        echo "   Removing micapipe_build directory..."
        rm -rf "$DOWNLOADS_DIR/micapipe_build"
        echo "   ✅ micapipe_build removed"
    else
        echo "   ℹ️  micapipe_build directory not found"
    fi
    echo ""
fi

read -p "4. Remove temporary download directories? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [[ -d "$DOWNLOADS_DIR" ]]; then
        echo "   Removing temp directories..."
        rm -rf "$DOWNLOADS_DIR"/temp_* 2>/dev/null || true
        rm -rf "$DOWNLOADS_DIR"/tmp_* 2>/dev/null || true
        rm -rf "$DOWNLOADS_DIR"/.cache 2>/dev/null || true
        echo "   ✅ Temp directories removed"
    fi
    echo ""
fi

read -p "5. List large files in downloads directory? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if [[ -d "$DOWNLOADS_DIR" ]]; then
        echo "   Large files (>100MB):"
        find "$DOWNLOADS_DIR" -maxdepth 1 -type f -size +100M -exec du -h {} \; | sort -h || true
    fi
    echo ""
fi

# Check /tmp and /var/tmp
echo "📁 Temporary Directories:"
echo "------------------------"
du -sh /tmp 2>/dev/null || echo "/tmp: Cannot access"
du -sh /var/tmp 2>/dev/null || echo "/var/tmp: Cannot access"
echo ""

read -p "6. Clean system temporary directories (/tmp, /var/tmp)? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "   ⚠️  This requires sudo permissions"
    sudo rm -rf /tmp/* 2>/dev/null || echo "Cannot clean /tmp"
    sudo rm -rf /var/tmp/* 2>/dev/null || echo "Cannot clean /var/tmp"
    echo "   ✅ System temp directories cleaned"
    echo ""
fi

# Check Docker's internal storage
if command -v docker >/dev/null 2>&1; then
    read -p "7. Check Docker's internal storage location? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        DOCKER_ROOT=$(docker info 2>/dev/null | grep "Docker Root Dir" | cut -d: -f2 | xargs)
        if [[ -n "$DOCKER_ROOT" ]]; then
            echo "   Docker Root: $DOCKER_ROOT"
            sudo du -sh "$DOCKER_ROOT" 2>/dev/null || echo "Cannot access Docker root"
            echo ""
            
            read -p "   Clean Docker root directory? (DANGEROUS - stops all containers) (y/N): " -n 1 -r
            echo
            if [[ $REPLY =~ ^[Yy]$ ]]; then
                echo "   Stopping all containers..."
                docker stop $(docker ps -aq) 2>/dev/null || true
                echo "   Removing all containers..."
                docker rm $(docker ps -aq) 2>/dev/null || true
                echo "   Removing all images..."
                docker rmi $(docker images -q) 2>/dev/null || true
                echo "   Pruning system..."
                docker system prune -a -f --volumes
                echo "   ✅ Docker completely cleaned"
            fi
        fi
        echo ""
    fi
fi

# Final disk usage report
echo ""
echo "📊 Final Disk Usage:"
echo "-------------------"
df -h "$SERVER_BASE" | tail -1
echo ""

# Show space used by downloads directory
if [[ -d "$DOWNLOADS_DIR" ]]; then
    echo "Downloads directory after cleanup:"
    du -sh "$DOWNLOADS_DIR"
    echo ""
fi

echo "✅ Cleanup complete!"
echo ""
echo "💡 Tips to free more space:"
echo "   - Remove unused Docker images: docker image prune -a"
echo "   - Remove old container layers: docker builder prune"
echo "   - Check large files: du -ah $DOWNLOADS_DIR | sort -h | tail -20"
echo "   - Monitor during build: watch -n 5 df -h"
