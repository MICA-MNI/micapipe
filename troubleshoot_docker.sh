#!/bin/bash

# Docker Troubleshooting Script for Micapipe Build
# This script helps diagnose and fix common Docker issues on servers

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

echo "🔧 Docker Troubleshooting for Micapipe Build"
echo "============================================="
echo

# Check 1: Docker installation and version
print_info "Step 1: Checking Docker installation..."
if command -v docker &> /dev/null; then
    DOCKER_VERSION=$(docker --version 2>/dev/null || echo "Unknown")
    print_success "Docker found: $DOCKER_VERSION"
else
    print_error "Docker is not installed or not in PATH"
    exit 1
fi

# Check 2: Docker daemon status
print_info "Step 2: Checking Docker daemon..."
if docker info &> /dev/null; then
    print_success "Docker daemon is running"
else
    print_error "Docker daemon is not accessible"
    print_info "Common issues:"
    echo "  1. Docker daemon not running"
    echo "  2. Permission issues (user not in docker group)"
    echo "  3. Docker socket permissions"
    echo "  4. TLS/certificate issues"
    echo
fi

# Check 3: User permissions
print_info "Step 3: Checking user permissions..."
if groups | grep -q docker; then
    print_success "User is in docker group"
else
    print_warning "User is not in docker group"
    print_info "Ask admin to run: sudo usermod -aG docker $USER"
    print_info "Then logout and login again"
fi

# Check 4: Docker socket
print_info "Step 4: Checking Docker socket..."
DOCKER_SOCK="/var/run/docker.sock"
if [ -S "$DOCKER_SOCK" ]; then
    print_success "Docker socket exists: $DOCKER_SOCK"
    if [ -r "$DOCKER_SOCK" ] && [ -w "$DOCKER_SOCK" ]; then
        print_success "Docker socket is readable and writable"
    else
        print_error "Docker socket permissions issue"
        ls -la "$DOCKER_SOCK"
    fi
else
    print_error "Docker socket not found: $DOCKER_SOCK"
fi

# Check 5: Try simple Docker command
print_info "Step 5: Testing basic Docker functionality..."
if docker version &> /dev/null; then
    print_success "Docker version command works"
else
    print_error "Docker version command failed"
    echo "Error details:"
    docker version 2>&1 | head -10
fi

# Check 6: Docker context and environment
print_info "Step 6: Checking Docker context..."
echo "Docker context: $(docker context show 2>/dev/null || echo 'Unknown')"
echo "DOCKER_HOST: ${DOCKER_HOST:-'Not set'}"
echo "DOCKER_CERT_PATH: ${DOCKER_CERT_PATH:-'Not set'}"
echo "DOCKER_TLS_VERIFY: ${DOCKER_TLS_VERIFY:-'Not set'}"

# Check 7: Alternative solutions
echo
print_info "Alternative Solutions:"
echo

echo "🐋 Option 1: Fix Docker (ask your admin)"
echo "  sudo systemctl start docker"
echo "  sudo usermod -aG docker $USER"
echo "  # Then logout/login"
echo

echo "🏗️ Option 2: Use Podman instead of Docker"
echo "  # Install podman (if available)"
echo "  alias docker=podman"
echo "  # Then try the build again"
echo

echo "📦 Option 3: Build on different machine"
echo "  # Build on your local machine or another server"
echo "  # Transfer the resulting .sif file"
echo

echo "☁️ Option 4: Use pre-built container"
echo "  # Download from Docker Hub or GitHub releases"
echo "  singularity pull docker://micapipe/micapipe:latest"
echo

echo "🔄 Option 5: Alternative build approach"
echo "  # Use Singularity to build directly from definition file"
echo "  # (requires creating a Singularity definition file)"
echo

# Quick fixes to try
echo
print_info "Quick fixes to try:"
echo

echo "1. Disable Docker Content Trust (common server issue):"
echo "   export DOCKER_CONTENT_TRUST=0"
echo "   docker version"
echo

echo "2. Restart Docker daemon (if you have permission):"
echo "   sudo systemctl restart docker"
echo

echo "3. Reset Docker environment:"
echo "   unset DOCKER_HOST DOCKER_CERT_PATH DOCKER_TLS_VERIFY"
echo "   export DOCKER_CONTENT_TRUST=0"
echo "   docker version"
echo

echo "3. Try with sudo (temporary fix):"
echo "   sudo docker version"
echo "   # If this works, it's a permissions issue"
echo

echo "4. Check if Docker is running in rootless mode:"
echo "   systemctl --user status docker"
echo

# Generate a simple test command
echo
print_info "Test command to verify Docker works:"
echo "docker run --rm hello-world"
echo

print_info "If Docker is working, retry the micapipe build:"
echo "./build_no_sudo.sh"