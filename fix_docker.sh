#!/bin/bash

# Quick Docker Fix for Micapipe Build
# Run this if you get Docker certificate or content trust errors

echo "🔧 Quick Docker Fix"
echo "==================="
echo

# Fix 1: Disable Docker Content Trust
echo "Setting DOCKER_CONTENT_TRUST=0..."
export DOCKER_CONTENT_TRUST=0
echo "✅ Docker Content Trust disabled"

# Fix 2: Test Docker
echo "Testing Docker connection..."
if docker version > /dev/null 2>&1; then
    echo "✅ Docker is working!"
else
    echo "❌ Docker still has issues"
    echo "Try these additional fixes:"
    echo "  1. Restart your terminal session"
    echo "  2. Ask admin to add you to docker group: sudo usermod -aG docker $USER"
    echo "  3. Restart Docker service: sudo systemctl restart docker"
    exit 1
fi

# Fix 3: Add to bashrc for persistence
echo "Adding to ~/.bashrc for future sessions..."
if ! grep -q "DOCKER_CONTENT_TRUST=0" ~/.bashrc; then
    echo "export DOCKER_CONTENT_TRUST=0" >> ~/.bashrc
    echo "✅ Added DOCKER_CONTENT_TRUST=0 to ~/.bashrc"
else
    echo "ℹ️ DOCKER_CONTENT_TRUST already in ~/.bashrc"
fi

echo
echo "🚀 Ready to build! Try running:"
echo "   ./build_no_sudo.sh"
echo
echo "Or test Docker first:"
echo "   docker run --rm hello-world"