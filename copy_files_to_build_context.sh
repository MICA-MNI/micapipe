#!/bin/bash
# ============================================================================
# COPY CODEBASE TO SERVER FOR BUILD
# Copy source code from home to server where pre-downloaded files exist
# ============================================================================

echo "� COPY CODEBASE TO SERVER"
echo "=========================="

# Directions: Home → Server (where you have unlimited space)
HOME_MICAPIPE="$HOME/CodeProj/micapipe"  # Your development location
SERVER_BUILD_DIR="/host/cassio/export03/data/enning/downloads"  # Server with space + pre-downloaded files

echo "📂 Source: $HOME_MICAPIPE"
echo "🎯 Target: $SERVER_BUILD_DIR"
echo ""

# Check source exists
if [[ ! -d "$HOME_MICAPIPE" ]]; then
    echo "❌ Source not found: $HOME_MICAPIPE"
    echo "   Please check your home micapipe location"
    exit 1
fi

# Check target exists and has space
if [[ ! -d "$SERVER_BUILD_DIR" ]]; then
    echo "❌ Server build directory not found: $SERVER_BUILD_DIR"
    echo "   Please ensure server is mounted"
    exit 1
fi

echo "✅ Source directory found"
echo "✅ Server directory accessible"
echo ""

# Copy source code to server build location
echo "📦 Copying source code to server..."
rsync -av --progress \
    --exclude='.git' \
    --exclude='*.tar.gz' \
    --exclude='*.tgz' \
    --exclude='build_logs/' \
    --exclude='Running/' \
    "$HOME_MICAPIPE/" "$SERVER_BUILD_DIR/"

echo ""
echo "✅ Codebase copied to server!"
echo "📁 Build location: $SERVER_BUILD_DIR"
echo ""
echo "� Pre-downloaded files already there:"
ls -la "$SERVER_BUILD_DIR"/*.tar.gz "$SERVER_BUILD_DIR"/*.tgz 2>/dev/null || echo "No pre-downloaded files visible"
echo ""
echo "🚀 Next steps:"
echo "   cd $SERVER_BUILD_DIR"
echo "   ./build_comprehensive_base_server.sh"