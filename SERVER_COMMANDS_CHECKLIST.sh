#!/bin/bash
# ============================================================================
# SERVER-SIDE COMMANDS CHECKLIST
# Copy-paste these commands on your server to build successfully
# ============================================================================

echo "🏗️ SERVER-SIDE BUILD CHECKLIST"
echo "==============================="
echo ""

# Step 1: Navigate to build directory
echo "📍 Step 1: Navigate to server build directory"
echo "   cd /host/cassio/export03/data/enning/downloads"
echo ""

# Step 2: Update code  
echo "📥 Step 2: Update code from git"
echo "   git pull origin comprehensive-base-image"
echo ""

# Step 3: Pre-build verification
echo "🔍 Step 3: Verify setup before building"
echo "   # Check required files:"
echo "   ls -la *.tar.gz *.tgz"
echo "   # Test Docker access:"
echo "   docker info"
echo ""

# Step 4: Build base image
echo "🐳 Step 4: Build comprehensive base image (one-time)"
echo "   docker build -f Dockerfile.mamba-base -t micapipe-comprehensive-base ."
echo ""

# Step 5: Build main image
echo "⚡ Step 5: Build main image (fast daily builds)"  
echo "   docker build -f Dockerfile.minimal -t micapipe:latest ."
echo ""

echo "🎯 Expected timeline:"
echo "   - Base image: 45-77 minutes (one-time setup)"
echo "   - Main image: 3-5 minutes (daily development)"
echo ""

# All-in-one copy-paste section
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📋 COPY-PASTE COMMANDS FOR SERVER:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

cat << 'EOF'
# Navigate to server directory
cd /host/cassio/export03/data/enning/downloads

# Update code
git pull origin comprehensive-base-image

# Verify setup
echo "🔍 Checking pre-downloaded files..."
ls -la *.tar.gz *.tgz 2>/dev/null | head -5
echo "🐳 Testing Docker access..."
docker info >/dev/null && echo "✅ Docker OK" || echo "❌ Docker issue"

# Build comprehensive base (one-time, ~45-77 min)
echo "🏗️ Building comprehensive base image..."
docker build -f Dockerfile.mamba-base -t micapipe-comprehensive-base .

# Build main image (fast, ~3-5 min)  
echo "⚡ Building main image..."
docker build -f Dockerfile.minimal -t micapipe:latest .

echo "✅ Build complete!"
docker images | grep micapipe
EOF

echo ""
echo "✅ Ready! Copy the commands above and run them on your server."