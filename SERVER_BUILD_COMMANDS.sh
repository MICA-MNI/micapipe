#!/bin/bash
# ============================================================================
# QUICK SERVER BUILD GUIDE
# Step-by-step commands to copy codebase and build on server
# ============================================================================

echo "🏗️  MICApipe Server Build Guide"
echo "=============================="
echo ""
echo "📋 STEP 1: Get codebase on server"
echo "   cd /host/cassio/export03/data/enning/downloads"
echo "   git clone https://github.com/MICA-MNI/micapipe.git ."
echo "   # OR if already exists:"
echo "   git checkout comprehensive-base-image"
echo "   git pull origin comprehensive-base-image"
echo ""
echo "📋 STEP 2: Build comprehensive base image (one-time, ~45-77 min)"
echo "   ./build_comprehensive_base_server.sh"
echo ""
echo "📋 STEP 3: Build fast main image (daily, ~3-5 min)"
echo "   ./build_fast_ci_server.sh"
echo ""
echo "🎯 Expected Results:"
echo "   - Base image: Uses pre-downloaded FreeSurfer, FSL, AFNI, etc."
echo "   - Main image: Ultra-fast builds with only code changes"
echo "   - 95% build time reduction achieved!"
echo ""
echo "📁 Pre-downloaded files location: /host/cassio/export03/data/enning/downloads"
echo "   ✅ freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
echo "   ✅ fsl-6.0.2-centos6_64.tar.gz"
echo "   ✅ afni-linux_openmp_64.tgz"
echo "   ✅ fix-1.068.tar.gz"
echo "   ✅ Miniconda3-py39_22.11.1-1-Linux-x86_64.sh"
echo ""
echo "🔧 Commands ready to copy-paste for server:"

cat << 'EOF'

# === COPY THESE COMMANDS TO SERVER ===

# Step 1: Navigate to downloads directory (unlimited space)
cd /host/cassio/export03/data/enning/downloads

# Step 2: Get latest codebase
git checkout comprehensive-base-image
git pull origin comprehensive-base-image

# Step 3: Build comprehensive base (one-time, uses pre-downloaded files)
./build_comprehensive_base_server.sh

# Step 4: Build fast main image (daily development)
./build_fast_ci_server.sh

# === END SERVER COMMANDS ===

EOF

echo ""
echo "✅ Ready! Copy the commands above and run them on your server."