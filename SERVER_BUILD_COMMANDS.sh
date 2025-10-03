#!/bin/bash
# ============================================================================
# CORRECTED SERVER BUILD GUIDE
# Copy codebase FROM home TO server (where you have unlimited space)
# ============================================================================

echo "🏗️  MICApipe Server Build Guide (CORRECTED)"
echo "==========================================="
echo ""
echo "📋 STEP 1: Copy codebase from home to server"
echo "   # Run this from your home directory:"
echo "   cd ~/CodeProj/micapipe"
echo "   rsync -av --exclude='.git' ./ /host/cassio/export03/data/enning/downloads/"
echo "   # OR use the script:"
echo "   ./copy_files_to_build_context.sh"
echo ""
echo "📋 STEP 2: Build on server where pre-downloaded files exist"
echo "   cd /host/cassio/export03/data/enning/downloads"
echo "   ./build_comprehensive_base_server.sh"
echo ""
echo "🎯 Why this works:"
echo "   - Server has unlimited space at /host/cassio/export03/data/enning"
echo "   - Pre-downloaded files are ALREADY there (FreeSurfer, FSL, etc.)"
echo "   - Just need to copy source code TO the server"
echo "   - Build uses existing downloads + fresh source code"
echo ""
echo "📁 Pre-downloaded files ALREADY on server:"
echo "   ✅ freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
echo "   ✅ fsl-6.0.2-centos6_64.tar.gz"
echo "   ✅ afni-linux_openmp_64.tgz"
echo "   ✅ fix-1.068.tar.gz (if available)"
echo "   ✅ Miniconda3-py39_22.11.1-1-Linux-x86_64.sh (if available)"
echo ""
echo "🔧 Commands to copy-paste:"

cat << 'EOF'

# === COPY THESE COMMANDS ===

# Step 1: Copy codebase FROM home TO server
cd ~/CodeProj/micapipe
rsync -av --progress --exclude='.git' --exclude='*.tar.gz' --exclude='*.tgz' \
    ./ /host/cassio/export03/data/enning/downloads/

# Step 2: Build on server using pre-downloaded files  
cd /host/cassio/export03/data/enning/downloads
./build_comprehensive_base_server.sh

# === END COMMANDS ===

EOF

echo ""
echo "✅ Ready! Copy the commands above and run them on your server."