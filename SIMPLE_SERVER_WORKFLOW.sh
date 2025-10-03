#!/bin/bash
# ============================================================================
# CORRECTED SERVER WORKFLOW (Based on docker-container-updates approach)
# Just git pull and build - that's it!
# ============================================================================

echo "✅ CORRECTED SERVER WORKFLOW"
echo "============================"
echo ""
echo "🎯 What you actually need to do on server:"
echo ""
echo "1️⃣ Navigate to your server directory:"
echo "   cd /host/cassio/export03/data/enning/downloads"
echo ""
echo "2️⃣ Update your code:"
echo "   git pull origin comprehensive-base-image"
echo ""  
echo "3️⃣ Build using pre-downloaded files:"
echo "   docker build -f Dockerfile.mamba-base -t micapipe-base ."
echo ""
echo "🎯 Why this works:"
echo "   - Code is ALREADY on server (you just git pull to update)"
echo "   - Pre-downloaded files are in same directory"
echo "   - Docker COPY . gets both source + downloads"
echo "   - Build context contains everything needed"
echo ""
echo "📁 Directory structure on server:"
echo "   /host/cassio/export03/data/enning/downloads/"
echo "   ├── Dockerfile.mamba-base          # ← Source code"
echo "   ├── build_comprehensive_base_server.sh"
echo "   ├── freesurfer-*.tar.gz          # ← Pre-downloaded"
echo "   ├── fsl-*.tar.gz                  # ← Pre-downloaded" 
echo "   ├── afni-*.tgz                    # ← Pre-downloaded"
echo "   └── ... (all other files)"
echo ""
echo "🔧 Copy-paste commands for server:"
echo ""

cat << 'EOF'
# === SERVER COMMANDS ===

cd /host/cassio/export03/data/enning/downloads
git pull origin comprehensive-base-image
docker build -f Dockerfile.mamba-base -t micapipe-base .

# === DONE! ===
EOF

echo ""
echo "✅ No complex scripts needed - just git pull and build!"