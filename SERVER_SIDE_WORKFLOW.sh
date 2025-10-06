#!/bin/bash
# ============================================================================
# SERVER-SIDE WORKFLOW: Copy from ~/micapipe to /host and build
# Run this script from ~/micapipe on the SERVER
# ============================================================================

echo "🚀 SERVER-SIDE MICAPIPE BUILD WORKFLOW"
echo "====================================="
echo ""
echo "📍 SERVER Setup:"
echo "   🏠 Development: ~/micapipe (where you git pull and develop)"  
echo "   🎯 Build location: /host/cassio/export03/data/enning/downloads (unlimited space)"
echo "   📦 Pre-downloaded files: Already at /host/cassio/export03/data/enning/downloads"
echo ""
echo "🔄 Complete SERVER workflow:"
echo ""
echo "1️⃣ Update your code (in ~/micapipe):"
echo "   cd ~/micapipe"
echo "   git pull origin comprehensive-base-image"
echo ""
echo "2️⃣ Copy code to build location with unlimited space:"
echo "   ./migrate_comprehensive_base_to_server.sh"
echo "   # This copies ~/micapipe → /host/cassio/export03/data/enning/downloads"
echo ""
echo "3️⃣ Build using copied code + pre-downloaded files:"
echo "   # Script automatically offers to start build after copying"
echo ""
echo "💡 Why this workflow:"
echo "   ✅ Develop in ~/micapipe (comfortable environment)"
echo "   ✅ Build in /host/... (unlimited space + pre-downloaded files)"
echo "   ✅ No space issues during Docker build"
echo "   ✅ Uses existing FreeSurfer, FSL, AFNI downloads"
echo ""
echo "🎯 Expected results:"
echo "   - Copy: < 1 minute"
echo "   - Base build: 45-90 minutes (one-time)"  
echo "   - Future builds: 3-5 minutes"
echo ""
echo "🔧 SERVER commands (copy-paste):"
echo ""

cat << 'EOF'
# === RUN ON SERVER ===

cd ~/micapipe
git pull origin comprehensive-base-image
./migrate_comprehensive_base_to_server.sh

# === DONE! ===
# Script handles copying and offers to build automatically
EOF

echo ""
echo "✅ Server-side 3-step process!"
echo "💡 migrate_comprehensive_base_to_server.sh does all the heavy lifting"