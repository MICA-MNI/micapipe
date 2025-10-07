#!/bin/bash
# Remote build script - run from local Mac
# This will SSH to venice and start the Docker build

set -e

echo "🚀 Starting remote Docker build on venice..."
echo ""
echo "This will:"
echo "  1. Pull latest code from comprehensive-base-image branch"
echo "  2. Migrate files to server build location"
echo "  3. Start Docker base image build (45-90 minutes)"
echo ""
read -p "Press Enter to continue or Ctrl+C to cancel..."

# SSH to venice via login node and execute build commands
ssh -t eyang@login.bic.mni.mcgill.ca "ssh -t venice '\
    cd ~/micapipe && \
    echo \"📥 Pulling latest changes...\" && \
    git pull origin comprehensive-base-image && \
    echo \"\" && \
    echo \"📦 Migrating files to server...\" && \
    ./migrate_comprehensive_base_to_server.sh && \
    echo \"\" && \
    echo \"🔨 Starting Docker build...\" && \
    cd /host/cassio/export03/data/enning/downloads && \
    ./build_base_image_server.sh
'"

echo ""
echo "✅ Build command executed!"
echo "Note: Build will continue on server even if you disconnect."
