#!/bin/bash
# ============================================================================
# GENERATE CONDA LOCK FILE (CROSS-PLATFORM REPRODUCIBLE)
# ============================================================================
# This uses conda-lock to create a fully reproducible lock file
# that works across different systems and guarantees exact versions.
#
# USAGE:
#   ./generate_conda_lock.sh
#
# This creates: conda-linux-64.lock (use in Dockerfile)
# ============================================================================

set -euo pipefail

echo "🔒 Generating conda lock file for linux-64..."
echo "=============================================="
echo ""

# Check if conda-lock is installed
if ! command -v conda-lock &> /dev/null; then
    echo "📦 Installing conda-lock..."
    pip install conda-lock
    echo ""
fi

# Generate lock file for linux-64 platform
echo "⏱️  Solving dependencies for linux-64 (this may take 5-10 minutes)..."
echo ""

if conda-lock lock \
    --file micapipe_environment.yml \
    --platform linux-64 \
    --kind explicit \
    --filename-template "micapipe-{platform}.lock"; then
    
    echo ""
    echo "✅ Lock file generated successfully!"
    echo ""
    echo "📄 Created: micapipe-linux-64.lock"
    ls -lh micapipe-linux-64.lock
    echo ""
    echo "🚀 This lock file ensures:"
    echo "   ✅ Exact package versions"
    echo "   ✅ Cross-platform reproducibility"
    echo "   ✅ Fast installation (no solving)"
    echo "   ✅ Deterministic builds"
    echo ""
    echo "📝 Use in Dockerfile:"
    echo "   COPY micapipe-linux-64.lock /tmp/"
    echo "   RUN conda create -n micapipe --file /tmp/micapipe-linux-64.lock"
    echo ""
    
else
    echo ""
    echo "❌ Failed to generate lock file!"
    echo "Check micapipe_environment.yml for conflicts"
    exit 1
fi

# Optional: Generate lock files for multiple platforms
read -p "Generate lock files for other platforms? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo ""
    echo "🔒 Generating lock files for all platforms..."
    
    conda-lock lock \
        --file micapipe_environment.yml \
        --platform linux-64 \
        --platform osx-64 \
        --platform osx-arm64 \
        --kind explicit \
        --filename-template "micapipe-{platform}.lock"
    
    echo ""
    echo "✅ Multi-platform lock files generated!"
    ls -lh micapipe-*.lock
fi
