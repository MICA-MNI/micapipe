#!/bin/bash
# ============================================================================
# PRE-COMPUTE CONDA ENVIRONMENT DEPENDENCY GRAPH
# ============================================================================
# This script pre-solves the conda environment and creates an explicit 
# specification file that can be installed instantly without solving.
#
# RUN THIS ONCE to generate the pre-solved environment spec:
#   ./precompute_conda_environment.sh
#
# The output file (micapipe_environment_explicit.txt) will contain:
#   - Exact package versions
#   - Direct download URLs
#   - Pre-computed dependency graph
#
# This makes Docker builds INSTANT for the conda environment step!
# ============================================================================

set -euo pipefail

echo "🧮 Pre-computing conda environment dependency graph..."
echo "========================================================"
echo ""

# Check if conda/mamba is available
if ! command -v mamba &> /dev/null; then
    echo "❌ mamba not found. Please install conda/mamba first."
    echo ""
    echo "Installation:"
    echo "  wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
    echo "  bash Mambaforge-Linux-x86_64.sh -b -p ~/mambaforge"
    echo "  ~/mambaforge/bin/conda init"
    exit 1
fi

# Create temporary environment to solve dependencies
TEMP_ENV="micapipe_solve_tmp"
echo "📦 Creating temporary environment to solve dependencies..."

# Remove if exists
mamba env remove -n "$TEMP_ENV" -y 2>/dev/null || true

# Create environment from yml file (this does the solving)
echo "⏱️  Solving dependencies (this may take 5-10 minutes)..."
echo ""
if mamba env create -n "$TEMP_ENV" -f micapipe_environment.yml; then
    echo ""
    echo "✅ Dependencies solved successfully!"
    echo ""
    
    # Export explicit specification
    echo "💾 Exporting explicit specification..."
    mamba list -n "$TEMP_ENV" --explicit > micapipe_environment_explicit.txt
    
    # Also export traditional spec for reference
    mamba env export -n "$TEMP_ENV" > micapipe_environment_solved.yml
    
    # Clean up temporary environment
    echo "🧹 Cleaning up..."
    mamba env remove -n "$TEMP_ENV" -y
    
    echo ""
    echo "🎉 SUCCESS! Pre-computed environment files created:"
    echo "   📄 micapipe_environment_explicit.txt  (use this in Dockerfile)"
    echo "   📄 micapipe_environment_solved.yml    (reference only)"
    echo ""
    echo "📊 File sizes:"
    ls -lh micapipe_environment_explicit.txt micapipe_environment_solved.yml
    echo ""
    echo "🚀 Now update Dockerfile.base to use the explicit spec:"
    echo "   COPY micapipe_environment_explicit.txt /tmp/"
    echo "   RUN conda create -n micapipe --file /tmp/micapipe_environment_explicit.txt"
    echo ""
    echo "⚡ This will make builds INSTANT (no solving needed)!"
    
else
    echo ""
    echo "❌ Failed to solve environment!"
    echo "Check micapipe_environment.yml for conflicts"
    exit 1
fi
