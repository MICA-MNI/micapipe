#!/bin/bash
# =============================================================================
# QUICK MICAPIPE SIF TESTER (SERVER VERSION)
# =============================================================================
# Simple script to test the current SIF file with FreeSurfer setup workaround
# Bypasses /neurodocker/startup.sh issues by setting up environment manually
# =============================================================================

SIF_PATH="/host/cassio/export03/data/enning/singularity/micapipe_v1_beta.sif"

echo "🧪 Testing MICApipe SIF with FreeSurfer workaround"
echo "📁 SIF: $SIF_PATH"
echo ""

if [[ ! -f "$SIF_PATH" ]]; then
    echo "❌ SIF file not found: $SIF_PATH"
    exit 1
fi

echo "🔧 Running micapipe -h with environment workaround..."
echo ""

# Use singularity exec with manual environment setup instead of the broken startup script
singularity exec "$SIF_PATH" bash -c '
    # Set up FreeSurfer environment manually (bypass broken startup script)
    export FREESURFER_HOME="/opt/freesurfer-7.4.1"
    export SUBJECTS_DIR="$FREESURFER_HOME/subjects"
    export PATH="/opt/freesurfer-7.4.1/bin:$PATH"
    
    # Set up other tools
    export AFNIDIR="/opt/afni-23.1.09"
    export ANTSPATH="/opt/ants-2.3.4/bin" 
    export FSLDIR="/opt/fsl-6.0.5.1"
    export PATH="${AFNIDIR}:${ANTSPATH}:${FREESURFER_HOME}/bin:${FSLDIR}/bin:${PATH}"
    
    # Run micapipe help
    echo "🚀 Running: /opt/micapipe/micapipe -h"
    echo ""
    /opt/micapipe/micapipe -h
'

echo ""
echo "✅ Test completed"
echo ""
echo "💡 This workaround bypasses the FreeSurfer startup script error."
echo "   For full functionality, rebuild the base Docker image with fixes."