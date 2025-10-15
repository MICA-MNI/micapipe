#!/bin/bash
# =============================================================================
# MICAPIPE SIF WORKAROUND SCRIPT
# =============================================================================
# This script allows you to use the current SIF file (micapipe_v1_beta.sif)
# even with the FreeSurfer setup error, by bypassing the startup script issues.
#
# BACKGROUND:
# The current SIF has a FreeSurfer setup error in /neurodocker/startup.sh:
#   "/opt/freesurfer-7.4.1/SetUpFreeSurfer.sh: No such file or directory"
#
# WORKAROUND:
# Instead of using the default entrypoint, we run commands directly in the
# container with proper environment setup.
#
# USAGE:
#   ./run_micapipe_workaround.sh -h                    # Show help
#   ./run_micapipe_workaround.sh --version             # Show version
#   ./run_micapipe_workaround.sh [other micapipe args] # Run micapipe commands
#
# REQUIREMENTS:
#   - SIF file at: /data_/mica1/01_programs/singularity/micapipe_v1_beta.sif
#   - Singularity installed
# =============================================================================

set -euo pipefail

# Configuration
SIF_PATH="/data_/mica1/01_programs/singularity/micapipe_v1_beta.sif"
SCRIPT_NAME=$(basename "$0")

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Functions
print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

show_usage() {
    cat << EOF
🔧 MICApipe SIF Workaround Script

DESCRIPTION:
    This script runs micapipe commands using the current SIF file, bypassing
    the FreeSurfer setup error by manually setting up the environment.

USAGE:
    $SCRIPT_NAME [micapipe arguments]

EXAMPLES:
    $SCRIPT_NAME -h                    # Show micapipe help
    $SCRIPT_NAME --version             # Show version info
    $SCRIPT_NAME --help                # Show micapipe help
    
TROUBLESHOOTING:
    If you see FreeSurfer-related errors, this is expected with the current SIF.
    For full functionality, rebuild the base Docker image with the latest fixes.

SIF FILE:
    $SIF_PATH

EOF
}

# Check if SIF file exists
if [[ ! -f "$SIF_PATH" ]]; then
    print_error "SIF file not found: $SIF_PATH"
    echo ""
    echo "Expected location: $SIF_PATH"
    echo "Please ensure the SIF file exists or update SIF_PATH in this script."
    exit 1
fi

# Show usage if no arguments or help requested
if [[ $# -eq 0 ]] || [[ "$1" == "--usage" ]] || [[ "$1" == "--script-help" ]]; then
    show_usage
    exit 0
fi

print_info "Using SIF: $(basename "$SIF_PATH")"
print_warning "Note: This SIF has FreeSurfer setup issues but basic commands work"

# Set up environment variables for FreeSurfer (workaround)
export FREESURFER_HOME="/opt/freesurfer-7.4.1"
export SUBJECTS_DIR="$FREESURFER_HOME/subjects"
export PATH="/opt/freesurfer-7.4.1/bin:$PATH"

# Create a temporary script to run inside the container
TEMP_SCRIPT=$(mktemp)
cat > "$TEMP_SCRIPT" << 'EOF'
#!/bin/bash
# Workaround script that runs inside the container

# Set up minimal FreeSurfer environment (bypass startup script issues)
export FREESURFER_HOME="/opt/freesurfer-7.4.1"
export SUBJECTS_DIR="$FREESURFER_HOME/subjects"
export FUNCTIONALS_DIR="$FREESURFER_HOME/sessions"  
export PATH="/opt/freesurfer-7.4.1/bin:$PATH"

# Set up other tool paths (from functions/init.sh Docker section)
export AFNIDIR="/opt/afni-23.1.09"
export ANTSPATH="/opt/ants-2.3.4/bin"
export workbench_path="/opt/workbench-1.4.2/bin"
export FSLDIR="/opt/fsl-6.0.5.1"
export FSL_DIR="/opt/fsl-6.0.5.1"
export FSL_BIN="${FSLDIR}/bin"
export mrtrixDir="/opt/mrtrix3-3.0.1"

# Update PATH with all tools
export PATH="${AFNIDIR}:${ANTSPATH}:${workbench_path}:${FREESURFER_HOME}/bin:${mrtrixDir}/bin:${mrtrixDir}/lib:${FSLDIR}:${FSL_BIN}:${PATH}"

# Run the micapipe command with passed arguments
exec /opt/micapipe/micapipe "$@"
EOF

chmod +x "$TEMP_SCRIPT"

print_info "Running micapipe with arguments: $*"
print_info "Environment: FreeSurfer, FSL, ANTs, MRtrix3 paths set manually"

# Run singularity with the workaround script
if singularity exec \
    --bind "$TEMP_SCRIPT:/tmp/run_micapipe.sh" \
    "$SIF_PATH" \
    /tmp/run_micapipe.sh "$@"; then
    
    print_success "Command completed successfully"
else
    EXIT_CODE=$?
    print_error "Command failed with exit code: $EXIT_CODE"
    
    if [[ $EXIT_CODE -eq 127 ]]; then
        print_warning "Command not found - check micapipe installation in SIF"
    elif [[ "$*" == *"-h"* ]] || [[ "$*" == *"--help"* ]]; then
        print_info "Help command may have succeeded despite error code"
    fi
fi

# Clean up
rm -f "$TEMP_SCRIPT"

print_info "Workaround script completed"
echo ""
print_warning "For full functionality without workarounds:"
echo "  1. Rebuild base Docker image with FreeSurfer fixes"
echo "  2. Create new SIF file from fixed Docker image"
echo "  3. Use standard singularity run commands"