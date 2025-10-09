#!/bin/bash
#
# Emergency: Kill stuck Singularity build and clean up
# USE ONLY IF BUILD IS CONFIRMED STUCK!
#

BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_PATH="${BASE_DIR}/singularity/micapipe_v1_beta.sif"
CACHE_DIR="${BASE_DIR}/.singularity_cache"
TEMP_DIR="${BASE_DIR}/.singularity_tmp"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🛑 EMERGENCY BUILD KILLER"
log "=========================="

# Show what will be killed
log "🔍 Current processes to kill:"
ps aux | grep -E "(singularity|docker save)" | grep -v grep | awk '{print "   PID:"$2" CMD:"$11" "$12" "$13}'

read -p "⚠️  Are you sure you want to KILL the build? (yes/no): " confirm
if [ "$confirm" != "yes" ]; then
    log "❌ Cancelled. Build continues running."
    exit 0
fi

log "🛑 Killing processes..."

# Kill singularity processes
KILLED_SING=$(pkill -f singularity 2>/dev/null && echo "killed" || echo "none")
log "🔫 Singularity processes: $KILLED_SING"

# Kill docker save processes  
KILLED_DOCKER=$(pkill -f "docker save" 2>/dev/null && echo "killed" || echo "none")
log "🔫 Docker save processes: $KILLED_DOCKER"

sleep 2

log "🧹 Cleaning up files..."

# Remove partial output
if [ -f "$OUTPUT_PATH" ]; then
    rm -f "$OUTPUT_PATH"
    log "🗑️  Removed partial SIF: $OUTPUT_PATH"
fi

# Clean temp directory
if [ -d "$TEMP_DIR" ]; then
    TEMP_SIZE_BEFORE=$(du -sh "$TEMP_DIR" 2>/dev/null | cut -f1 || echo "0")
    rm -rf "$TEMP_DIR"/* 2>/dev/null || true
    log "🗑️  Cleaned temp dir: was $TEMP_SIZE_BEFORE"
fi

# Remove any tar files
find "$BASE_DIR" -name "micapipe_docker_*.tar" -exec rm -f {} \; 2>/dev/null
log "🗑️  Removed any tar files"

log ""
log "✅ CLEANUP COMPLETE"
log "💡 You can now restart the build with: ./build_singularity.sh"
log ""