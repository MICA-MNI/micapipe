#!/bin/bash
#
# Check progress of running Singularity build - DON'T INTERRUPT THE BUILD!
#

BASE_DIR="/host/cassio/export03/data/enning"
OUTPUT_PATH="${BASE_DIR}/singularity/micapipe_v1_beta.sif"
CACHE_DIR="${BASE_DIR}/.singularity_cache"
TEMP_DIR="${BASE_DIR}/.singularity_tmp"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🔍 CHECKING BUILD PROGRESS (non-intrusive)"
log "================================================"

# Check for running processes
log "🔄 Running processes:"
SING_PROCS=$(ps aux | grep singularity | grep -v grep | wc -l)
DOCKER_PROCS=$(ps aux | grep "docker save" | grep -v grep | wc -l)

if [ "$SING_PROCS" -gt 0 ]; then
    log "✅ Found $SING_PROCS singularity processes (BUILD IS RUNNING):"
    ps aux | grep singularity | grep -v grep | awk '{print "   PID:"$2" CPU:"$3"% MEM:"$4"% TIME:"$10" CMD:"$11" "$12" "$13}'
else
    log "❌ No singularity processes found"
fi

if [ "$DOCKER_PROCS" -gt 0 ]; then
    log "✅ Found $DOCKER_PROCS docker save processes:"
    ps aux | grep "docker save" | grep -v grep | awk '{print "   PID:"$2" CPU:"$3"% MEM:"$4"% TIME:"$10}'
else
    log "ℹ️  No docker save processes (normal if using streaming)"
fi

log ""
log "📁 File sizes:"

# Check output file
if [ -f "$OUTPUT_PATH" ]; then
    OUTPUT_SIZE=$(du -h "$OUTPUT_PATH" | cut -f1)
    OUTPUT_BYTES=$(stat -f%z "$OUTPUT_PATH" 2>/dev/null || echo 0)
    log "✅ Output SIF: $OUTPUT_SIZE ($OUTPUT_BYTES bytes)"
    
    # Check if file is growing
    sleep 5
    NEW_BYTES=$(stat -f%z "$OUTPUT_PATH" 2>/dev/null || echo 0)
    if [ "$NEW_BYTES" -gt "$OUTPUT_BYTES" ]; then
        GROWTH=$((NEW_BYTES - OUTPUT_BYTES))
        log "📈 File is GROWING! +$GROWTH bytes in 5 seconds"
    else
        log "⚠️  File size unchanged (might be stuck or in different phase)"
    fi
else
    log "❌ No output file yet at: $OUTPUT_PATH"
fi

# Check temp directory
if [ -d "$TEMP_DIR" ]; then
    TEMP_SIZE=$(du -sh "$TEMP_DIR" 2>/dev/null | cut -f1 || echo "0")
    TEMP_FILES=$(find "$TEMP_DIR" -type f 2>/dev/null | wc -l)
    log "💾 Temp directory: $TEMP_SIZE ($TEMP_FILES files)"
else
    log "ℹ️  No temp directory"
fi

# Check cache directory  
if [ -d "$CACHE_DIR" ]; then
    CACHE_SIZE=$(du -sh "$CACHE_DIR" 2>/dev/null | cut -f1 || echo "0")
    log "🗄️  Cache directory: $CACHE_SIZE"
else
    log "ℹ️  No cache directory"
fi

# Check disk space
log ""
log "💽 Disk space:"
AVAILABLE=$(df -h "$BASE_DIR" | awk 'NR==2 {print $4}')
USED_PERCENT=$(df -h "$BASE_DIR" | awk 'NR==2 {print $5}')
log "📊 Available: $AVAILABLE (${USED_PERCENT} used)"

# Check for any tar files (backup method)
TAR_FILES=$(find "$BASE_DIR" -name "micapipe_docker_*.tar" 2>/dev/null)
if [ -n "$TAR_FILES" ]; then
    log ""
    log "📦 Found tar files (backup method active):"
    echo "$TAR_FILES" | while read tar_file; do
        if [ -f "$tar_file" ]; then
            TAR_SIZE=$(du -h "$tar_file" | cut -f1)
            log "   $tar_file: $TAR_SIZE"
        fi
    done
fi

log ""
log "🕐 Quick status check:"
if [ "$SING_PROCS" -gt 0 ] || [ "$DOCKER_PROCS" -gt 0 ]; then
    log "✅ BUILD IS RUNNING - be patient!"
    log "⏰ Typical build time: 15-30 minutes for streaming, 45-60 minutes for tar method"
    log "🔄 Run this script again in 10 minutes to check progress"
else
    log "❌ NO BUILD PROCESSES FOUND"
    log "💡 Either build completed, failed, or got stuck"
    log "🔍 Check build logs or restart the build"
fi

log ""
log "================================================"