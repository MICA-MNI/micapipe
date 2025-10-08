#!/bin/bash
#
# Debug Singularity build issues
#
# This script helps diagnose why Singularity build might be stuck

echo "============================================="
echo "🔍 SINGULARITY BUILD DEBUGGER"
echo "============================================="
echo ""

# Check if singularity process is running
echo "1️⃣  Checking for running Singularity processes..."
SING_PROCS=$(ps aux | grep -i singularity | grep -v grep)
if [ -z "$SING_PROCS" ]; then
    echo "❌ No singularity processes found"
else
    echo "✅ Singularity processes running:"
    echo "$SING_PROCS"
fi
echo ""

# Check Docker daemon
echo "2️⃣  Checking Docker daemon..."
if docker info &> /dev/null; then
    echo "✅ Docker daemon is running"
    echo "   Images available:"
    docker images | grep micapipe | head -5
else
    echo "❌ Docker daemon not accessible"
fi
echo ""

# Check temp directory space
echo "3️⃣  Checking temp directory space..."
TEMP_DIRS=(
    "/tmp"
    "$TMPDIR"
    "$SINGULARITY_TMPDIR"
    "/var/tmp"
    "/host/cassio/export03/data/enning"
)

for dir in "${TEMP_DIRS[@]}"; do
    if [ -n "$dir" ] && [ -d "$dir" ]; then
        SPACE=$(df -h "$dir" | awk 'NR==2 {print $4}')
        echo "   $dir: $SPACE available"
    fi
done
echo ""

# Check output directory permissions
echo "4️⃣  Checking output directory permissions..."
OUTPUT_DIR="/host/cassio/export03/data/enning/singularity"
if [ -d "$OUTPUT_DIR" ]; then
    echo "✅ Directory exists: $OUTPUT_DIR"
    ls -ld "$OUTPUT_DIR"
    
    # Check if writable
    if touch "$OUTPUT_DIR/test_write_$$.tmp" 2>/dev/null; then
        rm -f "$OUTPUT_DIR/test_write_$$.tmp"
        echo "✅ Directory is writable"
    else
        echo "❌ Directory is NOT writable"
    fi
else
    echo "❌ Directory does not exist: $OUTPUT_DIR"
fi
echo ""

# Check for partial files
echo "5️⃣  Checking for partial/temporary files..."
if [ -d "$OUTPUT_DIR" ]; then
    PARTIAL_FILES=$(find "$OUTPUT_DIR" -name "*.sif*" -o -name "*tmp*" 2>/dev/null)
    if [ -z "$PARTIAL_FILES" ]; then
        echo "   No partial files found"
    else
        echo "   Partial files:"
        find "$OUTPUT_DIR" -name "*.sif*" -o -name "*tmp*" 2>/dev/null | xargs ls -lh
    fi
fi
echo ""

# Check Singularity cache
echo "6️⃣  Checking Singularity cache..."
CACHE_DIRS=(
    "$HOME/.singularity/cache"
    "$HOME/.apptainer/cache"
    "/tmp/singularity-cache-$(whoami)"
)

for cache in "${CACHE_DIRS[@]}"; do
    if [ -d "$cache" ]; then
        SIZE=$(du -sh "$cache" 2>/dev/null | cut -f1)
        echo "   $cache: $SIZE"
    fi
done
echo ""

# Check for file locks
echo "7️⃣  Checking for file locks..."
if command -v lsof &> /dev/null; then
    LOCKS=$(lsof 2>/dev/null | grep -i "singularity\|docker" | grep -i "micapipe" || echo "None")
    if [ "$LOCKS" = "None" ]; then
        echo "   No locks found"
    else
        echo "   Active locks:"
        echo "$LOCKS" | head -10
    fi
else
    echo "   (lsof not available)"
fi
echo ""

# Check system resources
echo "8️⃣  Checking system resources..."
echo "   Memory:"
free -h | grep -E "Mem:|Swap:"
echo ""
echo "   Load average:"
uptime
echo ""

echo "============================================="
echo "💡 RECOMMENDATIONS"
echo "============================================="
echo ""

# Provide recommendations
if [ -z "$SING_PROCS" ]; then
    echo "⚠️  No Singularity process running - build may have failed silently"
    echo "   Try running: ./build_singularity.sh 2>&1 | tee singularity_build.log"
fi

if ! docker info &> /dev/null; then
    echo "⚠️  Docker daemon issue detected"
    echo "   Try: sudo systemctl restart docker"
fi

echo ""
echo "🔧 Troubleshooting steps:"
echo "   1. Check if build is actually running:"
echo "      ps aux | grep singularity"
echo ""
echo "   2. Monitor system resources:"
echo "      watch -n 2 'df -h /host/cassio/export03 && free -h'"
echo ""
echo "   3. Try with verbose output:"
echo "      singularity build --debug /path/to/output.sif docker-daemon://micapipe:latest"
echo ""
echo "   4. Clear Singularity cache and retry:"
echo "      rm -rf ~/.singularity/cache ~/.apptainer/cache"
echo "      ./build_singularity.sh"
echo ""
echo "   5. Check for 'nodev' mount issues:"
echo "      mount | grep nodev"
echo "      # If /host/cassio/export03 has nodev, try output to /tmp first:"
echo "      SINGULARITY_DIR=/tmp ./build_singularity.sh"
echo ""
