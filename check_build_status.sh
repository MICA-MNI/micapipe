#!/bin/bash
# ============================================================================
# BUILD STATUS CHECKER
# Quick status check for ongoing Docker builds
# ============================================================================

echo "🔍 MICApipe Build Status Check"
echo "=============================="

# Check if we're in the right location
if [[ "$PWD" != *"/host/cassio/export03/data/enning/downloads"* ]]; then
    echo "⚠️  Not in expected build directory"
    echo "   Current: $PWD"
    echo "   Expected: /host/cassio/export03/data/enning/downloads"
fi

# Check for ongoing Docker builds
echo ""
echo "🐳 Docker Build Status:"
DOCKER_PROCESSES=$(docker ps -a --format "table {{.Image}}\t{{.Status}}\t{{.Names}}" | grep -E "(micapipe|ubuntu)" 2>/dev/null || echo "No MICApipe containers found")
echo "$DOCKER_PROCESSES"

# Check recent build logs
echo ""
echo "📋 Recent Build Logs:"
ls -lt build_comprehensive_base_*.log 2>/dev/null | head -3 || echo "No build logs found"

# Check last few lines of most recent log
LATEST_LOG=$(ls -t build_comprehensive_base_*.log 2>/dev/null | head -1)
if [[ -f "$LATEST_LOG" ]]; then
    echo ""
    echo "📄 Last 10 lines of $LATEST_LOG:"
    echo "─────────────────────────────────────────────"
    tail -10 "$LATEST_LOG"
    echo "─────────────────────────────────────────────"
    
    # Check for common progress indicators
    echo ""
    echo "🎯 Build Progress Indicators:"
    STEPS_COMPLETED=$(grep -c "Step [0-9]*/[0-9]*" "$LATEST_LOG" 2>/dev/null || echo "0")
    echo "   Steps completed: $STEPS_COMPLETED"
    
    if grep -q "Successfully built" "$LATEST_LOG" 2>/dev/null; then
        echo "   ✅ Build completed successfully!"
    elif grep -q "returned a non-zero code" "$LATEST_LOG" 2>/dev/null; then
        echo "   ❌ Build failed - check log for details"
    elif grep -q "Running in" "$LATEST_LOG" 2>/dev/null; then
        echo "   🏗️  Build in progress..."
    else
        echo "   📝 Build status unclear - check log"
    fi
else
    echo "No build logs found"
fi

# Check disk space
echo ""
echo "💾 Disk Space:"
df -h /host/cassio/export03/data/enning | tail -1

echo ""
echo "💡 To monitor live: tail -f $LATEST_LOG"