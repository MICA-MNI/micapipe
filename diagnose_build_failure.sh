#!/bin/bash

# Build Failure Diagnostic Script
# Analyzes build logs to identify the exact failure point

echo "🔍 MICApipe Build Failure Diagnostic"
echo "===================================="
echo ""

LOG_DIR="./build_logs"

if [ ! -d "$LOG_DIR" ]; then
    echo "❌ Build logs directory not found: $LOG_DIR"
    echo "   Make sure you're in the micapipe directory where you ran the build"
    exit 1
fi

echo "📊 Build Log Analysis:"
echo "====================="

# Check main build log
if [ -f "$LOG_DIR/main_build.log" ]; then
    echo ""
    echo "🔍 Main Build Log Analysis:"
    echo "Size: $(du -h "$LOG_DIR/main_build.log" | cut -f1)"
    
    # Look for common error patterns
    echo ""
    echo "❌ ERROR SUMMARY:"
    echo "=================="
    
    # Check for download failures
    if grep -i "failed\|error\|timeout\|connection" "$LOG_DIR/main_build.log" | tail -10; then
        echo ""
    else
        echo "No obvious errors found in grep search"
    fi
    
    echo ""
    echo "📝 LAST 30 LINES OF BUILD LOG:"
    echo "=============================="
    tail -30 "$LOG_DIR/main_build.log"
    
    echo ""
    echo "🎯 LIKELY FAILURE POINT:"
    echo "========================"
    # Look for the step where it failed
    grep -n "Step\|RUN\|COPY\|FROM" "$LOG_DIR/main_build.log" | tail -5
    
else
    echo "❌ Main build log not found: $LOG_DIR/main_build.log"
fi

# Check session log for build configuration
if [ -f "$LOG_DIR/build_session.log" ]; then
    echo ""
    echo "⚙️  BUILD CONFIGURATION:"
    echo "======================="
    grep -E "CUDA|GPU|build|Starting" "$LOG_DIR/build_session.log" | tail -10
fi

echo ""
echo "💡 COMMON SOLUTIONS:"
echo "==================="
echo "1. NETWORK ISSUES:"
echo "   - If download failures: try ./test_fsl_download.sh"
echo "   - Retry during off-peak hours"
echo "   - Check proxy/firewall settings"
echo ""
echo "2. RESOURCE ISSUES:"
echo "   - Disk space: df -h"
echo "   - Memory: free -h"
echo "   - Docker system: docker system df"
echo ""
echo "3. DEPENDENCY ISSUES:"
echo "   - Try: docker system prune -a (clean up)"
echo "   - Restart Docker daemon"
echo "   - Check Docker version compatibility"
echo ""
echo "4. BUILD CACHE ISSUES:"
echo "   - Try building without --no-cache"
echo "   - Or clean build: docker build --no-cache"
echo ""
echo "🔧 NEXT STEPS:"
echo "============="
echo "1. Review the error details above"
echo "2. Try the suggested solutions"
echo "3. If network-related, run: ./test_fsl_download.sh"
echo "4. For quick retry: ./build_container.sh"
echo "5. For detailed retry: ./server_build_test.sh"
echo ""
echo "📞 Need help? Check the specific error message in the 'LAST 30 LINES' section above."