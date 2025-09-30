#!/bin/bash

# Quick FSL Download Test Script
# Tests FSL download outside Docker to diagnose network issues

set -e

echo "🔍 Testing FSL Download Connectivity"
echo "===================================="

# Test URLs
FSL_PRIMARY="https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.2-centos6_64.tar.gz"
FSL_MIRROR="https://www.nitrc.org/frs/download.php/11344/fsl-6.0.2-centos6_64.tar.gz"

echo "📡 Testing network connectivity..."

# Test basic connectivity
if ping -c 3 fsl.fmrib.ox.ac.uk > /dev/null 2>&1; then
    echo "✅ fsl.fmrib.ox.ac.uk is reachable"
else
    echo "❌ fsl.fmrib.ox.ac.uk is NOT reachable"
fi

if ping -c 3 www.nitrc.org > /dev/null 2>&1; then
    echo "✅ www.nitrc.org is reachable"
else
    echo "❌ www.nitrc.org is NOT reachable"
fi

echo ""
echo "🌐 Testing HTTP connectivity..."

# Test HTTP HEAD requests
echo "Testing primary FSL URL..."
if curl -I --connect-timeout 30 --max-time 60 "$FSL_PRIMARY" 2>/dev/null | head -1; then
    echo "✅ Primary FSL URL responds"
else
    echo "❌ Primary FSL URL failed"
fi

echo "Testing mirror URL..."
if curl -I --connect-timeout 30 --max-time 60 "$FSL_MIRROR" 2>/dev/null | head -1; then
    echo "✅ Mirror URL responds"
else
    echo "❌ Mirror URL failed"
fi

echo ""
echo "📥 Testing actual download (first 1MB)..."

# Test actual download
mkdir -p /tmp/fsl_test
cd /tmp/fsl_test

echo "Trying primary URL..."
if timeout 60 curl -fsSL --retry 3 --retry-delay 5 \
   -r 0-1048576 "$FSL_PRIMARY" -o fsl_test_primary.tar.gz 2>/dev/null; then
    echo "✅ Primary URL download test successful ($(du -h fsl_test_primary.tar.gz | cut -f1))"
else
    echo "❌ Primary URL download test failed"
fi

echo "Trying mirror URL..."
if timeout 60 curl -fsSL --retry 3 --retry-delay 5 \
   -r 0-1048576 "$FSL_MIRROR" -o fsl_test_mirror.tar.gz 2>/dev/null; then
    echo "✅ Mirror URL download test successful ($(du -h fsl_test_mirror.tar.gz | cut -f1))"
else
    echo "❌ Mirror URL download test failed"
fi

# Cleanup
cd - > /dev/null
rm -rf /tmp/fsl_test

echo ""
echo "🔧 Recommendations:"
echo "==================="

# Check if we're behind a proxy or firewall
if [ -n "$HTTP_PROXY" ] || [ -n "$HTTPS_PROXY" ] || [ -n "$http_proxy" ] || [ -n "$https_proxy" ]; then
    echo "📋 Proxy detected - this may cause Docker build issues"
    echo "   Consider adding proxy settings to Docker build"
fi

# Check Docker daemon proxy settings
if docker info 2>/dev/null | grep -q "HTTP Proxy\|HTTPS Proxy"; then
    echo "📋 Docker has proxy settings configured"
else
    echo "⚠️  Docker may need proxy configuration for downloads"
fi

echo ""
echo "💡 Solutions to try:"
echo "1. Increase Docker build timeout: docker build --build-arg BUILDKIT_TIMEOUT=3600"
echo "2. Use cached build if possible: docker build (without --no-cache)"
echo "3. Build during off-peak hours when network is less congested"
echo "4. Use the robust FSL download in updated Dockerfile (already implemented)"
echo ""
echo "🚀 The Dockerfile has been updated with:"
echo "   - Multiple retry attempts with longer timeouts"
echo "   - Alternative download methods (curl + wget)"
echo "   - Mirror URL fallback"
echo "   - Graceful degradation if FSL download fails"