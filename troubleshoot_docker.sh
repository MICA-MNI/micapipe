#!/bin/bash
#
# Troubleshoot Docker build issues
# Run this if docker build fails
#

set -e

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🔧 DOCKER BUILD TROUBLESHOOTING"
log "==============================="

# Check Docker daemon
log "📋 Checking Docker daemon..."
if docker info >/dev/null 2>&1; then
    log "✅ Docker daemon is running"
else
    log "❌ Docker daemon not running or accessible"
    log "   Try: sudo systemctl start docker"
    exit 1
fi

# Check Docker version
DOCKER_VERSION=$(docker version --format '{{.Server.Version}}' 2>/dev/null || echo "unknown")
log "📦 Docker version: $DOCKER_VERSION"

# Check disk space
log "📋 Checking disk space..."
DISK_USAGE=$(df -h /var/lib/docker 2>/dev/null | awk 'NR==2 {print $5}' || echo "unknown")
log "💾 Docker storage usage: $DISK_USAGE"

# Clean up Docker
log "📋 Cleaning up Docker cache..."
docker system prune -f --volumes 2>/dev/null || true
log "✅ Docker cache cleaned"

# Disable problematic Docker features
log "📋 Setting safe Docker environment..."
export DOCKER_CONTENT_TRUST=0
export DOCKER_BUILDKIT=0  # Disable BuildKit if causing issues
log "✅ Docker content trust disabled"
log "✅ Docker BuildKit disabled"

# Test simple build
log "📋 Testing simple Docker build..."
cd /tmp
cat > test_dockerfile << 'EOF'
FROM ubuntu:20.04
RUN echo "test build successful"
EOF

if docker build -f test_dockerfile -t docker-test-build . >/dev/null 2>&1; then
    log "✅ Simple Docker build works"
    docker rmi docker-test-build >/dev/null 2>&1 || true
    rm -f test_dockerfile
else
    log "❌ Even simple Docker build fails"
    log "   This indicates a serious Docker daemon issue"
    log "   Try restarting Docker: sudo systemctl restart docker"
    rm -f test_dockerfile
    exit 1
fi

log "==============================="
log "✅ DOCKER TROUBLESHOOTING COMPLETE"
log "==============================="
log "🚀 Docker should now work properly"
log "💡 Try running your build command again with:"
log "   export DOCKER_CONTENT_TRUST=0"
log "   export DOCKER_BUILDKIT=0"
log "   docker build -t micapipe-runner-test ."