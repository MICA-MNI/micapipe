#!/bin/bash
# ============================================================================
# DOCKER SPACE CLEANUP SCRIPT (NO SUDO REQUIRED)
# ============================================================================
# Run this BEFORE building Docker images to free up space
# Usage: ./cleanup_docker_space.sh
#
# IMPORTANT: This script works WITHOUT sudo access
# - Cleans Docker resources you own
# - Skips operations requiring root privileges
# - Safe to run on shared servers like venice/cassio
# ============================================================================

set -e

echo "🧹 Docker Space Cleanup Script (No Sudo Required)"
echo "========================================"
echo ""

# Show initial disk space
echo "📊 Initial disk usage:"
df -h /var/lib/docker 2>/dev/null || df -h /
echo ""

# Step 1: Remove stopped containers
echo "🗑️  Step 1: Removing stopped containers..."
STOPPED=$(docker ps -aq -f status=exited -f status=created 2>/dev/null | wc -l)
if [ "$STOPPED" -gt 0 ]; then
    docker rm $(docker ps -aq -f status=exited -f status=created) 2>/dev/null || true
    echo "   ✅ Removed $STOPPED stopped containers"
else
    echo "   ℹ️  No stopped containers to remove"
fi
echo ""

# Step 2: Remove dangling images
echo "🗑️  Step 2: Removing dangling images..."
DANGLING=$(docker images -qf dangling=true 2>/dev/null | wc -l)
if [ "$DANGLING" -gt 0 ]; then
    docker rmi $(docker images -qf dangling=true) 2>/dev/null || true
    echo "   ✅ Removed $DANGLING dangling images"
else
    echo "   ℹ️  No dangling images to remove"
fi
echo ""

# Step 3: Remove unused images (older than 24h)
echo "🗑️  Step 3: Removing unused images (older than 24h)..."
docker image prune -a --filter "until=24h" -f 2>/dev/null || true
echo "   ✅ Unused images pruned"
echo ""

# Step 4: Remove all stopped containers and unused networks
echo "🗑️  Step 4: Removing unused networks..."
docker network prune -f 2>/dev/null || true
echo "   ✅ Unused networks pruned"
echo ""

# Step 5: Remove build cache
echo "🗑️  Step 5: Removing Docker build cache..."
docker builder prune -f 2>/dev/null || true
echo "   ✅ Build cache cleared"
echo ""

# Step 6: Remove volumes (CAREFUL - only if you're sure!)
echo "⚠️  Step 6: Docker volumes (SKIPPED by default)"
echo "   To remove unused volumes, run: docker volume prune -f"
echo "   WARNING: This deletes data! Only do if you're sure."
echo ""

# Step 7: System-wide Docker prune
echo "🗑️  Step 7: System-wide Docker cleanup..."
docker system prune -f 2>/dev/null || true
echo "   ✅ System pruned"
echo ""

# Step 8: Clean apt cache on host (if accessible)
if [ -d "/var/cache/apt/archives" ] && [ -w "/var/cache/apt/archives" ]; then
    echo "🗑️  Step 8: Cleaning host apt cache..."
    rm -rf /var/cache/apt/archives/* 2>/dev/null || echo "   ⚠️  Cannot write to apt cache (no sudo)"
    echo "   ✅ Host apt cache cleaned"
else
    echo "🗑️  Step 8: Skipping apt cache (no write access - requires sudo)"
fi
echo ""

# Step 9: Clean tmp directories
echo "🗑️  Step 9: Cleaning temporary directories..."
rm -rf /tmp/docker-* 2>/dev/null || true
rm -rf /tmp/build-* 2>/dev/null || true
rm -rf /tmp/*.tar.gz 2>/dev/null || true
echo "   ✅ Temp directories cleaned"
echo ""

# Show final disk space
echo "📊 Final disk usage:"
df -h /var/lib/docker 2>/dev/null || df -h /
echo ""

# Show Docker disk usage
echo "🐳 Docker disk usage:"
docker system df
echo ""

echo "✅ Cleanup complete!"
echo ""
echo "💡 Tips to maintain space:"
echo "   1. Run this script before big builds"
echo "   2. Use multi-stage builds to reduce layer count"
echo "   3. Clean up in same RUN layer: && rm -rf /tmp/*"
echo "   4. Use .dockerignore to avoid large build contexts"
echo "   5. Regularly prune: docker system prune -a --volumes"
echo ""
