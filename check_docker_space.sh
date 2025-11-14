#!/bin/bash
# Quick script to check and optionally clean Docker space

echo "🔍 Checking Docker disk usage..."
echo "=================================="
echo ""

echo "📊 Docker system df:"
docker system df
echo ""

echo "📦 Docker images:"
docker images --format "table {{.Repository}}\t{{.Tag}}\t{{.Size}}"
echo ""

echo "💾 Disk space in /var/lib/docker (if accessible):"
if [ -d "/var/lib/docker" ]; then
    du -sh /var/lib/docker 2>/dev/null || echo "Cannot access /var/lib/docker (may need sudo)"
else
    echo "/var/lib/docker not found (Docker may use different location)"
fi
echo ""

echo "💾 Disk space on current filesystem:"
df -h . | tail -1
echo ""

echo "🧹 Cleanup options:"
echo "   1. docker system prune -a    # Remove all unused images"
echo "   2. docker system prune --volumes  # Remove unused volumes too"
echo "   3. docker builder prune -a  # Remove build cache"
echo ""

read -p "Would you like to run cleanup? (1/2/3/N): " -n 1 -r
echo ""

case $REPLY in
    1)
        echo "Running: docker system prune -a"
        docker system prune -a
        ;;
    2)
        echo "Running: docker system prune --volumes"
        docker system prune --volumes
        ;;
    3)
        echo "Running: docker builder prune -a"
        docker builder prune -a
        ;;
    *)
        echo "No cleanup performed"
        ;;
esac

echo ""
echo "✅ Done"
