#!/bin/bash
#
# Emergency Docker Recovery Script
# Force kills stuck Docker processes and cleans up system
#

set -e

echo "🚨 Emergency Docker Recovery"
echo "============================"

echo "🔍 Finding stuck Docker processes..."

# Find and kill Docker build processes
DOCKER_PIDS=$(ps aux | grep -E "(docker build|docker-proxy|dockerd)" | grep -v grep | awk '{print $2}' || true)

if [[ -n "$DOCKER_PIDS" ]]; then
    echo "📋 Found Docker processes: $DOCKER_PIDS"
    echo "💀 Force killing stuck Docker processes..."
    
    for PID in $DOCKER_PIDS; do
        echo "   Killing PID: $PID"
        sudo kill -9 $PID 2>/dev/null || true
    done
else
    echo "✅ No stuck Docker processes found"
fi

echo ""
echo "🐋 Restarting Docker daemon..."

# Method 1: Try systemctl (most Linux systems)
if command -v systemctl >/dev/null 2>&1; then
    echo "   Using systemctl..."
    sudo systemctl stop docker.socket docker.service || true
    sleep 5
    sudo systemctl start docker.service || true
    sleep 5
fi

# Method 2: Try service command (older systems)
if command -v service >/dev/null 2>&1; then
    echo "   Using service command..."
    sudo service docker stop || true
    sleep 5
    sudo service docker start || true
    sleep 5
fi

# Method 3: Direct daemon restart (if above fail)
if ! docker info >/dev/null 2>&1; then
    echo "   Trying direct daemon restart..."
    sudo pkill -f dockerd || true
    sleep 5
    sudo dockerd >/dev/null 2>&1 &
    sleep 10
fi

echo ""
echo "🧹 Cleaning up Docker resources..."

# Wait for Docker to be responsive
for i in {1..30}; do
    if docker info >/dev/null 2>&1; then
        echo "✅ Docker daemon is responsive"
        break
    else
        echo "   Waiting for Docker daemon... ($i/30)"
        sleep 2
    fi
done

# Clean up everything
if docker info >/dev/null 2>&1; then
    echo "   Stopping all containers..."
    docker ps -q | xargs -r docker kill || true
    
    echo "   Removing all containers..."
    docker ps -aq | xargs -r docker rm -f || true
    
    echo "   Pruning system..."
    docker system prune -af --volumes || true
    
    echo "   Removing unused images..."
    docker image prune -af || true
else
    echo "❌ Docker daemon still not responsive"
    echo "   You may need to reboot the system"
fi

echo ""
echo "📊 System Status:"
echo "Memory usage: $(free -h | grep Mem | awk '{print $3 "/" $2}')"
echo "Disk usage: $(df -h . | tail -1 | awk '{print $3 "/" $2 " (" $5 ")"}')"

if docker info >/dev/null 2>&1; then
    echo "Docker status: ✅ Running"
    echo "Docker containers: $(docker ps -q | wc -l) running"
    echo "Docker images: $(docker images -q | wc -l) total"
else
    echo "Docker status: ❌ Not responding"
fi

echo ""
echo "💡 Next steps:"
echo "1. If Docker is running, try building again with debug_build.sh"
echo "2. If Docker is still stuck, you may need to reboot the system"
echo "3. Consider increasing system memory or freeing disk space"