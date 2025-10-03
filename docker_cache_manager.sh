#!/bin/bash
#
# Docker Build Cache Manager
# Helps manage Docker build cache for faster rebuilds after failures
#

set -euo pipefail

echo "🗄️  Docker Build Cache Manager"
echo "============================="

# Functions
show_help() {
    cat << EOF
Usage: $0 COMMAND [OPTIONS]

Commands:
  status      Show cache status and available images
  clean       Clean old cache images
  prune       Prune all Docker cache
  resume      Resume build from last successful layer
  save        Save current build as cache image
  restore     Restore from cache image

Options:
  --all       Apply to all cache types
  --force     Force operation without confirmation

Examples:
  $0 status                    # Show cache information
  $0 clean                     # Clean old cache images
  $0 resume                    # Resume failed build
  $0 save micapipe:backup      # Save current state as backup
EOF
}

show_cache_status() {
    echo "📊 Docker Cache Status"
    echo "====================="
    
    echo ""
    echo "🐳 MICApipe Images:"
    docker images | grep -E "(micapipe|none)" | head -10 || echo "   No MICApipe images found"
    
    echo ""
    echo "💾 Cache Usage:"
    docker system df
    
    echo ""
    echo "📦 Build Cache:"
    docker builder du 2>/dev/null || echo "   BuildKit cache info not available"
    
    echo ""
    echo "🏷️  Recent Build Logs:"
    if [[ -d "build_logs" ]]; then
        ls -la build_logs/*.log 2>/dev/null | tail -5 || echo "   No build logs found"
    else
        echo "   No build_logs directory found"
    fi
}

clean_cache() {
    local force=${1:-false}
    
    echo "🧹 Cleaning Docker Cache"
    echo "========================"
    
    # Show what will be cleaned
    echo "Will remove:"
    echo "  - Dangling images"
    echo "  - Unused containers"  
    echo "  - Old MICApipe cache images (keep latest 3)"
    echo ""
    
    if [[ "$force" != "true" ]]; then
        read -p "Continue? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Cancelled"
            return 0
        fi
    fi
    
    echo "🗑️  Removing dangling images..."
    docker image prune -f
    
    echo "🗑️  Removing unused containers..."
    docker container prune -f
    
    echo "🗑️  Cleaning old MICApipe cache images..."
    # Keep latest 3 cache images, remove older ones
    docker images | grep "micapipe.*cache" | tail -n +4 | awk '{print $3}' | xargs -r docker rmi || true
    
    echo "✅ Cache cleanup complete"
}

resume_build() {
    echo "🔄 Resume Build from Cache"
    echo "========================="
    
    # Find the most recent failed build
    if [[ -d "build_logs" ]]; then
        LATEST_LOG=$(ls -t build_logs/*.log 2>/dev/null | head -1 || echo "")
        if [[ -n "$LATEST_LOG" ]]; then
            echo "📋 Last build log: $LATEST_LOG"
            
            # Check where the build failed
            LAST_SUCCESS=$(grep -n "Successfully built\|---> [a-f0-9]" "$LATEST_LOG" | tail -1 | cut -d: -f1 || echo "0")
            TOTAL_STEPS=$(grep -c "Step [0-9]" "$LATEST_LOG" || echo "unknown")
            
            echo "📊 Build progress: $LAST_SUCCESS / $TOTAL_STEPS steps completed"
            echo ""
            echo "🚀 Resuming build with cache..."
            
            # Resume with cache enabled
            ./build_container_server.sh --cache-tag micapipe:resume
        else
            echo "❌ No build logs found"
            echo "   Starting fresh build..."
            ./build_container_server.sh
        fi
    else
        echo "❌ No build_logs directory found"
        echo "   Starting fresh build..."
        ./build_container_server.sh
    fi
}

save_cache() {
    local cache_tag=${1:-"micapipe:cache-$(date +%Y%m%d_%H%M%S)"}
    
    echo "💾 Saving Build Cache"
    echo "===================="
    
    if docker image inspect micapipe:latest >/dev/null 2>&1; then
        echo "📦 Tagging micapipe:latest as $cache_tag"
        docker tag micapipe:latest "$cache_tag"
        echo "✅ Cache saved as: $cache_tag"
    else
        echo "❌ No micapipe:latest image found to save"
        return 1
    fi
}

# Parse command line
if [[ $# -eq 0 ]]; then
    show_help
    exit 0
fi

COMMAND=$1
shift

FORCE=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --force)
            FORCE=true
            shift
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

# Execute command
case $COMMAND in
    status)
        show_cache_status
        ;;
    clean)
        clean_cache $FORCE
        ;;
    prune)
        echo "🧹 Pruning all Docker cache..."
        if [[ "$FORCE" == "true" ]] || { echo "This will remove ALL Docker cache. Continue? (y/N)"; read -n 1 -r; [[ $REPLY =~ ^[Yy]$ ]]; }; then
            docker system prune -af --volumes
            echo "✅ All Docker cache pruned"
        fi
        ;;
    resume)
        resume_build
        ;;
    save)
        save_cache "$@"
        ;;
    restore)
        if [[ $# -eq 0 ]]; then
            echo "❌ Please specify cache image to restore from"
            echo "Usage: $0 restore <cache-image-tag>"
            exit 1
        fi
        CACHE_IMAGE=$1
        echo "📦 Restoring from cache: $CACHE_IMAGE"
        docker tag "$CACHE_IMAGE" micapipe:latest
        echo "✅ Restored micapipe:latest from $CACHE_IMAGE"
        ;;
    *)
        echo "❌ Unknown command: $COMMAND"
        show_help
        exit 1
        ;;
esac