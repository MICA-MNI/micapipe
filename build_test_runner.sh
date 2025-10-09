#!/bin/bash
#
# Build a test version of the action runner
# This creates micapipe-runner-test image for testing base + main
#

set -e

ACTIONS_RUNNER_DIR="/data_/mica1/03_projects/actions-runner"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "🔧 BUILDING TEST ACTION RUNNER"
log "=============================="

# Check if we're in the right directory
if [ ! -f "$ACTIONS_RUNNER_DIR/Dockerfile" ]; then
    log "❌ Action runner Dockerfile not found at: $ACTIONS_RUNNER_DIR"
    log "   Make sure you're running this from the correct location"
    exit 1
fi

log "📋 Building micapipe-runner-test image..."

cd "$ACTIONS_RUNNER_DIR"

# Build with test tag
docker build -t micapipe-runner-test .

log "✅ Test runner image built: micapipe-runner-test"

# Create test run script
cat > run_docker_test.sh << 'EOF'
#!/bin/bash
export DOCKER_CONTENT_TRUST=0

echo "🧪 Starting micapipe-runner-test..."

# Stop and remove if exists
docker stop micapipe-runner-test 2>/dev/null || true
docker rm micapipe-runner-test 2>/dev/null || true

# Run test runner
docker run -d --name micapipe-runner-test --restart always --privileged  \
    -v /data/mica1/03_projects/enning/BIDS_CI/sing_cache:/data/mica1/03_projects/enning/BIDS_CI/sing_cache \
    -v /export02/local/singularity_tmp/:/export02/local/singularity_tmp/ \
    -v /data/mica3/BIDS_CI/rawdata:/data/mica3/BIDS_CI/rawdata \
    -v /data_/mica3/BIDS_CI:/data_/mica3/BIDS_CI \
    -v /data/mica1/03_projects/enning/BIDS_CI:/data/mica1/03_projects/enning/BIDS_CI \
    -v /tmp:/tmp \
    -v /var/run/docker.sock:/var/run/docker.sock \
    micapipe-runner-test

echo "✅ micapipe-runner-test started"
echo "🔒 Your original micapipe-runner remains untouched"

# Check status
docker ps | grep micapipe-runner
EOF

chmod +x run_docker_test.sh

log "✅ Test run script created: $ACTIONS_RUNNER_DIR/run_docker_test.sh"

log "=============================="
log "✅ TEST RUNNER BUILD COMPLETE"
log "=============================="
log "📦 Image: micapipe-runner-test"
log "🚀 Start with:"
log "   cd $ACTIONS_RUNNER_DIR"
log "   ./run_docker_test.sh"
log ""
log "🔒 Your original runner is safe and unchanged"
log "🧪 Test the base + main structure with the test runner"