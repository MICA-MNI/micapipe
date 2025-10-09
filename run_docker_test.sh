#!/bin/bash
#
# Create a test action runner for base + main Docker structure
# This runs alongside your existing working runner
#

export DOCKER_CONTENT_TRUST=0

echo "🔧 Creating micapipe-runner-test for base + main structure..."

# Stop and remove if exists
docker stop micapipe-runner-test 2>/dev/null || true
docker rm micapipe-runner-test 2>/dev/null || true

# Run test runner on different port/name
docker run -d --name micapipe-runner-test --restart always --privileged  \
    -v /data/mica1/03_projects/enning/BIDS_CI/sing_cache:/data/mica1/03_projects/enning/BIDS_CI/sing_cache \
    -v /export02/local/singularity_tmp/:/export02/local/singularity_tmp/ \
    -v /data/mica3/BIDS_CI/rawdata:/data/mica3/BIDS_CI/rawdata \
    -v /data_/mica3/BIDS_CI:/data_/mica3/BIDS_CI \
    -v /data/mica1/03_projects/enning/BIDS_CI:/data/mica1/03_projects/enning/BIDS_CI \
    -v /tmp:/tmp \
    -v /var/run/docker.sock:/var/run/docker.sock \
    micapipe-runner

echo "✅ micapipe-runner-test started"
echo "📋 This runner will test the base + main structure"
echo "🔒 Your original micapipe-runner remains untouched"

# Check status
docker ps | grep micapipe-runner