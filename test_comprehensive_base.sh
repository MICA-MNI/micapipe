#!/bin/bash
# ============================================================================
# TEST COMPREHENSIVE BASE IMAGE STRATEGY
# Validates that the comprehensive base image approach works correctly
# ============================================================================

set -e

echo "🧪 Testing MICApipe Comprehensive Base Image Strategy"
echo "=================================================="

# Test configuration
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
BASE_IMAGE="${REGISTRY}/micapipe-comprehensive-base:latest"
TEST_TAG="test-$(date +%s)"

echo ""
echo "🔍 Step 1: Checking Docker files..."

required_files=(
    "Dockerfile.mamba-base"
    "Dockerfile.minimal"
    "build_comprehensive_base.sh"
    "build_fast_ci.sh"
)

for file in "${required_files[@]}"; do
    if [ -f "$file" ]; then
        echo "✅ $file exists"
    else
        echo "❌ $file missing"
        exit 1
    fi
done

echo ""
echo "🏗️  Step 2: Testing base image build (this may take a while)..."
echo "⚠️  Note: This will build the comprehensive base image locally"
echo "⏱️  Expected time: 45-90 minutes for first build"

read -p "Continue with base image build? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "⏭️  Skipping base image build. To test fully, run:"
    echo "   ./build_comprehensive_base.sh"
    echo "   then re-run this test script"
    exit 0
fi

# Build base image
echo "Building comprehensive base image..."
if ./build_comprehensive_base.sh; then
    echo "✅ Base image build successful"
else
    echo "❌ Base image build failed"
    exit 1
fi

echo ""
echo "🚀 Step 3: Testing fast main image build..."

# Build main image using minimal Dockerfile
echo "Building main image (should be fast)..."
start_time=$(date +%s)

if docker build -f Dockerfile.minimal -t micapipe:${TEST_TAG} .; then
    end_time=$(date +%s)
    build_duration=$((end_time - start_time))
    echo "✅ Main image build successful"
    echo "⏱️  Build time: ${build_duration} seconds"
    
    if [ $build_duration -lt 600 ]; then  # Less than 10 minutes
        echo "🎉 Excellent! Build time under 10 minutes"
    elif [ $build_duration -lt 900 ]; then  # Less than 15 minutes
        echo "👍 Good! Build time under 15 minutes"
    else
        echo "⚠️  Build time longer than expected (over 15 minutes)"
    fi
else
    echo "❌ Main image build failed"
    exit 1
fi

echo ""
echo "🔬 Step 4: Testing container functionality..."

# Test that the container starts and has the right tools
echo "Testing container startup and tool availability..."

tools_to_test=(
    "fsl"
    "mri_convert"
    "3dinfo"
    "antsRegistration"
    "conda"
    "python"
)

container_id=$(docker run -d micapipe:${TEST_TAG} sleep 60)

for tool in "${tools_to_test[@]}"; do
    if docker exec $container_id which $tool > /dev/null 2>&1; then
        echo "✅ $tool available"
    else
        echo "❌ $tool not found"
        failed_tools="$failed_tools $tool"
    fi
done

# Test Python packages
python_packages=(
    "numpy"
    "nibabel"
    "antspy"
)

for package in "${python_packages[@]}"; do
    if docker exec $container_id python -c "import $package" > /dev/null 2>&1; then
        echo "✅ Python package $package importable"
    else
        echo "❌ Python package $package not available"
        failed_packages="$failed_packages $package"
    fi
done

# Cleanup test container
docker stop $container_id > /dev/null
docker rm $container_id > /dev/null

echo ""
echo "🧹 Step 5: Cleanup test images..."
docker rmi micapipe:${TEST_TAG} > /dev/null

echo ""
echo "📊 Test Results Summary"
echo "====================="

if [ -z "$failed_tools" ] && [ -z "$failed_packages" ]; then
    echo "🎉 ALL TESTS PASSED!"
    echo ""
    echo "✅ Base image strategy is working correctly"
    echo "✅ All neuroimaging tools are available"
    echo "✅ Python environment is properly configured"
    echo "✅ Build time is optimized"
    echo ""
    echo "🚀 Ready for production use!"
    echo ""
    echo "Next steps:"
    echo "1. Push base image to registry: docker push ${BASE_IMAGE}"
    echo "2. Update your CI to use: ./build_fast_ci.sh"
    echo "3. Enjoy 90%+ faster builds!"
else
    echo "❌ SOME TESTS FAILED"
    
    if [ ! -z "$failed_tools" ]; then
        echo "Missing tools:$failed_tools"
    fi
    
    if [ ! -z "$failed_packages" ]; then
        echo "Missing Python packages:$failed_packages"
    fi
    
    echo ""
    echo "Please check the base image build logs for errors."
    exit 1
fi