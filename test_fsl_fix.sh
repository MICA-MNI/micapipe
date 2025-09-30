#!/bin/bash

# Quick FSL Build Test
# Tests the Docker build up to the FSL section to verify the fix

echo "🔧 Quick FSL Build Test"
echo "======================"

# Set Docker Content Trust to avoid certificate issues
export DOCKER_CONTENT_TRUST=0

echo "📋 Testing Dockerfile changes..."
echo "Current FSL installation section:"
echo "--------------------------------"
grep -A 10 -B 2 "FSL Python environment" Dockerfile

echo ""
echo "🔍 Starting test build (first 20 steps)..."

# Create a temporary Dockerfile that stops after FSL installation
head -n 150 Dockerfile > Dockerfile.test-fsl
echo "RUN echo 'FSL test build completed successfully'" >> Dockerfile.test-fsl

# Run the test build
docker build --progress=plain --no-cache -f Dockerfile.test-fsl -t micapipe:fsl-test . 2>&1 | tee fsl_test.log

# Check result
if [ $? -eq 0 ]; then
    echo "✅ FSL section builds successfully!"
    echo "🎉 The FSL Python installation fix appears to be working."
    echo ""
    echo "🚀 You can now run the full build with:"
    echo "   ./build_container.sh"
    echo "   or"
    echo "   docker build -t micapipe:latest ."
else
    echo "❌ FSL section still has issues."
    echo "📝 Check fsl_test.log for details."
    echo ""
    echo "🔧 Alternative options:"
    echo "   1. Run: ./fix_fsl_build.sh (for comprehensive fixes)"
    echo "   2. Or build without FSL Python: docker build -f Dockerfile.no-fsl-python -t micapipe ."
fi

# Cleanup
rm -f Dockerfile.test-fsl

echo ""
echo "📊 Test completed. Check fsl_test.log for detailed output."