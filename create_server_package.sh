#!/bin/bash

# Package MICApipe for Server Deployment
# Creates a deployment package with all necessary files and scripts

echo "📦 MICApipe Server Deployment Package Creator"
echo "============================================"

PACKAGE_NAME="micapipe-server-deployment-$(date +%Y%m%d_%H%M%S)"
PACKAGE_DIR="./$PACKAGE_NAME"

echo "Creating deployment package: $PACKAGE_NAME"

# Create package directory
mkdir -p "$PACKAGE_DIR"

# Core files
echo "📋 Copying core files..."
cp Dockerfile "$PACKAGE_DIR/"
cp micapipe.py "$PACKAGE_DIR/"
cp pyproject.toml "$PACKAGE_DIR/"
cp config.yaml "$PACKAGE_DIR/"
cp README.md "$PACKAGE_DIR/"
cp LICENSE "$PACKAGE_DIR/"

# Copy essential directories
echo "📁 Copying directories..."
cp -r functions "$PACKAGE_DIR/"
cp -r parcellations "$PACKAGE_DIR/"
cp -r surfaces "$PACKAGE_DIR/"
cp -r MNI152Volumes "$PACKAGE_DIR/"
cp -r MICs60_T1-atlas "$PACKAGE_DIR/"
cp -r fsl_conf "$PACKAGE_DIR/"

# Build scripts
echo "🔧 Copying build scripts..."
cp build_no_sudo.sh "$PACKAGE_DIR/"
cp server_build_test.sh "$PACKAGE_DIR/"
cp fix_fsl_build.sh "$PACKAGE_DIR/"
cp test_fsl_fix.sh "$PACKAGE_DIR/"

# Documentation
echo "📚 Copying documentation..."
cp SERVER_DEPLOYMENT.md "$PACKAGE_DIR/"
cp BUILD_GUIDE.md "$PACKAGE_DIR/" 2>/dev/null || echo "BUILD_GUIDE.md not found, skipping"

# Create a deployment script inside the package
cat > "$PACKAGE_DIR/deploy.sh" << 'EOF'
#!/bin/bash

# MICApipe Server Deployment Script
# Run this script on your server to deploy MICApipe

echo "🚀 MICApipe Server Deployment"
echo "============================"

# Check if Docker is available
if ! command -v docker &> /dev/null; then
    echo "❌ Docker not found. Please install Docker first."
    echo "   Visit: https://docs.docker.com/get-docker/"
    exit 1
fi

echo "✅ Docker found: $(docker --version)"

# Make scripts executable
chmod +x *.sh

echo ""
echo "📋 Available deployment options:"
echo "1. Complete build and test (recommended): ./server_build_test.sh"
echo "2. Quick build without sudo: ./build_no_sudo.sh"
echo "3. FSL-specific testing: ./server_build_test.sh fsl-only"
echo "4. Build without testing: ./server_build_test.sh no-test"
echo ""
echo "📖 For detailed instructions, see: SERVER_DEPLOYMENT.md"
echo ""

# Ask user for deployment choice
read -p "Choose deployment option (1-4) or press Enter for option 1: " choice
choice=${choice:-1}

case $choice in
    1)
        echo "🚀 Running complete build and test..."
        ./server_build_test.sh
        ;;
    2)
        echo "🚀 Running quick build..."
        ./build_no_sudo.sh
        ;;
    3)
        echo "🧪 Running FSL-specific test..."
        ./server_build_test.sh fsl-only
        ;;
    4)
        echo "🏗️ Running build without testing..."
        ./server_build_test.sh no-test
        ;;
    *)
        echo "❌ Invalid choice. Please run one of the scripts manually."
        exit 1
        ;;
esac

echo ""
echo "🎉 Deployment completed!"
echo "📊 Check build_logs/ directory for detailed logs and reports."
EOF

chmod +x "$PACKAGE_DIR/deploy.sh"

# Create README for the package
cat > "$PACKAGE_DIR/README_DEPLOYMENT.txt" << EOF
MICApipe Server Deployment Package
==================================

This package contains everything needed to deploy MICApipe on your server.

QUICK START:
1. Upload this entire directory to your server
2. Run: ./deploy.sh
3. Follow the prompts

MANUAL DEPLOYMENT:
1. Make scripts executable: chmod +x *.sh
2. Run deployment script: ./server_build_test.sh

FILES INCLUDED:
- Dockerfile: Updated container definition
- build scripts: Various build options for different server configurations
- SERVER_DEPLOYMENT.md: Comprehensive deployment guide
- deploy.sh: Interactive deployment script

REQUIREMENTS:
- Docker installed and running
- At least 10GB free disk space
- At least 4GB RAM
- Internet connection for downloading dependencies

SUPPORT:
- Check SERVER_DEPLOYMENT.md for troubleshooting
- Review build logs in build_logs/ directory
- Use fix_fsl_build.sh if FSL installation issues occur

Updated tool versions included:
- MRtrix 3.0.7
- FreeSurfer 7.4.1  
- FastSurfer v2.4.2
- DESIGNER pipeline
- Synb0-DISCO & SynBOLD-DisCo
- LAMAReg
- SWM analysis tools
EOF

# Create archive
echo "📦 Creating deployment archive..."
tar -czf "${PACKAGE_NAME}.tar.gz" "$PACKAGE_DIR"

# Summary
echo ""
echo "✅ Deployment package created successfully!"
echo ""
echo "📦 Package details:"
echo "   Name: $PACKAGE_NAME"
echo "   Archive: ${PACKAGE_NAME}.tar.gz"
echo "   Size: $(du -sh "${PACKAGE_NAME}.tar.gz" | cut -f1)"
echo ""
echo "🚀 To deploy on server:"
echo "   1. Upload ${PACKAGE_NAME}.tar.gz to your server"
echo "   2. Extract: tar -xzf ${PACKAGE_NAME}.tar.gz"
echo "   3. Enter directory: cd $PACKAGE_NAME"
echo "   4. Run deployment: ./deploy.sh"
echo ""
echo "📚 For detailed instructions, see SERVER_DEPLOYMENT.md in the package."

# Cleanup
rm -rf "$PACKAGE_DIR"

echo "📁 Package ready for server deployment!"