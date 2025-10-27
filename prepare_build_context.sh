#!/bin/bash
#
# Prepare build context with pre-downloaded files for Docker build
# This script copies pre-downloaded neuroimaging tools to the build directory
# to speed up Docker builds by avoiding downloads during build time.
#

set -e

echo "📦 Preparing build context with pre-downloaded files"
echo "=================================================="

# Source and destination directories
SOURCE_DIR="/export02/data/enning"
BUILD_DIR="$(pwd)"

echo "📁 Copying pre-downloaded files from ${SOURCE_DIR} to build context..."

# Copy pre-downloaded files to build context
# These will be available as build context for Docker COPY commands

# FreeSurfer
if [ -f "${SOURCE_DIR}/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" ]; then
    echo "✅ Copying freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz..."
    cp "${SOURCE_DIR}/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" "${BUILD_DIR}/"
fi

# FSL
if [ -f "${SOURCE_DIR}/fsl-6.0.2-centos6_64.tar.gz" ]; then
    echo "✅ Copying fsl-6.0.2-centos6_64.tar.gz..."
    cp "${SOURCE_DIR}/fsl-6.0.2-centos6_64.tar.gz" "${BUILD_DIR}/"
fi

# AFNI
if [ -f "${SOURCE_DIR}/afni-linux_openmp_64.tgz" ]; then
    echo "✅ Copying afni-linux_openmp_64.tgz..."
    cp "${SOURCE_DIR}/afni-linux_openmp_64.tgz" "${BUILD_DIR}/"
fi

# FSL FIX
if [ -f "${SOURCE_DIR}/fix-1.068.tar.gz" ]; then
    echo "✅ Copying fix-1.068.tar.gz..."
    cp "${SOURCE_DIR}/fix-1.068.tar.gz" "${BUILD_DIR}/"
fi

# Miniconda
if [ -f "${SOURCE_DIR}/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh" ]; then
    echo "✅ Copying Miniconda3-py39_22.11.1-1-Linux-x86_64.sh..."
    cp "${SOURCE_DIR}/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh" "${BUILD_DIR}/"
fi

echo "✅ Build context preparation complete!"
echo "📦 Files copied to: ${BUILD_DIR}"