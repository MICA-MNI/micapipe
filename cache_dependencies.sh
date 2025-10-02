#!/bin/bash
#
# Cache Management Script for Large Dependencies
# Downloads and caches FSL and FreeSurfer locally to speed up future builds
#

set -e

# Configuration
CACHE_DIR="/data_/mica1/01_programs/micapipe_cache"
FSL_URL="https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.2-centos6_64.tar.gz"
FREESURFER_URL="ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"

echo "🗂️  MICApipe Dependency Cache Manager"
echo "=================================="

# Create cache directory
mkdir -p "$CACHE_DIR"
echo "📁 Cache directory: $CACHE_DIR"

# Function to download and validate large files
download_and_cache() {
    local url="$1"
    local filename="$2"
    local expected_min_size="$3"
    local filepath="$CACHE_DIR/$filename"
    
    if [ -f "$filepath" ] && [ $(stat -c%s "$filepath") -gt $expected_min_size ]; then
        echo "✅ $filename already cached ($(du -h "$filepath" | cut -f1))"
        return 0
    fi
    
    echo "⬇️  Downloading $filename..."
    echo "   URL: $url"
    
    # Try multiple download methods
    if curl -fsSL --retry 3 --retry-delay 10 --connect-timeout 30 --max-time 3600 \
            "$url" -o "$filepath.tmp"; then
        echo "✅ Download completed"
    elif wget --timeout=300 --tries=3 --retry-connrefused --waitretry=30 \
            "$url" -O "$filepath.tmp"; then
        echo "✅ Download completed with wget"
    else
        echo "❌ Download failed for $filename"
        rm -f "$filepath.tmp"
        return 1
    fi
    
    # Validate file size
    if [ ! -f "$filepath.tmp" ] || [ $(stat -c%s "$filepath.tmp") -lt $expected_min_size ]; then
        echo "❌ Download incomplete or corrupted"
        rm -f "$filepath.tmp"
        return 1
    fi
    
    # Move to final location
    mv "$filepath.tmp" "$filepath"
    echo "✅ $filename cached successfully ($(du -h "$filepath" | cut -f1))"
    return 0
}

# Download FSL (expected ~600MB compressed)
echo ""
echo "🧠 Caching FSL 6.0.2..."
download_and_cache "$FSL_URL" "fsl-6.0.2-centos6_64.tar.gz" 500000000

# Download FreeSurfer (expected ~3GB compressed)
echo ""
echo "🧠 Caching FreeSurfer 7.4.1..."
download_and_cache "$FREESURFER_URL" "freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz" 2000000000

echo ""
echo "📊 Cache Summary:"
ls -lh "$CACHE_DIR"

echo ""
echo "🎯 Cache ready! Build times will be much faster now."
echo "💡 To use cached files, set environment variable:"
echo "   export MICAPIPE_CACHE_DIR=$CACHE_DIR"