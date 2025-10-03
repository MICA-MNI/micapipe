#!/bin/bash
#
# Download Dependencies Script for MICApipe
# Downloads FSL and FreeSurfer to a specified directory before Docker build
# This allows pre-downloading large files to avoid timeouts during builds
#

set -e

# Default configuration
DOWNLOAD_DIR="/host/cassio/export03/data/enning/downloads"
FSL_URL="https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.2-centos6_64.tar.gz"
FREESURFER_URL="ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
AFNI_URL="https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz"
FSL_FIX_URL="https://git.fmrib.ox.ac.uk/fsl/fix/-/archive/1.068/fix-1.068.tar.gz"
FSL_FILENAME="fsl-6.0.2-centos6_64.tar.gz"
FREESURFER_FILENAME="freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
AFNI_FILENAME="afni-linux_openmp_64.tgz"
FSL_FIX_FILENAME="fix-1.068.tar.gz"
FSL_MIN_SIZE=$((500 * 1024 * 1024))      # 500MB
FREESURFER_MIN_SIZE=$((2 * 1024 * 1024 * 1024))  # 2GB
AFNI_MIN_SIZE=$((300 * 1024 * 1024))     # 300MB
FSL_FIX_MIN_SIZE=$((30 * 1024 * 1024))   # 30MB

# Options
DOWNLOAD_FSL=true
DOWNLOAD_FREESURFER=true
DOWNLOAD_AFNI=true
DOWNLOAD_FSL_FIX=true
FORCE_DOWNLOAD=false
VERIFY_ONLY=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --download-dir)
            DOWNLOAD_DIR="$2"
            shift 2
            ;;
        --fsl-only)
            DOWNLOAD_FSL=true
            DOWNLOAD_FREESURFER=false
            DOWNLOAD_AFNI=false
            DOWNLOAD_FSL_FIX=false
            shift
            ;;
        --freesurfer-only)
            DOWNLOAD_FSL=false
            DOWNLOAD_FREESURFER=true
            DOWNLOAD_AFNI=false
            DOWNLOAD_FSL_FIX=false
            shift
            ;;
        --afni-only)
            DOWNLOAD_FSL=false
            DOWNLOAD_FREESURFER=false
            DOWNLOAD_AFNI=true
            DOWNLOAD_FSL_FIX=false
            shift
            ;;
        --fix-only)
            DOWNLOAD_FSL=false
            DOWNLOAD_FREESURFER=false
            DOWNLOAD_AFNI=false
            DOWNLOAD_FSL_FIX=true
            shift
            ;;
        --no-afni)
            DOWNLOAD_AFNI=false
            shift
            ;;
        --no-fix)
            DOWNLOAD_FSL_FIX=false
            shift
            ;;
        --force)
            FORCE_DOWNLOAD=true
            shift
            ;;
        --verify)
            VERIFY_ONLY=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Downloads large dependencies for MICApipe container builds"
            echo ""
            echo "Options:"
            echo "  --download-dir DIR    Directory to download files (default: /host/cassio/export03/data/enning/downloads)"
            echo "  --fsl-only           Download only FSL (~600MB)"
            echo "  --freesurfer-only    Download only FreeSurfer (~3GB)"
            echo "  --afni-only          Download only AFNI (~500MB)"
            echo "  --fix-only           Download only FSL FIX (~50MB)"
            echo "  --no-afni            Skip AFNI download"
            echo "  --no-fix             Skip FSL FIX download"
            echo "  --force              Force re-download even if files exist"
            echo "  --verify             Only verify existing files, don't download"
            echo "  --help, -h           Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                                    # Download all dependencies to server location"
            echo "  $0 --download-dir ./downloads         # Use local directory instead"
            echo "  $0 --fsl-only                        # Download only FSL"
            echo "  $0 --no-afni --no-fix                # Download only FSL and FreeSurfer"
            echo "  $0 --verify                          # Check existing downloads"
            echo "  $0 --force                           # Re-download all files"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "🔽 MICApipe Dependencies Downloader"
echo "===================================="

# Check if default server directory is accessible, fallback to local if not
if [[ "$DOWNLOAD_DIR" == "/host/cassio/export03/data/enning/downloads" ]]; then
    if [[ ! -d "/host/cassio/export03/data/enning" ]] || [[ ! -w "/host/cassio/export03/data/enning" ]]; then
        echo "⚠️  Server directory not accessible, using local fallback"
        DOWNLOAD_DIR="./downloads"
    fi
fi

echo "📁 Download directory: $DOWNLOAD_DIR"
echo ""

# Create download directory
if [[ "$VERIFY_ONLY" == "false" ]]; then
    mkdir -p "$DOWNLOAD_DIR"
    if [[ ! -w "$DOWNLOAD_DIR" ]]; then
        echo "❌ Cannot write to download directory: $DOWNLOAD_DIR"
        exit 1
    fi
fi

# Function to format file size
format_size() {
    local size=$1
    if [[ $size -gt $((1024 * 1024 * 1024)) ]]; then
        echo "$(( size / 1024 / 1024 / 1024 ))GB"
    elif [[ $size -gt $((1024 * 1024)) ]]; then
        echo "$(( size / 1024 / 1024 ))MB"
    else
        echo "$(( size / 1024 ))KB"
    fi
}

# Function to verify file integrity
verify_file() {
    local filepath="$1"
    local min_size="$2"
    local filename=$(basename "$filepath")
    
    if [[ ! -f "$filepath" ]]; then
        echo "❌ $filename: File not found"
        return 1
    fi
    
    local actual_size=$(stat -c%s "$filepath" 2>/dev/null || stat -f%z "$filepath" 2>/dev/null || echo "0")
    local formatted_size=$(format_size $actual_size)
    
    if [[ $actual_size -lt $min_size ]]; then
        echo "❌ $filename: File too small ($formatted_size, expected ≥$(format_size $min_size))"
        return 1
    else
        echo "✅ $filename: Valid ($formatted_size)"
        return 0
    fi
}

# Function to download with multiple methods
download_file() {
    local url="$1"
    local filepath="$2"
    local min_size="$3"
    local filename=$(basename "$filepath")
    local temp_file="${filepath}.tmp"
    
    echo "🔄 Downloading $filename..."
    echo "   URL: $url"
    
    # Remove temporary file if it exists
    rm -f "$temp_file"
    
    # Try different download methods
    local download_success=false
    
    # Method 1: wget with progress bar
    if command -v wget >/dev/null 2>&1; then
        echo "   Using wget..."
        if wget --progress=bar --timeout=300 --tries=3 -O "$temp_file" "$url" 2>&1; then
            download_success=true
        else
            echo "   wget failed, trying curl..."
            rm -f "$temp_file"
        fi
    fi
    
    # Method 2: curl with progress bar
    if [[ "$download_success" == "false" ]] && command -v curl >/dev/null 2>&1; then
        echo "   Using curl..."
        if curl -L --progress-bar --connect-timeout 30 --max-time 3600 -o "$temp_file" "$url" 2>&1; then
            download_success=true
        else
            echo "   curl failed"
            rm -f "$temp_file"
        fi
    fi
    
    if [[ "$download_success" == "false" ]]; then
        echo "❌ Failed to download $filename"
        return 1
    fi
    
    # Verify downloaded file
    local actual_size=$(stat -c%s "$temp_file" 2>/dev/null || stat -f%z "$temp_file" 2>/dev/null || echo "0")
    if [[ $actual_size -lt $min_size ]]; then
        echo "❌ Downloaded file is too small ($(format_size $actual_size))"
        rm -f "$temp_file"
        return 1
    fi
    
    # Move to final location
    mv "$temp_file" "$filepath"
    echo "✅ Successfully downloaded $filename ($(format_size $actual_size))"
    return 0
}

# Function to process a dependency
process_dependency() {
    local name="$1"
    local url="$2"
    local filename="$3"
    local min_size="$4"
    local filepath="$DOWNLOAD_DIR/$filename"
    
    echo "📦 Processing $name..."
    
    # Check if file exists and is valid
    if verify_file "$filepath" "$min_size"; then
        if [[ "$FORCE_DOWNLOAD" == "true" && "$VERIFY_ONLY" == "false" ]]; then
            echo "🔄 Force download requested, re-downloading..."
            download_file "$url" "$filepath" "$min_size"
        else
            echo "✅ $name is already downloaded and valid"
        fi
    else
        if [[ "$VERIFY_ONLY" == "true" ]]; then
            echo "❌ $name verification failed"
            return 1
        else
            download_file "$url" "$filepath" "$min_size"
        fi
    fi
    echo ""
}

# Check network connectivity
if [[ "$VERIFY_ONLY" == "false" ]]; then
    echo "🌐 Checking network connectivity..."
    if curl -I --max-time 10 https://fsl.fmrib.ox.ac.uk >/dev/null 2>&1; then
        echo "✅ Network connectivity OK"
    else
        echo "⚠️  Network issues detected - downloads may fail"
    fi
    echo ""
fi

# Process FSL
FSL_SUCCESS=true
if [[ "$DOWNLOAD_FSL" == "true" ]]; then
    process_dependency "FSL 6.0.2" "$FSL_URL" "$FSL_FILENAME" "$FSL_MIN_SIZE" || FSL_SUCCESS=false
fi

# Process FreeSurfer
FREESURFER_SUCCESS=true
if [[ "$DOWNLOAD_FREESURFER" == "true" ]]; then
    process_dependency "FreeSurfer 7.4.1" "$FREESURFER_URL" "$FREESURFER_FILENAME" "$FREESURFER_MIN_SIZE" || FREESURFER_SUCCESS=false
fi

# Process AFNI
AFNI_SUCCESS=true
if [[ "$DOWNLOAD_AFNI" == "true" ]]; then
    process_dependency "AFNI Latest" "$AFNI_URL" "$AFNI_FILENAME" "$AFNI_MIN_SIZE" || AFNI_SUCCESS=false
fi

# Process FSL FIX
FSL_FIX_SUCCESS=true
if [[ "$DOWNLOAD_FSL_FIX" == "true" ]]; then
    process_dependency "FSL FIX 1.068" "$FSL_FIX_URL" "$FSL_FIX_FILENAME" "$FSL_FIX_MIN_SIZE" || FSL_FIX_SUCCESS=false
fi

# Summary
echo "📋 Download Summary"
echo "=================="
echo "📁 Directory: $DOWNLOAD_DIR"
if [[ -d "$DOWNLOAD_DIR" ]]; then
    echo "📊 Total size: $(du -sh "$DOWNLOAD_DIR" 2>/dev/null | cut -f1 || echo "0B")"
else
    echo "📊 Total size: Directory does not exist"
fi
echo ""
echo "Files:"
if [[ "$DOWNLOAD_FSL" == "true" ]]; then
    if [[ -f "$DOWNLOAD_DIR/$FSL_FILENAME" ]]; then
        verify_file "$DOWNLOAD_DIR/$FSL_FILENAME" "$FSL_MIN_SIZE" || true
    else
        echo "❌ $FSL_FILENAME: Not downloaded"
    fi
fi
if [[ "$DOWNLOAD_FREESURFER" == "true" ]]; then
    if [[ -f "$DOWNLOAD_DIR/$FREESURFER_FILENAME" ]]; then
        verify_file "$DOWNLOAD_DIR/$FREESURFER_FILENAME" "$FREESURFER_MIN_SIZE" || true
    else
        echo "❌ $FREESURFER_FILENAME: Not downloaded"
    fi
fi
if [[ "$DOWNLOAD_AFNI" == "true" ]]; then
    if [[ -f "$DOWNLOAD_DIR/$AFNI_FILENAME" ]]; then
        verify_file "$DOWNLOAD_DIR/$AFNI_FILENAME" "$AFNI_MIN_SIZE" || true
    else
        echo "❌ $AFNI_FILENAME: Not downloaded"
    fi
fi
if [[ "$DOWNLOAD_FSL_FIX" == "true" ]]; then
    if [[ -f "$DOWNLOAD_DIR/$FSL_FIX_FILENAME" ]]; then
        verify_file "$DOWNLOAD_DIR/$FSL_FIX_FILENAME" "$FSL_FIX_MIN_SIZE" || true
    else
        echo "❌ $FSL_FIX_FILENAME: Not downloaded"
    fi
fi

echo ""
echo "💡 Next steps:"
echo "   1. Use these files in your Dockerfile by copying them to the container"
echo "   2. Or mount the directory during docker build:"
echo "      docker build --build-arg DOWNLOADS_DIR=/path/to/downloads ."
echo "   3. Update your build_container.sh to use: --downloads-dir $DOWNLOAD_DIR"

# Exit with error if verification mode and any files failed
if [[ "$VERIFY_ONLY" == "true" ]]; then
    if [[ "$DOWNLOAD_FSL" == "true" && "$FSL_SUCCESS" == "false" ]] || \
       [[ "$DOWNLOAD_FREESURFER" == "true" && "$FREESURFER_SUCCESS" == "false" ]] || \
       [[ "$DOWNLOAD_AFNI" == "true" && "$AFNI_SUCCESS" == "false" ]] || \
       [[ "$DOWNLOAD_FSL_FIX" == "true" && "$FSL_FIX_SUCCESS" == "false" ]]; then
        exit 1
    fi
fi