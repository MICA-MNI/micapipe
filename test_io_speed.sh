#!/bin/bash
#
# Test I/O speed of different drives
# Usage: ./test_io_speed.sh /path/to/test1 /path/to/test2 ...
#

set -e

# Default test size (16MB)
TEST_SIZE_MB=16
TEST_FILE="io_test_file_$$.tmp"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log() { echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $*"; }
success() { echo -e "${GREEN}✅ $*${NC}"; }
warn() { echo -e "${YELLOW}⚠️  $*${NC}"; }
error() { echo -e "${RED}❌ $*${NC}"; }

# Function to test a single path
test_path() {
    local TEST_PATH="$1"
    local TEST_NAME="${2:-$TEST_PATH}"
    
    echo ""
    echo "============================================="
    echo -e "${BLUE}Testing: ${TEST_NAME}${NC}"
    echo "Path: ${TEST_PATH}"
    echo "============================================="
    
    # Check if path exists
    if [[ ! -d "$TEST_PATH" ]]; then
        error "Path does not exist: $TEST_PATH"
        return 1
    fi
    
    # Check if writable
    if [[ ! -w "$TEST_PATH" ]]; then
        error "Path is not writable: $TEST_PATH"
        return 1
    fi
    
    local FULL_TEST_FILE="${TEST_PATH}/${TEST_FILE}"
    
    # Get filesystem info
    echo ""
    echo "📁 Filesystem Info:"
    df -h "$TEST_PATH" | head -2
    
    # Check mount type
    echo ""
    echo "🔗 Mount Type:"
    MOUNT_INFO=$(df "$TEST_PATH" | tail -1 | awk '{print $1}')
    if [[ "$MOUNT_INFO" == *":"* ]] || [[ "$MOUNT_INFO" == *"nfs"* ]]; then
        warn "This appears to be a NETWORK mount (NFS/CIFS)"
    elif [[ "$MOUNT_INFO" == "/dev/"* ]]; then
        success "This appears to be a LOCAL drive"
    else
        echo "Mount device: $MOUNT_INFO"
    fi
    
    # ============================================
    # Write Test
    # ============================================
    echo ""
    echo "📤 WRITE TEST (${TEST_SIZE_MB}MB)..."
    
    # Clear cache (if possible)
    sync
    
    # Write test with dd
    WRITE_START=$(date +%s.%N)
    dd if=/dev/zero of="$FULL_TEST_FILE" bs=1M count=$TEST_SIZE_MB conv=fdatasync 2>&1 | tail -1
    WRITE_END=$(date +%s.%N)
    
    WRITE_TIME=$(echo "$WRITE_END - $WRITE_START" | bc)
    WRITE_SPEED=$(echo "scale=2; $TEST_SIZE_MB / $WRITE_TIME" | bc)
    
    echo -e "   Write Speed: ${GREEN}${WRITE_SPEED} MB/s${NC}"
    echo "   Write Time: ${WRITE_TIME}s"
    
    # ============================================
    # Read Test (with cache cleared)
    # ============================================
    echo ""
    echo "📥 READ TEST (${TEST_SIZE_MB}MB)..."
    
    # Try to clear cache
    sync
    if [[ -w /proc/sys/vm/drop_caches ]]; then
        echo 3 > /proc/sys/vm/drop_caches 2>/dev/null || true
    fi
    
    # Read test with dd
    READ_START=$(date +%s.%N)
    dd if="$FULL_TEST_FILE" of=/dev/null bs=1M 2>&1 | tail -1
    READ_END=$(date +%s.%N)
    
    READ_TIME=$(echo "$READ_END - $READ_START" | bc)
    READ_SPEED=$(echo "scale=2; $TEST_SIZE_MB / $READ_TIME" | bc)
    
    echo -e "   Read Speed: ${GREEN}${READ_SPEED} MB/s${NC}"
    echo "   Read Time: ${READ_TIME}s"
    
    # ============================================
    # Random I/O Test (small files)
    # ============================================
    echo ""
    echo "🎲 RANDOM I/O TEST (1000 x 4KB files)..."
    
    RANDOM_DIR="${TEST_PATH}/io_test_random_$$"
    mkdir -p "$RANDOM_DIR"
    
    RANDOM_START=$(date +%s.%N)
    for i in $(seq 1 1000); do
        dd if=/dev/urandom of="${RANDOM_DIR}/file_${i}" bs=4K count=1 2>/dev/null
    done
    sync
    RANDOM_END=$(date +%s.%N)
    
    RANDOM_TIME=$(echo "$RANDOM_END - $RANDOM_START" | bc)
    IOPS=$(echo "scale=2; 1000 / $RANDOM_TIME" | bc)
    
    echo -e "   IOPS (approx): ${GREEN}${IOPS} ops/s${NC}"
    echo "   Time for 1000 files: ${RANDOM_TIME}s"
    
    # Cleanup random files
    rm -rf "$RANDOM_DIR"
    
    # ============================================
    # Latency Test
    # ============================================
    echo ""
    echo "⏱️  LATENCY TEST (single 4KB write)..."
    
    LATENCY_TOTAL=0
    for i in $(seq 1 10); do
        LAT_START=$(date +%s.%N)
        dd if=/dev/zero of="${FULL_TEST_FILE}.lat" bs=4K count=1 conv=fdatasync 2>/dev/null
        LAT_END=$(date +%s.%N)
        LAT=$(echo "$LAT_END - $LAT_START" | bc)
        LATENCY_TOTAL=$(echo "$LATENCY_TOTAL + $LAT" | bc)
        rm -f "${FULL_TEST_FILE}.lat"
    done
    
    AVG_LATENCY=$(echo "scale=4; $LATENCY_TOTAL / 10 * 1000" | bc)
    echo -e "   Avg Latency: ${GREEN}${AVG_LATENCY} ms${NC}"
    
    # Cleanup
    rm -f "$FULL_TEST_FILE"
    
    # ============================================
    # Summary
    # ============================================
    echo ""
    echo "📊 SUMMARY for ${TEST_NAME}:"
    echo "   Sequential Write: ${WRITE_SPEED} MB/s"
    echo "   Sequential Read:  ${READ_SPEED} MB/s"
    echo "   Random IOPS:      ${IOPS} ops/s"
    echo "   Avg Latency:      ${AVG_LATENCY} ms"
    
    # Classification
    echo ""
    if (( $(echo "$WRITE_SPEED > 500" | bc -l) )); then
        success "FAST - Likely LOCAL SSD/NVMe or fast RAID"
    elif (( $(echo "$WRITE_SPEED > 100" | bc -l) )); then
        warn "MODERATE - Could be local HDD or fast network"
    else
        error "SLOW - Likely NETWORK mounted drive"
    fi
    
    # Return results for comparison
    echo "${TEST_NAME}|${WRITE_SPEED}|${READ_SPEED}|${IOPS}|${AVG_LATENCY}"
}

# ============================================
# Main
# ============================================

echo ""
echo "🔬 I/O SPEED TESTER"
echo "==================="
echo "Test size: ${TEST_SIZE_MB}MB"
echo ""

# Default paths if none provided
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 /path1 [/path2] [/path3] ..."
    echo ""
    echo "Example paths to test:"
    echo "  $0 /export03/data/enning ~/micapipe /tmp"
    echo ""
    echo "Testing default paths..."
    
    PATHS_TO_TEST=(
        "/export03/data/enning|Server Export"
        "$HOME|Home Directory"
        "/tmp|Temp (usually local)"
    )
else
    PATHS_TO_TEST=()
    for path in "$@"; do
        PATHS_TO_TEST+=("$path|$path")
    done
fi

# Store results for comparison
declare -a RESULTS

for path_info in "${PATHS_TO_TEST[@]}"; do
    IFS='|' read -r path name <<< "$path_info"
    if [[ -d "$path" ]]; then
        RESULT=$(test_path "$path" "$name" | tail -1)
        RESULTS+=("$RESULT")
    else
        warn "Skipping non-existent path: $path"
    fi
done

# ============================================
# Comparison Table
# ============================================
echo ""
echo ""
echo "============================================="
echo "📊 COMPARISON TABLE"
echo "============================================="
printf "%-25s %12s %12s %12s %12s\n" "PATH" "WRITE MB/s" "READ MB/s" "IOPS" "LATENCY ms"
echo "---------------------------------------------"

for result in "${RESULTS[@]}"; do
    IFS='|' read -r name write read iops latency <<< "$result"
    printf "%-25s %12s %12s %12s %12s\n" "$name" "$write" "$read" "$iops" "$latency"
done

echo ""
echo "============================================="
echo "🏆 RECOMMENDATIONS:"
echo "============================================="
echo "• For Docker builds: Use the path with highest WRITE speed"
echo "• For Singularity: Use the path with most FREE SPACE + good speed"
echo "• For temp files: Use /tmp or fastest local drive"
echo ""
