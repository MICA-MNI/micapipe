#!/usr/bin/env bash
#
# Unified Singularity/Apptainer SIF builder for micapipe that works on
# workstations and HPC nodes without sudo access.
#
# Features:
#   * Auto-detects whether to use local Docker daemon, docker-archive tarball,
#     or a remote OCI registry source (GHCR by default)
#   * Works with either `apptainer` or `singularity`
#   * Non-interactive (safe for batch jobs)
#   * Stores cache/tmp data outside $HOME to avoid quota limits
#   * Handles nodev mountpoints by staging the build in a safe location
#
# Usage examples:
#   ./build_singularity_auto.sh                   # auto, uses ghcr.io/mica-mni/micapipe:latest
#   ./build_singularity_auto.sh --tag v0.2.3      # explicit tag (local docker preferred)
#   ./build_singularity_auto.sh --oci docker://ghcr.io/mica-mni/micapipe:v0.2.3
#   ./build_singularity_auto.sh --tar /path/micapipe.tar --output /data/micapipe.sif
#   ./build_singularity_auto.sh --force --output-dir /data_/mica1/01_programs/singularity
#
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
DOCKER_IMAGE_DEFAULT="ghcr.io/mica-mni/micapipe"
DEFAULT_TAG="latest"
DEFAULT_OUTPUT_DIR="${SINGULARITY_DIR:-/host/cassio/export03/data/enning/singularity}"
DEFAULT_OUTPUT_NAME="micapipe_v1_beta.sif"
CACHE_ROOT_DEFAULT="${CACHE_ROOT:-${DEFAULT_OUTPUT_DIR}}"

DOCKER_IMAGE="$DOCKER_IMAGE_DEFAULT"
DOCKER_TAG="$DEFAULT_TAG"
OUTPUT_DIR="$DEFAULT_OUTPUT_DIR"
OUTPUT_NAME="$DEFAULT_OUTPUT_NAME"
OUTPUT_PATH=""
SOURCE_MODE="auto"         # auto | docker-daemon | docker-archive | oci
OCI_URI=""
TAR_FILE=""
FORCE_OVERWRITE=0
DEBUG_FLAG=0

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [TAG] [options]

Options:
  --tag TAG              Docker tag (default: ${DEFAULT_TAG})
  --image IMAGE          Docker image name (default: ${DOCKER_IMAGE_DEFAULT})
  --oci URI              Build/pull from OCI registry URI (forces source=oci)
  --source MODE          Source mode: auto | docker-daemon | docker-archive | oci
  --tar FILE             Use docker-archive FILE (forces source=docker-archive)
  --output FILE          Output SIF path (default: <output-dir>/${DEFAULT_OUTPUT_NAME})
  --output-dir DIR       Destination directory (default: ${DEFAULT_OUTPUT_DIR})
  --force                Overwrite existing output without prompt (backs up old file)
  --debug                Enable verbose debug mode on the container CLI
  --help                 Show this help message

Environment overrides:
  SINGULARITY_DIR        Default output directory (same as --output-dir)
  CACHE_ROOT             Base directory for caches/temp (default: output dir)
EOF
}

log() { printf '%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
fail() { echo "❌ $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tag)
            [[ $# -lt 2 ]] && fail "Missing value for --tag"
            DOCKER_TAG="$2"; shift 2 ;;
        --image)
            [[ $# -lt 2 ]] && fail "Missing value for --image"
            DOCKER_IMAGE="$2"; shift 2 ;;
        --oci)
            [[ $# -lt 2 ]] && fail "Missing value for --oci"
            OCI_URI="$2"; SOURCE_MODE="oci"; shift 2 ;;
        --source)
            [[ $# -lt 2 ]] && fail "Missing value for --source"
            SOURCE_MODE="$2"; shift 2 ;;
        --tar)
            [[ $# -lt 2 ]] && fail "Missing value for --tar"
            TAR_FILE="$2"; SOURCE_MODE="docker-archive"; shift 2 ;;
        --output)
            [[ $# -lt 2 ]] && fail "Missing value for --output"
            OUTPUT_PATH="$2"; shift 2 ;;
        --output-dir)
            [[ $# -lt 2 ]] && fail "Missing value for --output-dir"
            OUTPUT_DIR="$2"; shift 2 ;;
        --force)
            FORCE_OVERWRITE=1; shift ;;
        --debug)
            DEBUG_FLAG=1; shift ;;
        --help|-h)
            usage; exit 0 ;;
        --)
            shift; break ;;
        -* )
            fail "Unknown option: $1" ;;
        *)
            if [[ "$DOCKER_TAG" == "$DEFAULT_TAG" && "$SOURCE_MODE" == "auto" ]]; then
                DOCKER_TAG="$1"; shift
            else
                fail "Unexpected argument: $1"
            fi ;;
    esac
done

CACHE_ROOT="${CACHE_ROOT:-$CACHE_ROOT_DEFAULT}"
FULL_DOCKER_IMAGE="${DOCKER_IMAGE}:${DOCKER_TAG}"
[[ -z "$OUTPUT_PATH" ]] && OUTPUT_PATH="${OUTPUT_DIR}/${OUTPUT_NAME}"
[[ -z "$OCI_URI" ]] && OCI_URI="docker://${FULL_DOCKER_IMAGE}"

log "🔄 Source image: ${FULL_DOCKER_IMAGE}"
log "📍 Output path: ${OUTPUT_PATH}"
log "🗂️  Cache root: ${CACHE_ROOT}"

# ---------------------------------------------------------------------------
# Detect container CLI (apptainer preferred, fallback to singularity)
# ---------------------------------------------------------------------------
if command -v apptainer >/dev/null 2>&1; then
    SING_CMD="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    SING_CMD="singularity"
else
    fail "Neither apptainer nor singularity found in PATH"
fi

log "🛠️  Using CLI: ${SING_CMD} $(${SING_CMD} --version 2>/dev/null)"

# ---------------------------------------------------------------------------
# Prepare directories and env
# ---------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$CACHE_ROOT"
CACHE_DIR="${CACHE_ROOT}/.singularity_cache"
TMP_DIR="${CACHE_ROOT}/.singularity_tmp"
mkdir -p "$CACHE_DIR" "$TMP_DIR"

export SINGULARITY_CACHEDIR="$CACHE_DIR"
export SINGULARITY_TMPDIR="$TMP_DIR"
export APPTAINER_CACHEDIR="$CACHE_DIR"
export APPTAINER_TMPDIR="$TMP_DIR"

log "📦 Cache dir: ${CACHE_DIR}"
log "📁 Temp dir:  ${TMP_DIR}"

if command -v df >/dev/null 2>&1; then
    log "💾 Space @output: $(df -h "$OUTPUT_DIR" | awk 'NR==2 {print $4}')"
    log "💾 Space @cache:  $(df -h "$CACHE_ROOT" | awk 'NR==2 {print $4}')"
fi

# ---------------------------------------------------------------------------
# Decide source mode
# ---------------------------------------------------------------------------
if [[ "$SOURCE_MODE" == "auto" ]]; then
    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1 && \
       docker image inspect "$FULL_DOCKER_IMAGE" >/dev/null 2>&1; then
        SOURCE_MODE="docker-daemon"
    elif [[ -n "$TAR_FILE" && -f "$TAR_FILE" ]]; then
        SOURCE_MODE="docker-archive"
    else
        SOURCE_MODE="oci"
    fi
fi

case "$SOURCE_MODE" in
    docker-daemon)
        log "📚 Source mode: docker-daemon"
        docker image inspect "$FULL_DOCKER_IMAGE" >/dev/null 2>&1 || \
            fail "Docker image ${FULL_DOCKER_IMAGE} not present locally"
        BUILD_SOURCE="docker-daemon://${FULL_DOCKER_IMAGE}"
        ;;
    docker-archive)
        log "📚 Source mode: docker-archive"
        if [[ -z "$TAR_FILE" ]]; then
            TAR_FILE="${CACHE_ROOT}/micapipe_${DOCKER_TAG}_$(date +%Y%m%d_%H%M%S).tar"
            log "📦 Exporting docker image to ${TAR_FILE}"
            docker image inspect "$FULL_DOCKER_IMAGE" >/dev/null 2>&1 || \
                fail "Docker image ${FULL_DOCKER_IMAGE} not present locally"
            docker save "$FULL_DOCKER_IMAGE" -o "$TAR_FILE"
            CLEANUP_TAR=1
        else
            [[ -f "$TAR_FILE" ]] || fail "Specified tar file not found: ${TAR_FILE}"
            CLEANUP_TAR=0
        fi
        BUILD_SOURCE="docker-archive://${TAR_FILE}"
        ;;
    oci)
        log "📚 Source mode: oci"
        BUILD_SOURCE="$OCI_URI"
        ;;
    *)
        fail "Unknown source mode: ${SOURCE_MODE}"
        ;;
esac

# ---------------------------------------------------------------------------
# Handle pre-existing output
# ---------------------------------------------------------------------------
if [[ -f "$OUTPUT_PATH" ]]; then
    if [[ $FORCE_OVERWRITE -eq 1 ]]; then
        log "♻️  Removing existing output (force): ${OUTPUT_PATH}"
        rm -f "$OUTPUT_PATH"
    else
        BACKUP_PATH="${OUTPUT_PATH%.sif}_$(date +%Y%m%d_%H%M%S).bak.sif"
        log "ℹ️  Existing output detected → moving to ${BACKUP_PATH}"
        mv "$OUTPUT_PATH" "$BACKUP_PATH"
    fi
fi

# Choose safe temp location (handles nodev mounts)
OUTPUT_DEVICE=$(df "$OUTPUT_DIR" | awk 'NR==2 {print $1}')
if mount | grep -q "${OUTPUT_DEVICE}.*nodev"; then
    TMP_OUTPUT="${CACHE_ROOT}/micapipe_build_$$.sif"
    log "⚠️  nodev filesystem detected → using staged output ${TMP_OUTPUT}"
else
    TMP_OUTPUT="${OUTPUT_PATH}.tmp"
fi

CLEANUP_LIST=()
cleanup() {
    local status=$?
    if [[ $status -ne 0 ]]; then
        log "⚠️  Build aborted (exit code ${status})"
    fi
    for item in "${CLEANUP_LIST[@]}"; do
        [[ -f "$item" ]] && rm -f "$item"
    done
    if [[ "${CLEANUP_TAR:-0}" -eq 1 && -f "$TAR_FILE" ]]; then
        log "🧹 Removing temporary tar ${TAR_FILE}"
        rm -f "$TAR_FILE"
    fi
}
trap cleanup EXIT

CLEANUP_LIST+=("$TMP_OUTPUT")

log "🚀 Starting build (mode=${SOURCE_MODE})"
START_TIME=$(date +%s)

cmd=("${SING_CMD}")
if [[ "$SOURCE_MODE" == "oci" ]]; then
    cmd+=(pull --force)
else
    cmd+=(build --force)
fi
[[ $DEBUG_FLAG -eq 1 ]] && cmd+=(--debug)
cmd+=("$TMP_OUTPUT" "$BUILD_SOURCE")

log "🔧 Command: ${cmd[*]}"
if ! "${cmd[@]}"; then
    fail "Container CLI reported a failure"
fi

mv "$TMP_OUTPUT" "$OUTPUT_PATH"
trap - EXIT
cleanup

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
log "✅ Build complete in $((DURATION/60))m $((DURATION%60))s"
log "📦 Output size: $(du -h "$OUTPUT_PATH" | cut -f1)"

if "${SING_CMD}" inspect "$OUTPUT_PATH" >/dev/null 2>&1; then
    log "🔍 Inspect OK"
else
    log "⚠️  Inspect reported issues (image may still be usable)"
fi

echo ""
echo "Next steps:"
echo "  ${SING_CMD} run ${OUTPUT_PATH} --help"
echo "  ${SING_CMD} inspect ${OUTPUT_PATH}"
echo "  ${SING_CMD} exec ${OUTPUT_PATH} micapipe --help"
