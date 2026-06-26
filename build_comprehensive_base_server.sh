#!/bin/bash
# Build the comprehensive MICApipe base image from a staged context.

set -euo pipefail

export DOCKER_CONTENT_TRUST=0
export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain

SERVER_ROOT="${MICAPIPE_SERVER_ROOT:-/export03/data/enning}"
DOWNLOADS_SOURCE="${MICAPIPE_DOWNLOADS_SOURCE:-$SERVER_ROOT/downloads}"
STAGING_ROOT="${MICAPIPE_STAGING_ROOT:-$(cd .. && pwd)/staging}"
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
BASE_IMAGE_NAME="micapipe-comprehensive-base"
BASE_IMAGE_TAG="${MICAPIPE_BASE_TAG:-$(date +%Y%m%d)}"
FULL_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:${BASE_IMAGE_TAG}"
LATEST_BASE_IMAGE="${REGISTRY}/${BASE_IMAGE_NAME}:latest"

REQUIRED_DOWNLOADS=(
    "fsl-6.0.2-centos6_64.tar.gz"
    "freesurfer-linux-ubuntu18_amd64-7.4.1.tar.gz"
    "afni-linux_openmp_64.tgz"
    "fix-1.068.tar.gz"
)

if [[ "$PWD" != "$SERVER_ROOT"* ]]; then
    echo "Run this script from a server workspace under $SERVER_ROOT" >&2
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "Cannot access Docker daemon" >&2
    exit 1
fi

if ! docker buildx version >/dev/null 2>&1; then
    echo "docker buildx is required for mounted server downloads" >&2
    exit 1
fi

if [[ ! -x "./prepare_build_context.sh" ]]; then
    echo "prepare_build_context.sh not found or not executable" >&2
    exit 1
fi

echo "🐳 MICApipe Comprehensive Base Image Builder"
echo "📍 Workspace: $PWD"
echo "🎯 Target image: ${FULL_BASE_IMAGE}"
echo "📦 Archive source: ${DOWNLOADS_SOURCE}"
echo
echo "Archive availability:"
for file in "${REQUIRED_DOWNLOADS[@]}"; do
    if [[ -f "$PWD/$file" ]]; then
        echo "  ✅ ${file} (workspace)"
    elif [[ -f "${DOWNLOADS_SOURCE}/${file}" ]]; then
        echo "  ✅ ${file} (${DOWNLOADS_SOURCE})"
    else
        echo "  ⚠️  ${file} missing, Dockerfile will download it"
    fi
done

CONTEXT_DIR="$(./prepare_build_context.sh --target base --source "$PWD" --downloads-source "$DOWNLOADS_SOURCE" --output-root "$STAGING_ROOT" --print-path)"
BUILD_LOG="build_comprehensive_base_$(date +%Y%m%d_%H%M%S).log"
export DOCKER_TMPDIR="${SERVER_ROOT}/docker_tmp"
mkdir -p "$DOCKER_TMPDIR"
mkdir -p "$STAGING_ROOT"

echo
echo "🧱 Staged base context: ${CONTEXT_DIR}"
du -sh "$CONTEXT_DIR"
echo "📦 Mounted downloads dir: ${DOWNLOADS_SOURCE}"
echo "🗂️  Docker temp dir: ${DOCKER_TMPDIR}"
echo "📝 Build log: ${BUILD_LOG}"
echo

CACHE_FROM_IMAGES=(
    "micapipe-comprehensive-base:latest"
    "ubuntu:bionic-20201119"
    "kaczmarj/ants:2.3.4"
)

CACHE_ARGS=()
for image in "${CACHE_FROM_IMAGES[@]}"; do
    CACHE_ARGS+=(--cache-from "$image")
done

if docker buildx build --load \
    --file "${CONTEXT_DIR}/Dockerfile.mamba-base" \
    --build-context "downloads=${DOWNLOADS_SOURCE}" \
    --build-arg CUSTOM_TMPDIR="${SERVER_ROOT}" \
    "${CACHE_ARGS[@]}" \
    --tag "${FULL_BASE_IMAGE}" \
    --tag "${LATEST_BASE_IMAGE}" \
    "${CONTEXT_DIR}" 2>&1 | tee "$BUILD_LOG"; then
    echo
    echo "✅ Base image build finished"
    docker images "${REGISTRY}/${BASE_IMAGE_NAME}" --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"
else
    echo
    echo "❌ Base image build failed. Check ${BUILD_LOG}" >&2
    exit 1
fi
