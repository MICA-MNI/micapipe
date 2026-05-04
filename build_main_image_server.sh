#!/bin/bash
# Build the fast-changing MICApipe image layer from a staged context.

set -euo pipefail

export DOCKER_CONTENT_TRUST=0

SERVER_ROOT="${MICAPIPE_SERVER_ROOT:-/export03/data/enning}"
STAGING_ROOT="${MICAPIPE_STAGING_ROOT:-$(cd .. && pwd)/staging}"
BASE_IMAGE="${MICAPIPE_BASE_IMAGE:-ghcr.io/mica-mni/micapipe-comprehensive-base:latest}"
REGISTRY="${MICAPIPE_REGISTRY:-ghcr.io/mica-mni}"
CUSTOM_TAG=""
PUSH_TO_REGISTRY="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --base-image)
            BASE_IMAGE="$2"
            shift 2
            ;;
        --tag)
            CUSTOM_TAG="$2"
            shift 2
            ;;
        --registry)
            REGISTRY="$2"
            shift 2
            ;;
        --push)
            PUSH_TO_REGISTRY="true"
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ "$PWD" != "$SERVER_ROOT"* ]]; then
    echo "Run this script from a server workspace under $SERVER_ROOT" >&2
    exit 1
fi

if ! docker info >/dev/null 2>&1; then
    echo "Cannot access Docker daemon" >&2
    exit 1
fi

if [[ ! -x "./prepare_build_context.sh" ]]; then
    echo "prepare_build_context.sh not found or not executable" >&2
    exit 1
fi

if docker image inspect "$BASE_IMAGE" >/dev/null 2>&1; then
    :
elif ! docker pull "$BASE_IMAGE" >/dev/null 2>&1; then
    echo "Base image not found: $BASE_IMAGE" >&2
    exit 1
fi

if [[ -z "$CUSTOM_TAG" ]]; then
    IMAGE_TAG="v0.2.3-$(date +%Y%m%d)"
else
    IMAGE_TAG="$CUSTOM_TAG"
fi

MAIN_IMAGE_NAME="micapipe"
FULL_MAIN_IMAGE="${REGISTRY}/${MAIN_IMAGE_NAME}:${IMAGE_TAG}"
LATEST_MAIN_IMAGE="${REGISTRY}/${MAIN_IMAGE_NAME}:latest"
mkdir -p "$STAGING_ROOT"
CONTEXT_DIR="$(./prepare_build_context.sh --target main --source "$PWD" --output-root "$STAGING_ROOT" --print-path)"
BUILD_LOG="build_main_${IMAGE_TAG}_$(date +%H%M%S).log"

echo "🚀 MICApipe Main Image Builder"
echo "📍 Workspace: $PWD"
echo "🎯 Base image: ${BASE_IMAGE}"
echo "🧱 Staged main context: ${CONTEXT_DIR}"
du -sh "$CONTEXT_DIR"
echo "📝 Build log: ${BUILD_LOG}"
echo

START_TIME=$(date +%s)
if docker build \
    --file "${CONTEXT_DIR}/Dockerfile.main" \
    --build-arg "BASE_IMAGE=${BASE_IMAGE}" \
    --tag "${FULL_MAIN_IMAGE}" \
    --tag "${LATEST_MAIN_IMAGE}" \
    "${CONTEXT_DIR}" 2>&1 | tee "$BUILD_LOG"; then
    END_TIME=$(date +%s)
    BUILD_DURATION=$((END_TIME - START_TIME))
    BUILD_MINUTES=$((BUILD_DURATION / 60))
    BUILD_SECONDS=$((BUILD_DURATION % 60))

    echo
    echo "✅ Main image build finished in ${BUILD_MINUTES}m ${BUILD_SECONDS}s"
    docker images "${REGISTRY}/${MAIN_IMAGE_NAME}" --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}\t{{.CreatedAt}}"

    if [[ "$PUSH_TO_REGISTRY" == "true" ]]; then
        docker push "${FULL_MAIN_IMAGE}"
        docker push "${LATEST_MAIN_IMAGE}"
    fi
else
    echo
    echo "❌ Main image build failed. Check ${BUILD_LOG}" >&2
    exit 1
fi
