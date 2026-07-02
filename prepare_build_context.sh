#!/bin/bash
#
# Prepare explicit Docker build contexts for MICApipe server experiments.
# The base context keeps static data only.
# Large installer archives stay in the server downloads directory and are
# mounted into the Docker build directly to avoid re-copying multi-GB files.
# The main context keeps only frequently changing code.
#

set -euo pipefail

TARGET=""
SOURCE_DIR="${MICAPIPE_SOURCE_DIR:-$PWD}"
OUTPUT_ROOT=""
DOWNLOADS_SOURCE="${MICAPIPE_DOWNLOADS_SOURCE:-/export03/data/enning/downloads}"
PRINT_PATH_ONLY="false"

usage() {
    cat <<'EOF'
Usage: ./prepare_build_context.sh --target {base|main} [options]

Options:
  --source DIR          Source repository/workspace (default: current directory)
  --output-root DIR     Output root for staged contexts (default: SOURCE/.build-context)
  --downloads-source    Directory containing pre-downloaded installers
  --print-path          Print only the staged context path
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --target)
            TARGET="$2"
            shift 2
            ;;
        --source)
            SOURCE_DIR="$2"
            shift 2
            ;;
        --output-root)
            OUTPUT_ROOT="$2"
            shift 2
            ;;
        --downloads-source)
            DOWNLOADS_SOURCE="$2"
            shift 2
            ;;
        --print-path)
            PRINT_PATH_ONLY="true"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$TARGET" ]]; then
    echo "--target is required" >&2
    usage >&2
    exit 1
fi

if [[ ! -d "$SOURCE_DIR" ]]; then
    echo "Source directory not found: $SOURCE_DIR" >&2
    exit 1
fi

if [[ -z "$OUTPUT_ROOT" ]]; then
    OUTPUT_ROOT="${SOURCE_DIR}/.build-context"
fi

log() {
    if [[ "$PRINT_PATH_ONLY" != "true" ]]; then
        echo "$@"
    fi
}

copy_file() {
    local rel="$1"
    local dest_root="$2"

    if [[ ! -f "${SOURCE_DIR}/${rel}" ]]; then
        echo "Required file missing: ${SOURCE_DIR}/${rel}" >&2
        exit 1
    fi

    mkdir -p "${dest_root}/$(dirname "$rel")"
    cp "${SOURCE_DIR}/${rel}" "${dest_root}/${rel}"
}

copy_dir() {
    local rel="$1"
    local dest_root="$2"

    if [[ ! -d "${SOURCE_DIR}/${rel}" ]]; then
        echo "Required directory missing: ${SOURCE_DIR}/${rel}" >&2
        exit 1
    fi

    cp -R "${SOURCE_DIR}/${rel}" "${dest_root}/${rel}"
}

stage_base_context() {
    local context_dir="${OUTPUT_ROOT}/base"
    rm -rf "$context_dir"
    mkdir -p "$context_dir/functions"

    copy_file "Dockerfile.mamba-base" "$context_dir"
    copy_file "fix_settings.sh" "$context_dir"
    copy_file "functions/MICAMTL_training_15HC_15PX.RData" "$context_dir"
    copy_dir "R_config" "$context_dir"
    copy_dir "parcellations" "$context_dir"
    copy_dir "surfaces" "$context_dir"
    copy_dir "MNI152Volumes" "$context_dir"
    copy_dir "MICs60_T1-atlas" "$context_dir"
    copy_dir "fsl_conf" "$context_dir"

    if [[ "$PRINT_PATH_ONLY" == "true" ]]; then
        printf '%s\n' "$context_dir"
    else
        log "Prepared base build context: ${context_dir}"
        du -sh "$context_dir" | awk '{print "Context size: " $1}'
        printf '%s\n' "$context_dir"
    fi
}

stage_main_context() {
    local context_dir="${OUTPUT_ROOT}/main"
    rm -rf "$context_dir"
    mkdir -p "$context_dir"

    copy_file "Dockerfile.main" "$context_dir"
    copy_file "micapipe" "$context_dir"
    copy_file "micapipe.py" "$context_dir"
    copy_dir "functions" "$context_dir"

    rm -f "${context_dir}/functions/MICAMTL_training_15HC_15PX.RData"

    if [[ "$PRINT_PATH_ONLY" == "true" ]]; then
        printf '%s\n' "$context_dir"
    else
        log "Prepared main build context: ${context_dir}"
        du -sh "$context_dir" | awk '{print "Context size: " $1}'
        printf '%s\n' "$context_dir"
    fi
}

case "$TARGET" in
    base)
        stage_base_context
        ;;
    main)
        stage_main_context
        ;;
    *)
        echo "Unknown target: $TARGET" >&2
        exit 1
        ;;
esac
