#!/bin/bash
set -euo pipefail
if [[ "${MICAPIPE_FIX_BACKEND:-legacy}" == legacy ]]; then
    export FSL_FIXDIR=/opt/fix-legacy FSL_FIX_R_CMD=/usr/bin/R
    exec /opt/fix-legacy/fix "$@"
fi
logdir="${TMPDIR:-/tmp}"
has_log=0
has_jobs=0
for arg in "$@"; do
    [[ "$arg" == -lf || "$arg" == --logfile ]] && has_log=1
    [[ "$arg" == -nj || "$arg" == --n_jobs ]] && has_jobs=1
    if [[ -d "$arg" ]]; then logdir="$arg"; fi
done
extra=()
if [[ "$has_log" == 0 ]]; then extra=(-lf "$logdir/pyfix.log"); fi
if [[ "$has_jobs" == 0 ]]; then extra+=(-nj "${MICAPIPE_FIX_JOBS:-${threads:-1}}"); fi
# pyFIX 0.6 treats its first argument as a stage selector: keep it first.
exec /opt/fsl-6.0.7.23/bin/fix "$@" "${extra[@]}"
