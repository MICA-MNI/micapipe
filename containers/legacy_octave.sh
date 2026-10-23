#!/bin/bash
# Load the packages required by the original FIX implementation.
set -euo pipefail
args=("$@")
for ((i=0;i<${#args[@]};i++)); do
    if [[ "${args[$i]}" == --eval ]]; then
        ((i+=1))
        args[$i]="warning('off','Octave:shadowed-function'); pkg load signal; pkg load statistics; ${args[$i]}"
    fi
done
exec /usr/bin/octave "${args[@]}"
