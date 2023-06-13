#!/bin/bash
#
# MICAPIPE  function that renames multiple acquisitions log files

log_file_str="$1"
labelStr="$2"
randStr="$3"
if [ -f "${log_file_str}".txt ]; then
    new_log_str="${log_file_str/-${randStr}/}"
    tagMRI=$(grep "${labelStr}" "${log_file_str}".txt | awk -F ':' 'NR==1 {print $2}'| tr -d '[:space:]' | sed -r "s/\x1B\[(([0-9]+)(;[0-9]+)*)?[m,K,H,f,J]//g")
    mv "${log_file_str}".txt "${new_log_str}_${tagMRI}".txt
    rm "${new_log_str}"*_.txt 2>/dev/null
    if [ -f "${log_file_str}".e ]; then mv "${log_file_str}".e "${new_log_str}_${tagMRI}".e; fi
    rm "${new_log_str}"*_.e 2>/dev/null
fi
