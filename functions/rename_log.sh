#!/bin/bash
#
# MICAPIPE  function that renames multiple acquisitions log files

log_file_str="$1"
labelStr="$2"
if [ -f ${log_file_str}.txt ]; then
    tagMRI=$(grep ${labelStr} ${log_file_str}.txt | awk -F ':' 'NR==1 {print $2}'| tr -d '[:space:]' | sed -r "s/\x1B\[(([0-9]+)(;[0-9]+)*)?[m,K,H,f,J]//g")
    mv -v ${log_file_str}.txt ${log_file_str}_${tagMRI}.txt
    if [ -f ${log_file_str}.e ]; then mv ${log_file_str}.e ${log_file_str}_${tagMRI}.e; fi
fi
