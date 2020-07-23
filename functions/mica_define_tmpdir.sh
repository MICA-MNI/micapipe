#!/bin/bash
# Define a tmporal directory for mica-bids-proc

# Written by Reinder Vos de Wael (Oct 2020).

tmp=""
while getopts ":t:" opt; do
    case $opt in
        t)
            tmp=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done
​
if [[ -z $tmp ]]; then
    if [[ -z $TMP ]]; then
        tmp=$(mktmp -d)
    else
        tmp=$TMP
    fi
fi
