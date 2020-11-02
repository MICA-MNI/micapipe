#!/bin/bash
#
# MICA pipe Quality Check script
#
# This script will create a basic html file for QC of the processing
#

#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#
# ONLY for scripting and debugging:
# TEST=ON

#---------------- FUNCTION: HELP ----------------#
help() {
echo -e "\033[38;5;141m
Usage:  `basename $0`\033[0m  \033[38;5;197m-sub\033[0m <subject_id> \033[38;5;197m-out\033[0m <outputDirectory> \033[38;5;197m-bids\033[0m <BIDS-directory>\n
\t\033[38;5;197m-sub\033[0m 	         : Subject identification (no 'sub-')
\t\033[38;5;197m-out\033[0m 	         : Output directory for the processed files <derivatives>.
\t\033[38;5;197m-bids\033[0m 	         : Path to BIDS directory

\t\033[38;5;197m-ses\033[0m 	         : OPTIONAL flag that indicates the session name (DEFAULT is ses-pre)


McGill University, MNI, MICA-lab, May-September 2020
https://github.com/MICA-MNI/micapipe
http://mica-mni.github.io/
"
}

# Check MICAPIPE
if [ -z "${MICAPIPE}" ]; then
echo -e "\033[38;5;1m\n---------------------------------------------------------\n
[ERROR]... MICAPIPE must be define in your enviroment\033[0m
           TRY: export MICAPIPE=<github_Directory>/micasoft\n
\033[38;5;1m---------------------------------------------------------\033[0m\n"; exit 1
fi

if [ ! -f ${MICAPIPE}/functions/utilities.sh ]; then
echo -e "\033[38;5;1m\n---------------------------------------------------------\n
[ERROR]... MICAPIPE is defined, but the PATH is wrong,
           it should match /micasoft directory\033[0m
           CHECK PATH to MICAPIPE:
           $MICAPIPE\n
\033[38;5;1m---------------------------------------------------------\033[0m\n"; exit 1
fi

# Source utilities functions from MICAPIPE
source ${MICAPIPE}/functions/utilities.sh

# -----------------------------------------------------------------------------------------------#
#			ARGUMENTS
# Create VARIABLES
for arg in "$@"
do
  case "$arg" in
  -h|-help)
    help
    exit 1
  ;;
  -sub)
    id=$2
    shift;shift
  ;;
  -out)
    out=$2
    shift;shift
  ;;
  -bids)
    BIDS=$2
    shift;shift
  ;;
  -ses)
    SES=$2
    shift
  ;;
  -*)
    Error "Unknown option ${2}"
    help
    exit 1
  ;;
    esac
done

# argument check out & WARNINGS
arg=($id $out $BIDS)
if [ "${#arg[@]}" -lt 3 ]; then
Error "One or more mandatory arguments are missing:"
Note "-sub " $id
Note "-out " "$out"
Note "-bids " "$BIDS"
help; exit 1; fi


# Get the real path of the Inputs
id=${id/sub-/}
out=`realpath $out`
BIDS=`realpath $BIDS`

# Number of session (Default is "ses-pre")
if [ -z ${SES} ]; then SES="ses-pre"; else SES="ses-${SES/ses-/}"; fi
# Exit if subject is not found
if [ ! -d ${out}/sub-${id}/${SES} ]; then Error "$id was not found on the OUTPUT directory\n\t Check ls ${out}/sub-${id}/${SES}"; exit 1; fi

# Assigns variables names
bids_variables $BIDS $id $out $SES
