#!/bin/bash
#
# DWI structural processing with bash:
#
# Preprocessing workflow for diffusion MRI.
#
# This workflow makes use of MRtrix3
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out parcDirectory
#   $4 : Temporal directory (default /tmp)
#
# ONLY for scripting and debugging
#TEST=ON
#

# Source utilities functions from MICAPIPE
source ${MICAPIPE}/functions/utilities.sh

function usage() {
    echo "MICA Diffusion Processing."
    echo "mica_diffusionProcessing [OPTION] finalDirectory warpDirectory algorithm subject"
    echo ""
    echo Options:
    echo "    -f          freesurfer directory."
    echo "    -d DWI      DWI scans, can be supplied multiple times for multishell processing"
    echo "    -r DWI      reverse phase DWI scan, optional"
    echo "    -h          bring up this help menu"
    echo "    -n          enable network processing (requires -p)"
    echo "    -p          path to processing directory, can be set to SETDEFAULT for default location"
    echo "    -b          debug mode; do not delete temporary directory"
}

export SUBJECTS_DIR=/

debug=false
while getopts ":d:r:f:np:b:" opt; do
    case $opt in
        d)
            dwiFiles+=("$OPTARG") ;;
        f)
            fsDirectory="$OPTARG" ;;
        r)
            reversePhaseScan="$OPTARG" ;;
        n)
            networkProcessing=true ;;
        p)
            processingDirectoryArg="$OPTARG" ;;
        b)
            debug=true ;;
        h)
            usage
            exit 1 ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1 ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1 ;;
    esac
done
shift $((OPTIND-1))
# Set basic parameters.
finalDirectory=$1
warpDirectory=$2
dwi2fodAlgorithm=$3
subject=$4

if [[ ! $dwi2fodAlgorithm == msmt_csd ]] && [[ ! $dwi2fodAlgorithm == csd ]]; then
    echo "Error: dwi2fodAlgorithm must be msmt_csd or csd."
    exit 1
fi

if [[ ! -z $networkProcessing ]] && [[ -z $processingDirectoryArg ]]; then
    echo "ERROR: If using network processing, processing directory must be specified."
    echo "Consider setting -p SETDEFAULT"
    exit
fi

[[ -z $NSLOTS ]] && NSLOTS=0
[[ -z $processingDirectoryArg ]] && processingDirectory=$(mktemp -d -p $finalDirectory) || processingDirectory=$(setProcessingDirectory $processingDirectoryArg)

# TRAP: Cleanup processing directory on exit.
function cleanup {
    # Clear temporary files. Add some checks to be sure there's little chance of unintentional deletes.
    if [[ $networkProcessing == true ]]; then
        find "$processingDirectory"/network/finalDirectory -maxdepth 1 -regextype posix-extended -regex ".*[.txt|.mif|.mif.gz]" -delete
        find "$processingDirectory"/network/warpDirectory/ -maxdepth 1 -regextype posix-extended -regex ".*[.dat|.dat~|.log|.mincost|.param|.sum|.lta]" -delete
        find "$processingDirectory"/network -maxdepth 1 -regextype posix-extended -regex ".*[.json|.bval|.bvec|.nii.gz|.nii|Directory]" -delete
        find "$processingDirectory" -empty -delete
    else
        find "$processingDirectory" -mindepth 1 -maxdepth 1 -type f -regextype posix-extended -regex ".*[.txt|.mif|.nii.gz|.bvec|.bval]" -delete
        find "$processingDirectory" -mindepth 0 -maxdepth 0 -type d -empty -delete
    fi
}
if [[ $debug == false ]] ; then
    trap cleanup EXIT
fi

# Set up network processing.
if [[ $networkProcessing == true ]]; then
    finalDirectoryReal=$finalDirectory
    finalDirectory=$processingDirectory/network/finalDirectory
    warpDirectoryReal=$warpDirectory
    warpDirectory=$processingDirectory/network/warpDirectory

    mkdir -p $processingDirectory/network/
    for i in $(seq 1 ${#dwiFiles[@]}); do
        cp -v ${dwiFiles[$i-1]}               $processingDirectory/network/$(basename ${dwiFiles[$i-1]})
        cp -v ${dwiFiles[$i-1]/.nii.gz/.bvec} $processingDirectory/network/$(basename ${dwiFiles[$i-1]/.nii.gz/.bvec})
        cp -v ${dwiFiles[$i-1]/.nii.gz/.bval} $processingDirectory/network/$(basename ${dwiFiles[$i-1]/.nii.gz/.bval})
        cp -v ${dwiFiles[$i-1]/.nii.gz/.json} $processingDirectory/network/$(basename ${dwiFiles[$i-1]/.nii.gz/.json})
        dwiFiles[$i-1]=$processingDirectory/network/$(basename ${dwiFiles[$i-1]})
    done
fi
for x in "$finalDirectory" "$warpDirectory" ; do
    [[ ! -d "$x" ]] && mkdir -p "$x"
done

# Concatenate shells -if only one shell then just convert to mif and rename.
if [[ ! -f "$finalDirectory"/"$subject"_dwi_residuals.mif.gz ]] && [[ ! -f "$finalDirectoryReal"/"$subject"_dwi_residuals.mif.gz ]] ; then
    if [[ ${#dwiFiles[@]} -gt 1 ]]; then
        do_cmd mrcat -force -nthreads $NSLOTS $(echo ${dwiFiles[@]}   | sort)   "$processingDirectory"/concatenateDwi.mif
        paste -d " " $(echo ${dwiFiles[@]/.nii.gz/.bvec} | sort) > "$processingDirectory"/concatenateBvec.bvec
        paste -d " " $(echo ${dwiFiles[@]/.nii.gz/.bval} | sort) > "$processingDirectory"/concatenateBval.bval
    else
        mrconvert ${dwiFiles[0]}        "$processingDirectory"/concatenateDwi.mif
        cp ${dwiFiles[0]/.nii.gz/.bvec} "$processingDirectory"/concatenateBvec.bvec
        cp ${dwiFiles[0]/.nii.gz/.bval} "$processingDirectory"/concatenateBval.bval
    fi

    # Get b-zero mean.
    dwiextract -force -nthreads $NSLOTS -fslgrad "$processingDirectory"/concatenateBvec.bvec "$processingDirectory"/concatenateBval.bval "$processingDirectory"/concatenateDwi.mif - -bzero | mrmath -force - mean "$processingDirectory"/meanMainPhase.mif -axis 3
    if [[ ! -z "$reversePhaseScan" ]]; then
        dwiextract -force -nthreads $NSLOTS -fslgrad ${reversePhaseScan/.nii.gz/.bvec} ${reversePhaseScan/.nii.gz/.bval} "$reversePhaseScan" - -bzero | mrmath - mean "$processingDirectory"/meanReversePhase.mif -axis 3
        mrcat -force -nthreads $NSLOTS "$processingDirectory"/meanMainPhase.mif "$processingDirectory"/meanReversePhase.mif "$processingDirectory"/reversePhasePair.nii.gz

        # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc).
        dimensions=$(fsl5.0-fslhd "$processingDirectory"/reversePhasePair.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}')
        newDimensions=$(echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}')
        do_cmd fsl5.0-fslroi "$processingDirectory"/reversePhasePair.nii.gz "$processingDirectory"/reversePhasePairCut.nii.gz $(echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /') 0 -1
        opt="-rpe_pair -se_epi $processingDirectory/reversePhasePairCut.nii.gz"
    else
        opt='-rpe_none'
    fi

    # Temporary fix for segmentation fault in topup.
    T=$(mktemp -d)
    cp /export01/local/fsl/bin/topup $T/fsl5.0-topup
    PATH=$T:$PATH

    # Preprocess each shell.
    for x in ${dwiFiles[@]}; do
        # Remove slices to make an even number of slices in all directions (requisite for dwi_preproc).
        dimensions=$(fsl5.0-fslhd "$x" | grep -E "^dim[1-3]" | awk '{print $NF}')
        newDimensions=$(echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}')
        do_cmd fsl5.0-fslroi "$x" "$processingDirectory"/tempDwiCut.nii.gz $(echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /') 0 -1

        tmp=$(basename $x)
        output="$processingDirectory"/${tmp/.nii.gz/preproc.mif}
        phaseEncodingDirection=$(readBids phaseEncodingDirection ${x/.nii.gz/.json})
        tempDirectory=$(mktemp -d -p ${processingDirectory})


        do_cmd dwipreproc $opt \
            -pe_dir "$phaseEncodingDirection" \
            -fslgrad ${x/.nii.gz/.bvec} ${x/.nii.gz/.bval} \
            -nocleanup \
            -tempdir "$tempDirectory" \
            -nthreads $NSLOTS \
            -force \
            "$processingDirectory"/tempDwiCut.nii.gz "$output"
        find "$tempDirectory" -type f -name "*mov*" -exec mv {} "$finalDirectory" \;
        nShells=$(echo $(cat "$finalDirectory"/dwi_post_eddy.eddy_movement_rms | wc -l) -1 | bc)
        mv "$finalDirectory"/dwi_post_eddy.eddy_movement_rms "$finalDirectory"/"$subject"_"$nShells"shells_dwi_post_eddy.eddy_movement_rms
        [[ -e "$finalDirectory"/field_movpar.txt ]] && mv "$finalDirectory"/field_movpar.txt "$finalDirectory"/"$subject"_"$nShells"shells_field_movpar.txt # Only exists if there is a fieldmap
        rm -rf "$tempDirectory"
    done
    rm $T/fsl5.0-topup
    find $T -empty -delete

    # Mask and denoise all shells, compute residuals for quality control
    if [[ $(find "$processingDirectory" -name "*preproc.mif" | wc -l) -gt 1 ]]; then
        do_cmd mrcat -force -nthreads $NSLOTS $(find "$processingDirectory" -name "*preproc.mif" | sort) "$processingDirectory"/concatenateDwiPreproc.mif
    else
        mv "$processingDirectory"/*preproc.mif "$processingDirectory"/concatenateDwiPreproc.mif
    fi
    do_cmd dwi2mask   -force -nthreads $NSLOTS "$processingDirectory"/concatenateDwiPreproc.mif "$finalDirectory"/"$subject"_dwiMask.mif
    do_cmd dwidenoise -force -nthreads $NSLOTS "$processingDirectory"/concatenateDwiPreproc.mif "$finalDirectory"/"$subject"_dwi.mif
    do_cmd mrcalc     -force -nthreads $NSLOTS "$processingDirectory"/concatenateDwiPreproc.mif "$finalDirectory"/"$subject"_dwi.mif -subtract "$finalDirectory"/"$subject"_dwi_residuals.mif
elif [[ $networkProcessing == true ]] ; then
    cp -rv $finalDirectoryReal/"$subject"_dwi.mif.gz $finalDirectoryReal/"$subject"_dwiMask.mif.gz $finalDirectory
fi

function uncompress() {
    [[ -f "$finalDirectory"/"$subject"_dwi.mif.gz ]] && do_cmd gunzip "$finalDirectory"/"$subject"_dwi.mif.gz
    [[ -f "$finalDirectory"/"$subject"_dwiMask.mif.gz ]] && do_cmd gunzip "$finalDirectory"/"$subject"_dwiMask.mif.gz
}

# Compute the alignment between structural and DWI
if [[ ! -f "$warpDirectoryReal"/"$subject"_fs2dti_bbr.lta ]] && [[ ! -f "$warpDirectory"/"$subject"_fs2dti_bbr.lta ]]; then
    uncompress
    do_cmd mrconvert "$finalDirectory"/"$subject"_dwi.mif "$processingDirectory"/"$subject"_dwi_forwarp.nii
    do_cmd bbregister --s "$fsDirectory" --mov "$processingDirectory"/"$subject"_dwi_forwarp.nii --reg "$warpDirectory"/"$subject"_dti2fs_bbr.lta --init-fsl --dti
    do_cmd lta_convert --inlta  "$warpDirectory"/"$subject"_dti2fs_bbr.lta --outlta  "$warpDirectory"/"$subject"_fs2dti_bbr.lta --invert
fi

# Get some basic metrics.
if [[ ! -f $finalDirectory/"$subject"_FA.mif.gz ]] || [[ ! -f $finalDirectory/"$subject"_ADC.mif.gz ]]; then
    if [[ ! -f $finalDirectoryReal/"$subject"_FA.mif.gz ]] || [[ ! -f $finalDirectoryReal/"$subject"_ADC.mif.gz ]]; then
        uncompress
        do_cmd dwi2tensor -nthreads $NSLOTS "$finalDirectory"/"$subject"_dwi.mif $processingDirectory/"$subject"_tensor.mif
        do_cmd tensor2metric -nthreads $NSLOTS -fa $finalDirectory/"$subject"_FA.mif.gz -adc $finalDirectory/"$subject"_ADC.mif.gz $processingDirectory/"$subject"_tensor.mif
    fi
fi

# Prepare tractography
if [[ ! -f $finalDirectory/"$subject"_FOD_WM_norm.mif.gz ]] && [[ ! -f $finalDirectory/"$subject"_FOD_norm.mif.gz ]]; then
    if [[ ! -f $finalDirectoryReal/"$subject"_FOD_WM_norm.mif.gz ]] && [[ ! -f $finalDirectoryReal/"$subject"_FOD_norm.mif.gz ]]; then
        uncompress
        if [[ $dwi2fodAlgorithm == msmt_csd ]]; then
            do_cmd dwi2response dhollander -nthreads $NSLOTS "$finalDirectory"/"$subject"_dwi.mif \
                "$finalDirectory"/"$subject"_dhollander_response_wm.txt \
                "$finalDirectory"/"$subject"_dhollander_response_gm.txt \
                "$finalDirectory"/"$subject"_dhollander_response_csf.txt \
                -mask "$finalDirectory"/"$subject"_dwiMask.mif
            do_cmd dwi2fod -nthreads $NSLOTS $dwi2fodAlgorithm "$finalDirectory"/"$subject"_dwi.mif \
                "$finalDirectory"/"$subject"_dhollander_response_wm.txt \
                "$processingDirectory"/wmFOD.mif \
                -mask "$finalDirectory"/"$subject"_dwiMask.mif
            do_cmd mtnormalise -nthreads $NSLOTS -mask "$finalDirectory"/"$subject"_dwiMask.mif \
                "$processingDirectory"/wmFOD.mif $finalDirectory/"$subject"_FOD_WM_norm.mif.gz
        else
            do_cmd dwi2response tournier -nthreads $NSLOTS \
                "$finalDirectory"/"$subject"_dwi.mif \
                "$finalDirectory"/"$subject"_tournier_response.txt \
                -mask "$finalDirectory"/"$subject"_dwiMask.mif
            do_cmd dwi2fod -nthreads $NSLOTS $dwi2fodAlgorithm \
                "$finalDirectory"/"$subject"_dwi.mif \
                "$finalDirectory"/"$subject"_tournier_response.txt  \
                "$processingDirectory"/FOD.mif \
                -mask "$finalDirectory"/"$subject"_dwiMask.mif
            do_cmd mtnormalise -nthreads $NSLOTS -mask "$finalDirectory"/"$subject"_dwiMask.mif \
                "$processingDirectory"/FOD.mif $finalDirectory/"$subject"_FOD_norm.mif.gz

        fi
    fi
fi

# Compress some files.
for x in $(find $finalDirectory -maxdepth 1 -mindepth 1 -name "*.mif"); do
    do_cmd gzip $x
done

if [[ $networkProcessing == true ]]; then
    mkdir -p $finalDirectoryReal
    mkdir -p $warpDirectoryReal
    do_cmd cp -rv $finalDirectory/* $finalDirectoryReal
    do_cmd cp -rv $warpDirectory/* $warpDirectoryReal
fi
