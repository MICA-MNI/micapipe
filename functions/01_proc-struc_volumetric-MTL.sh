#!/bin/bash
# Performs structural preprocessing, freesurfer preprocessing, and surfpatch segmentation.
# Will be the base for TEMPORAL proc
# Written by Reinder Vos de Wael (Oct 2018).

# source ~/.sge_profile
source $MICASOFT_DIR/pipelines/08_micaProcessing/mica_processingSupportFunctions.sh
source $MICASOFT_DIR/pipelines/09_bids_micaProcessing/utilities.sh

Title "Running structural processing: Volumetric-MTL"

while getopts ":np:t:" opt; do
    case $opt in
        n)
            networkProcessing=true ;;
        p)
            processingDirectory=$OPTARG ;;
        t)
            T1+=($OPTARG ) ;;
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

shift $((OPTIND-1))

# Set basic parameters.
mmTemplates=$1       # "0.8 2"
finalDirectory=$2    # $proc_struct
warpDirectory=$3     # $dir_warp
templateDirectory=$4 # $util_MNIvolumes
firstDirectory=$5    # $dir_first
subject=$6           # $id

# Set processing directory
if [[ ! -z $networkProcessing ]] && [[ -z $processingDirectory ]]; then
    echo "ERROR: If using network processing, processing directory must be specified."
    echo "Consider setting -p SETDEFAULT"
    exit
fi

if [[ $processingDirectory == SETDEFAULT ]]; then
    case $HOSTNAME in
        fladgate|yeatman)
            processingDirectory=$(mktemp -d -p /host/$HOSTNAME/local_raid/temporaryLocalProcessing) ;;
        cassio|varro)
            processingDirectory=$(mktemp -d -p /host/$HOSTNAME/export02/data/temporaryLocalProcessing) ;;
        *)
            processingDirectory=$(mktemp -d -p /data/mica2/temporaryNetworkProcessing/) ;;
    esac
elif [[ ! -z $processingDirectory ]]; then
    processingDirectory=$(mktemp -d -p $processingDirectory)
else
    processingDirectory=$(mktemp -d -p $finalDirectory)
fi

# Switch temporary folder depending on /tmp/ data capacity.
case $HOSTNAME in
    emilia|montano)
        tmp5ttgen=/data/mica2/temporaryNetworkProcessing ;;
    *)
        tmp5ttgen=$processingDirectory/ ;;
esac
# TRAP: Cleanup processing directory on exit.
function cleanup {
# Clear temporary files. Add some checks to be sure there's little chance of unintentional deletes.
if [[ $networkProcessing == true ]]; then
    find ${processingDirectory}/network/finalDirectory -maxdepth 1 -regextype posix-extended -regex ".*[.omat|.mif|.mif.gz|.nii|.nii.gz|.BRIK|.HEAD|.log|.mnc|.imp]" -delete
    find ${processingDirectory}/network/warpDirectory -maxdepth 1 -regextype posix-extended -regex ".*[.omat|.mif|.mif.gz|.nii|.nii.gz|.BRIK|.HEAD|.log|.mnc|.imp]" -delete
    find ${processingDirectory}/network/firstDirectory/${subject}_first_native.logs -mindepth 1 -maxdepth 1 -delete
    find ${processingDirectory}/network/firstDirectory -maxdepth 1 -regextype posix-extended -regex ".*[.omat|.mif|.mif.gz|.nii|.nii.gz|.BRIK|.HEAD|.log|.mnc|.imp|.bvars|.vtk|.com|.com2|.com3]" -delete
    find ${processingDirectory} -maxdepth 2 -regextype posix-extended -regex ".*[.omat|.mif|.mif.gz|.nii|.nii.gz|.BRIK|.HEAD|.log|.mnc|.imp]" -delete
    find ${processingDirectory} -empty -delete
else
    find ${processingDirectory} -type f -regextype posix-extended -regex ".*[.omat|.mif|.nii.gz|.BRIK|.HEAD|.log|.mnc|.imp]" -delete
    find ${processingDirectory} -mindepth 0 -maxdepth 0 -type d -empty -delete
fi
}
trap cleanup EXIT

if [[ $networkProcessing == true ]]; then
    finalDirectoryReal=$finalDirectory
    finalDirectory=$processingDirectory/network/finalDirectory
    warpDirectoryReal=$warpDirectory
    warpDirectory=$processingDirectory/network/warpDirectory
    firstDirectoryReal=$firstDirectory
    firstDirectory=$processingDirectory/network/firstDirectory

    mkdir -p $processingDirectory/network
    for i in $(seq 1 ${#T1[@]}); do
        cp ${T1[$i-1]} $processingDirectory/network/$(basename ${T1[$i-1]})
        T1[$i-1]=$processingDirectory/network/$(basename ${T1[$i-1]})
    done
fi
for x in ${finalDirectory} ${warpDirectory} "$firstDirectory"; do
    [[ ! -d "$x" ]] && mkdir -p "$x"
done
echo $processingDirectory
echo ${T1[@]}

# Reorient to FSL space
for rawNifti in ${T1[@]}; do
    voxelsize=`fslinfo "$rawNifti" | grep pixdim2 | awk '{printf("%1.1f\n",$2)}'` # Assumes isotropic voxels.
    voxelsize=${voxelsize%.0}
    Do_cmd 3dresample -orient RPI -prefix ${processingDirectory}/${subject}_t1w_"$RANDOM"_${voxelsize}mm_native.nii.gz -inset "$rawNifti"
done

####################################################################################################################################################################################
####################################################################################################################################################################################
# NOEL tools
####################################################################################################################################################################################
# NUC replaced by N4BiasFieldCorrection  then rescale intensity
# CLAMP wasn't changed, only deprecaited and was kept only for temporal lobe processing (SURFPATCH)
for x in `find ${processingDirectory} -name "*mm_native.nii.gz"`; do
    mica_niigz2mnc ${x}
    Do_cmd "$NOELSOFT_DIR"/SPIP/spip3T-nuc_and_clamp.sh ${x/.nii.gz/.mnc} -o ${processingDirectory}
done

# If multiple scans were provided, average all of them.
# This is replaced with FSL flirt
if [[ `find ${processingDirectory} -name "*clamp.mnc.gz" | wc -l` -gt 1 ]]; then
    Do_cmd "$NOELSOFT_DIR"/SPIP/spip-average.sh ${processingDirectory}/*clamp.mnc.gz ${processingDirectory}/${subject}_native_merge.mnc.gz
    Do_cmd "$NOELSOFT_DIR"/SPIP/spip3T-nuc_and_clamp.sh ${processingDirectory}/${subject}_native_merge.mnc.gz -o ${processingDirectory}
else
    # If only one scan, just rename it.
    Do_cmd find ${processingDirectory} -name "*clamp.mnc.gz" -exec mv {} ${processingDirectory}/${subject}_native_merge_clamp.mnc.gz \;
fi
Do_cmd mica_mncgz2niigz ${processingDirectory}/${subject}_native_merge_clamp.mnc.gz
rm -f ${processingDirectory}/${subject}_native_merge_clamp.mnc # FSL gets annoying when the .mnc file is still there.

#rm -f ${processingDirectory}/*.mnc* #Remove because FSL doesn't like having same-name .mnc and .nii files.
voxelsize=`fslinfo ${processingDirectory}/${subject}_native_merge_clamp.nii.gz | grep pixdim2 | awk '{printf("%1.1f\n",$2)}'` # Assumes isotropic voxels.
voxelsize=${voxelsize%.0}

####################################################################################################################################################################################
####################################################################################################################################################################################

# T1w native output names
T1str_nat=${finalDirectory}/${subject}_t1w_${voxelsize}mm_nativepro
T1nativepro=${finalDirectory}/${T1str_nat}.nii.gz
T1nativepro_brain=${T1nativepro/.nii.gz/_brain.nii.gz}
T1nativepro_first=${T1nativepro_brain/.nii.gz/_first.nii.gz}
T1nativepro_5tt=${T1nativepro/.nii.gz/_5TT.nii.gz}

# Noel processing to output directory
# NUC replaced by N4BiasFieldCorrection  then rescale intensity
# CLAMP wasn't changed, only deprecaited and was kept for temporal lobe processing
mv ${processingDirectory}/${subject}_native_merge_clamp.nii.gz $T1nativepro

# FSL tries to submit to SGE. We don't want that. ??
unset FSLPARALLEL

# Brainmask - CHANGE split this functions into 2 steps bet | MASK only with only ONE step
Do_cmd brainExtraction $T1nativepro $T1nativepro_brain ${templateDirectory} ${voxelsize}

# FSL first
Do_cmd run_first_all -b -i $T1nativepro_brain -o $T1nativepro_first

# FSL FAST DEACTIVATE THE NUC with -N
Do_cmd fast $T1nativepro_brain # add -N here

# Loop over all requested templates - could use a check for whether the template and the mask exists.
# mmTemplates is a fixed value=(0.8 2)
for mm in $mmTemplates; do
  # Only runs if the output doesn't exist
  if [[ ! -e ${finalDirectory}/${subject}_t1w_${mm}mm_nlMNI152_brain_pve_0.nii.gz ]]; then
      # MNI152 templates
      MNI152_brain=${templateDirectory}/MNI152_T1_${mm}mm_brain.nii.gz
      MNI152_mask=${templateDirectory}/MNI152_T1_${mm}mm_brain_mask.nii.gz

      # T1 on MNI152 templates
      T1_MNI152_l=${finalDirectory}/${subject}_t1w_${mm}mm_MNI152_l.nii.gz
      T1_MNI152_l_brain=${finalDirectory}/${subject}_t1w_${mm}mm_MNI152_l_brain.nii.gz
      T1_MNI152_nl=${finalDirectory}/${subject}_t1w_${mm}mm_MNI152_nl.nii.gz
      T1_MNI152_nl_brain=${finalDirectory}/${subject}_t1w_${mm}mm_MNI152_nl_brain.nii.gz

      # T1 to MNI152 transformation matrix and warp fields
      mat_MNI152_l=${warpDirectory}/${subject}_t1w_nativepro_brain_to_${mm}mm_MNI152_l_brain.omat
      T1_MNI152_l_warp=${warpDirectory}/${subject}_${mm}mm_MNI152_l_to_${mm}mm_MNI152_nl_warp.nii.gz
      T1_MNI152_nl_warp=${warpDirectory}/${subject}_t1w_nativepro_brain_to_${mm}mm_MNI152_nl_warp.nii.gz
      T1_MNI152_nl_Invwarp=${warpDirectory}/${subject}_${mm}mm_MNI152_nl_to_t1w_nativepro_brain_warp.nii.gz

      # Register to MNI space and store warps
      # Rigid registration: T1w native to MNI152
      Do_cmd flirt  -ref $MNI152_brain -in $T1nativepro_brain -out $T1_MNI152_l_brain -omat $mat_MNI152_l -cost mutualinfo -searchcost mutualinfo -dof 12 -interp trilinear
      # Deformable registration: T1w native to MNI152
      Do_cmd fnirt --ref=$MNI152_brain --in=$T1_MNI152_l_brain --fout=$T1_MNI152_l_warp --interp=linear --refmask=$MNI152_mask

      # Create warps
      Do_cmd convertwarp -m $mat_MNI152_l -w $T1_MNI152_l_warp -r $MNI152_brain -o $T1_MNI152_nl_warp
      Do_cmd invwarp -w $T1_MNI152_nl_warp -o $T1_MNI152_nl_Invwarp -r $T1nativepro_brain

      # Warp the T1 native to MNI152
      Do_cmd flirt -in $T1nativepro -out $T1_MNI152_l  -ref $MNI152_brain -applyxfm -init $mat_MNI152_l
      Do_cmd applywarp -i $T1nativepro -o $T1_MNI152_nl -r $MNI152_brain -w $T1_MNI152_nl_warp

      # Tissue classification
      # Tissue classification
      Do_cmd fslmaths $T1_MNI152_nl -mul $MNI152_mask $T1_MNI152_nl_brain
      # Why to apply FAST? why don't register nativepro-fast to MNI152_nl?
      Do_cmd fast -N -v $T1_MNI152_nl_brain
  fi
done

# Generate a five-tissue-type image for anatomically constrained tractography
if [[ ! -e $T1nativepro_5tt ]]; then
    Do_cmd 5ttgen fsl $T1nativepro_brain $T1nativepro_5tt -premasked -tempdir ${tmp5ttgen}
fi

if [[ $networkProcessing == true ]]; then
    mkdir -p $finalDirectoryReal
    mkdir -p $warpDirectoryReal
    mkdir -p $firstDirectoryReal
    Do_cmd cp -rv $finalDirectory/* $finalDirectoryReal
    Do_cmd cp -rv $warpDirectory/* $warpDirectoryReal
    Do_cmd cp -rv $firstDirectory/* $firstDirectoryReal
fi
