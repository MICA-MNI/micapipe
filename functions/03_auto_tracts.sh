#!/bin/bash
#
export FSLOUTPUTTYPE=NIFTI_GZ
dir_functions=`dirname $(realpath $0)`
MICAPIPE=`dirname $(realpath $dir_functions)`
source ${MICAPIPE}/functions/utilities.sh

help() {
echo -e "
COMMAND:
    `basename $0`

        Perform automatic virtual dissection of a full-brain tractogram.
        Provide a pre-computed tck file and it will be dissected.

        This is an adaptation of AutoPtx and XTRACT to work with MRtrix3.
        See these links:
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/AutoPtx
        https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/XTRACT
        - The adaptation is only partial, as in AutoPtx and XTRACT the seeding of Streamlines
          is performed for each bundle. Here, a full tractogram is provided, where the user
          had the option to seed with whatever strategy was preferred. It is assumed, however,
          that a full-brain seeding approach was used (white matter mask, GM/WM border, etc.)
        - Another difference is the use of STOP ROIs, which only make sense if seeding
          is performed per bundle, and not if filtering a full tractogram.
          STOP ROIs are used, nonetheless, but rather as termination criteria that will truncate
          the streamlines. Thus, STOP ROIs should be much larger than usual, to avoid
          the appearance of multiple short-length truncated streamlines.
        - The quality of the full brain tractogram will determine the quality of bundle separation.
          It is highly recommended to provide a tractogram with more than one million streamlines,
          and one that has been checked for errors. Strategies such as anatomically-contstrained
          tractography (ACT) and spherical deconvolution informed filtering of tractograms (SIFT),
          both available in MRTrix3 should aid in obtaining such high-quality tractograms.

          Requires ANTs v2.3.3   (https://github.com/ANTsX/ANTs)


ARGUMENTS:
    Compulsory:
      -tck <file>       : .tck file (SIFTED, preferraly.)
      -outbase <string> : base name for all your outputs.
      -mask <file>      : Binary mask in subject dwi space.
      -fa <file>        : FA map in subject dwi space. Used for registration to template.


OPTIONS:
      -h|-help
      -keep_tmp         : Do not delete temporary directory
      -tmpDir <path>    : Specify location of temporary directory.
      -minStreamlinesPerVoxel <int>  Streamlines are truncated if voxel contains
                                     less than this number of streamlines.
                                     Default is $minStreamlinesPerVoxel
      -robust           : This option to runs a more ROBUST SyN registration ( More computation time )
      -weights <file>   : Add this option if you calculated obtained a weights file from SIFT2
      -threads <int>    : Number of threads (Default is 6)



USAGE:
    \033[38;5;141m`basename $0`\033[0m  \033[38;5;197m-tck\033[0m <file> \033[38;5;197m-outbase\033[0m <string> \033[38;5;197m-mask\033[0m <file> \033[38;5;197m-fa\033[0m <file>\n


Created by Luis Concha, INB, UNAM, lconcha@unam.mx (https://github.com/lconcha/auto_tracto)
Modified by rcruces, McGill University, MNI, MICA-lab Nov 2020
"
}

if [ "$#" -lt 2 ]; then Error "Not enough arguments\033[0m"; help; exit 2; fi

for arg in "$@"
do
  case "$arg" in
  -h|-help)
    help
    exit 1
  ;;
  -tck)
    tckIN=$2
    shift;shift
  ;;
  -outbase)
    outbase=$2
    shift;shift
  ;;
  -mask)
    mask=$2
    shift;shift
  ;;
  -fa)
    fa=$2
    shift;shift
  ;;
  -weights)
    tck_weights=$2
    shift;shift
  ;;
  -tmpDir)
    tmp=$2
    shift;shift;
  ;;
  -keep_tmp)
    keep_tmp=1
    shift
  ;;
  -robust)
    robust=TRUE
    shift
  ;;
  -threads)
    threads=$2
    shift;shift
  ;;
  -minStreamlinesPerVoxel)
    minStreamlinesPerVoxel=$2
    shift;shift
  ;;
  -*)
    echo "Unknown option ${2}"
    help
    exit 1
  ;;
  esac
done

# -----------------------------------------------------------------------------------------------
## Argument checks
arg=($tckIN $mask $outbase $fa)
if [ "${#arg[@]}" -lt 4 ]; then
  Error "One or more mandatory arguments are missing:"
  Note "-tck     :" $tckIN
  Note "-mask    :" "$mask"
  Note "-outbase :" "$outbase"
  Note "-fa      :" "$fa"
help; exit 1; fi

if [ ! -f $tckIN ]; then Error "File not found -tck:\n\t\t $tckIN"; help; exit 2; fi
if [ ! -f $mask ]; then Error "File not found -mask:\n\t\t $mask"; help; exit 2; fi
if [ ! -f $fa ]; then Error "File not found  -fa:\n\t\t $fa"; help; exit 2; fi
if [[ $fa != *".nii"* ]]; then Error "-fa is not NIFTI or NIFTI_GZ"; help; exit 2; fi

# Inputs and variables
autoPtx=$MICAPIPE/MNI152Volumes/protocols
atlas=$MICAPIPE/MNI152Volumes/FMRIB58_FA_1mm.nii.gz
Info "Checking inputs and variables"
Note "-tck      :" $tckIN
Note "-outbase  :" $outbase
Note "-mask     :" $mask
Note "-fa       :" $fa
Note "Protocols :" $fa
Note "FA atlas  :" $fa

tckIN=`realpath $tckIN`
mask=`realpath $mask`
outbase=`realpath $outbase`_
fa=`realpath $fa`
here=`pwd`

if [ -z ${minStreamlinesPerVoxel} ]; then minStreamlinesPerVoxel=1; fi
if [ -z ${tck_weights} ]; then tck_weights=""; else tck_weights="-tck_weights_in $tck_weights"; fi
structures=`ls -1 $autoPtx`

#------------------------------------------------------------------------------#
Title "Running Auto-Tractography segmentation"

#	Timer
aloita=$(date +%s)

# Create temp directory
if [ -z ${tmp} ]; then tmp=tmp_autotract_$$; else tmp=$tmp/autotract_$$; fi
Do_cmd mkdir -p $tmp

# TRAP in case the script fails
trap cleanup INT TERM

cd $tmp

# -----------------------------------------------------------------------------------------------
## Registration of atlas to subject DWI space
str_fa2atlas=$tmp/fa2atlas_
mat_fa2atlas=${str_fa2atlas}0GenericAffine.mat
mat_fa2atlas_warp=${str_fa2atlas}1Warp.nii.gz
mat_fa2atlas_Invwarp=${str_fa2atlas}1InverseWarp.nii.gz

Info "Calculating transformations: FA to FMRIB58_FA_1mm"
if [[ ${robust} == TRUE ]]; then
    Do_cmd antsRegistrationSyN.sh -d 3 -f $atlas -m $fa -o $str_fa2atlas -t s -n 6
else
    Do_cmd antsRegistrationSyNQuick.sh -d 3 -f $atlas -m $fa -o $str_fa2atlas -t s -n 6
fi

Do_cmd antsApplyTransforms -d 3 -e 3 -i $atlas -r $fa -n NearestNeighbor -t [$mat_fa2atlas,1] -t $mat_fa2atlas_Invwarp -o ${outbase}atlas2fa.nii.gz -v -u int

# -----------------------------------------------------------------------------------------------
filter(){
    if [ -f ${autoPtx}/${st}/skip ]
    then
      Info "Skipping structure $st because folder it has a skip file: ${autoPtx}/${st}/skip"
      return 0
    fi
    st=$1
    seed=${autoPtx}/${st}/seed.nii.gz
    target=${autoPtx}/${st}/target.nii.gz
    target2=${autoPtx}/${st}/target_02.nii.gz
    exclude=${autoPtx}/${st}/exclude.nii.gz
    stop=${autoPtx}/${st}/stop.nii.gz
    nat_seed=${tmp}/${st}_nat_seed.nii.gz
    nat_target=${tmp}/${st}_nat_target.nii.gz
    nat_exclude=${tmp}/${st}_nat_exclude.nii.gz
    nat_stop=${tmp}/${st}_nat_stop.nii.gz
    nat_mask=${tmp}/${st}_nat_mask.nii.gz
    summary=${outbase}summary.txt

    inc_seed=${tmp}/${st}_seed.nii.gz
    inc_target=${tmp}/${st}_target.nii.gz
    inc_target2=${tmp}/${st}_target_02.nii.gz
    exclude_interp=${tmp}/${st}_exclude_interp.nii.gz

    # Apply transformations
    Do_cmd antsApplyTransforms -r $fa -i $target -d 3 -e 3 -n GenericLabel -t [$mat_fa2atlas,1] -t $mat_fa2atlas_Invwarp -o $inc_target -v -u int
    Do_cmd antsApplyTransforms -r $fa -i $seed -d 3 -e 3 -n GenericLabel -t [$mat_fa2atlas,1] -t $mat_fa2atlas_Invwarp -o $inc_seed -v -u int
    Do_cmd antsApplyTransforms -r $fa -i $exclude -d 3 -e 3 -n GenericLabel -t [$mat_fa2atlas,1] -t $mat_fa2atlas_Invwarp -o $exclude_interp -v -u int
    if [ -f $target2 ]; then
        Do_cmd antsApplyTransforms -r $fa -i $target2 -d 3 -e 3 -n GenericLabel -t [$mat_fa2atlas,1] -t $mat_fa2atlas_Invwarp -o $inc_target2 -v -u int
        include="$inc_target -include $inc_target2"
    else
        include="$inc_target"
    fi
    mrcalc -quiet $exclude_interp 0 -gt - | maskfilter -force -quiet - dilate $nat_exclude

    # -----------------------------------------------------------------------------------------------
    # modify the mask according to the stop criterion
    if [ `imtest $stop` -eq 1 ]; then
         Info "Found $stop"
         Do_cmd antsApplyTransforms -r $fa -i $stop -d 3 -e 3 -n GenericLabel -t [$mat_fa2atlas,1] -t $mat_fa2atlas_Invwarp -o ${tmp}/${st}_stop_interp.nii.gz -v -u int
         mrcalc -quiet ${tmp}/${st}_stop_interp.nii.gz 0 -gt - | maskfilter -force -quiet - dilate $nat_stop
         Do_cmd mrcalc -force -quiet $fa 0 -gt ${tmp}/mask.nii.gz
         Do_cmd mrcalc -force -quiet ${tmp}/mask.nii.gz $nat_stop -subtract $nat_mask
     else
          Do_cmd mrcalc -force $fa 0 -gt $nat_mask
     fi

     Do_cmd tckedit -force $tck_weights\
                $tckIN \
                ${outbase}${st}.tck \
                -include $inc_seed\
                -include $include\
                -mask $nat_mask \
                -exclude $nat_exclude

     if [ $minStreamlinesPerVoxel -gt 1 ]; then
          Info "Truncating streamlines if streamlines in voxel is less than $minStreamlinesPerVoxel"
          Do_cmd tckmap -force -quiet -template $fa ${outbase}${st}.tck ${tmp}/${st}_n.nii
          Do_cmd mrcalc ${tmp}/${st}_n.nii.gz $minStreamlinesPerVoxel -ge ${tmp}/${st}_streamlinesMask.nii
          Do_cmd tckedit -force $tck_weights\
                    -mask ${tmp}/${st}_streamlinesMask.nii.gz \
                    ${outbase}${st}.tck \
                    ${outbase}${st}_masked.tck
     fi

     nTcks=`tckinfo -count ${outbase}${st}.tck | grep "actual count" | awk '{print $NF}'`
     echo "$st $nTcks" >> $summary
}
# End of filter

# -----------------------------------------------------------------------------------------------
## Main loop
structures=`ls -1 $autoPtx`
for st in $structures; do
  Info "Working on $st"; filter $st
done

# -----------------------------------------------------------------------------------------------
## Clean up
if [ $keep_tmp -eq 0 ]; then
    Info "Auto-tract: Deleting tmpDir $tmp"
    Do_cmd rm -fR $tmp
else
    Info "Auto-tract: Not deleting tmpDir $tmp"
fi

# -----------------------------------------------------------------------------------------------
Do_cmd cd $here

# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "Auto-Tracto ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m"
bids_variables_unset
