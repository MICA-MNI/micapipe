#!/bin/bash
# Resting state preprocessing
# Written by Casey Paquola and Reinder Vos De Wael (Oct 2018).
# and a tiny bit from Sara (Feb 2019)... 
# and a whole lot from Sara (August 2019) 
# ... YOLO

source ~/.sge_profile
source $MICASOFT_DIR/pipelines/08_micaProcessing/mica_processingSupportFunctions.sh
export MATLABPATH=$MATLABPATH:$MICASOFT_DIR/pipelines/08_micaProcessing/

echo "-------------------------------------------------------------"
echo "Warning: Topup correction can currently only use AP-PA scans!"
echo "-------------------------------------------------------------"
export OMP_NUM_THREADS=1 # Sets wb_command to only use one thread. 

# Deal with flags
while getopts ":a:b:hn:p:" opt; do
    case $opt in
        a) 
            mainPhaseScan=$OPTARG ;; 
        b)
            reversePhaseScan=$OPTARG ;;
        n)
            hostname=$OPTARG 
            networkProcessing=true ;;
        p)
            processingDirectoryArg=$OPTARG ;;
        h)
            usage 
            exit 1 ;;
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
subject=$1
mainScan=$2
outputDirectory=$3
warpOutputDirectory=$4
structuralDirectory=$5
struct2fs=$6
baseDirectory=${7}
rsTag=${8}
rsTagRegex=${9}
mainPhaseScan=${10}
reversePhaseScan=${11}
mnitemplate3mm=${12}
postrsfmri=${13}

if [[ -z $mainPhaseScan ]] && [[ ! -z $reversePhaseScan ]]; then
    echo "You must either specify two phase maps, or none at all."
    exit
fi
if [[ ! -z $networkProcessing ]] && [[ -z $processingDirectoryArg ]]; then
    echo "Must specify a processing directory if network processing is set to on."
    exit
fi

[[ -z $processingDirectoryArg ]] && processingDirectory=$(mktemp -d -p $outputDirectory) || processingDirectory=$(setProcessingDirectory $processingDirectoryArg)


surfaceTemplateDirectory=$MICASOFT_DIR/pipelines/08_micaProcessing/utilities/resample_fsaverage/

# define directories
freesurferDirectory="$structuralDirectory"/surfaces/
export SUBJECTS_DIR=$freesurferDirectory
conte69Directory="$freesurferDirectory"/conte69/
volumetricOutputDirectory="$outputDirectory"/volumetric/
surfaceOutputDirectory="$outputDirectory"/surfaces/
ICAOutputDirectory="$outputDirectory"/ICA_MELODIC/


if [[ $networkProcessing == true ]]; then
    mkdir $processingDirectory/network
    outputDirectoryReal=$outputDirectoryReal
    outputDirectoryReal=$processingDirectory/network/$(basename $outputDirectory)
    warpOutputDirectoryReal=$warpOutputDirectory
    warpOutputDirectory=$processingDirectory/network/$(basename $warpOutputDirectory)
    surfaceOutputDirectoryReal=$surfaceOutputDirectory
    volumetricOutputDirectoryReal=$volumetricOutputDirectory
    ICAOutputDirectoryReal=$ICAOutputDirectory
    
    [[ -f $struct2fs ]] && directAccess=true || directAccess=false
    if [[ $directAccess == true ]]; then
        do_cmd cp -v $struct2fs $processingDirectory/network/struct2fs.lta
        do_cmd cp -v $mainScan $processingDirectory/network/mainScan.nii.gz
        do_cmd cp -v ${mainScan/.nii.gz/}.json $processingDirectory/network/mainScan.json
    else
        do_cmd scp $USER@$hostname:$struct2fs $processingDirectory/network/struct2fs.lta
        do_cmd scp $USER@$hostname:$mainScan $processingDirectory/network/mainScan.nii.gz
        do_cmd scp $USER@$hostname:${mainScan/.nii.gz/}.json $processingDirectory/network/mainScan.json
    fi
    struct2fs=$processingDirectory/network/struct2fs.lta
    mainScan=$processingDirectory/network/mainScan.nii.gz
    if [[ $directAccess == true ]]; then
        do_cmd cp -vr \
            $structuralDirectory \
            $processingDirectory/network
    else
        do_cmd scp -r \
            $USER@$hostname:$structuralDirectory \
            $processingDirectory/network
    fi
    freesurferDirectory=$processingDirectory/network/$(basename $structuralDirectory)/$(basename $freesurferDirectory)
    conte69Directory=$freesurferDirectory/$(basename $conte69Directory)
    surfaceOutputDirectory=$processingDirectory/network/$(basename $surfaceOutputDirectory)   
     
    if [ "$mainPhaseScan" == "NOMAINPHASESCAN" ]; then
        mainPhaseScan=""
    else
        do_cmd scp $USER@$hostname:$mainPhaseScan $processingDirectory/network/mainPhaseScan.nii.gz
        mainPhaseScan=$processingDirectory/network/mainPhaseScan.nii.gz
    fi
    if [[ $reversePhaseScan == "NOREVERSESCAN" ]]; then
        reversePhaseScan=""
    else
        do_cmd scp $USER@$hostname:$reversePhaseScan $processingDirectory/network/reversePhaseScan.nii.gz
        reversePhaseScan=$processingDirectory/network/reversePhaseScan.nii.gz
    fi
fi

# Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script). 
for x in "$warpOutputDirectory" "$surfaceOutputDirectory" "$volumetricOutputDirectory" "$ICAOutputDirectory"; do
    [[ ! -d "$x" ]] && mkdir -p "$x"
done
[[ -z $processingDirectory ]] && processingDirectory=$(mktemp -d -p $outputDirectory) 

# gettin dat readout time
readoutTime=`cat ${mainScan/.nii.gz/}.json | grep "TotalReadoutTime" | grep -Eo [0-9].[0-9]+` 

toProcess="$reversePhaseScan $mainScan $mainPhaseScan"
# Loop over all scans for everything before motion correction across scans.
if [[ ! -f "$volumetricOutputDirectoryReal"/"$subject"_singleecho_fmrispace.nii.gz ]]; then
    step=-1
    tags=(mainScan mainPhaseScan reversePhaseScan)
    
    for rawNifti in ${toProcess}; do
        # Get basic parameters
        tag=${tags[$step]}
        step=$(echo $step + 1 | bc)
        
        echo "TAG"
        echo $tag
        
        echo "RAWNIFTI"
        echo $rawNifti
        
        # Drop first five TRs and reorient. 
        if [ "$tag" == "mainScan" ]; then
            do_cmd nifti_tool -cbl -prefix "$processingDirectory"/"$tag"_trDrop.nii.gz -infiles "$rawNifti"'[5..$]'
            do_cmd 3dresample -orient RPI -prefix "$processingDirectory"/"$tag"_reorient.nii.gz -inset "$processingDirectory"/"$tag"_trDrop.nii.gz
        else
            do_cmd 3dresample -orient RPI -prefix "$processingDirectory"/"$tag"_reorient.nii.gz -inset "$rawNifti"
        fi    
            
        # Remove slices to make an even number of slices in all directions (requisite for topup).
        dimensions=`fsl5.0-fslhd "$processingDirectory"/"$tag"_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
        newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
        do_cmd fsl5.0-fslroi "$processingDirectory"/"$tag"_reorient.nii.gz \
                             "$processingDirectory"/"$tag"_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` -1 -1
        
        # Motion correction within scans
        do_cmd fsl5.0-fslmaths "$processingDirectory"/"$tag"_sliceCut.nii.gz \
                                -Tmean "$processingDirectory"/"$tag"_sliceCutMean.nii.gz
        do_cmd 3dvolreg -Fourier -twopass -base "$processingDirectory"/"$tag"_sliceCutMean.nii.gz \
                        -zpad 4 -prefix "$processingDirectory"/"$tag"_mc.nii.gz \
                        -1Dfile "$volumetricOutputDirectory"/"$subject"_"$tag".1D \
                        "$processingDirectory"/"$tag"_sliceCut.nii.gz
        do_cmd fsl5.0-fslmaths "$processingDirectory"/"$tag"_mc.nii.gz \
                                -Tmean "$processingDirectory"/"$tag"_mcMean.nii.gz
    done
    
    do_cmd fsl_motion_outliers -i "$processingDirectory"/mainScan_sliceCut.nii.gz \
                               -o "$volumetricOutputDirectory"/"$subject"_singleecho_spikeRegressors_FD.1D \
                               -s "$volumetricOutputDirectory"/"$subject"_singleecho_metric_FD.1D --fd
    mv "$volumetricOutputDirectory"/"$subject"_mainScan.1D "$volumetricOutputDirectory"/"$subject"_singleecho.1D
    
    # Only do distortion correction if field maps were provided, if not then rename the scan to distortionCorrected (just to make the next lines of code easy). 
    if [[ $mainPhaseScan == "" ]] || [[ $reversePhaseScan == "" ]]; then         
        echo "++++++++++++++++++++"
        echo "skipped TOPUP!!!!!!!"
        echo "++++++++++++++++++++"
        mv -v "$processingDirectory"/mainScan_mc.nii.gz "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz
    else
        mainPhaseScanMean=`find "$processingDirectory"    -maxdepth 1 -name "*mainPhaseScan*_mcMean.nii.gz"`
        mainPhaseScan=`find "$processingDirectory"        -maxdepth 1 -name "*mainPhaseScan*_mc.nii.gz"`
        reversePhaseScanMean=`find "$processingDirectory" -maxdepth 1 -name "*reversePhaseScan*_mcMean.nii.gz"`
        reversePhaseScan=`find "$processingDirectory"     -maxdepth 1 -name "*reversePhaseScan*_mc.nii.gz"`
        mainScan=`find "$processingDirectory"             -maxdepth 1 -name "*mainScan*_mc.nii.gz"`
        
        do_cmd fsl5.0-flirt -in "$mainPhaseScanMean" -ref "$processingDirectory"/mainScan_mcMean.nii.gz -omat "$processingDirectory"/singleecho_tmpXfmMain.omat
        do_cmd fsl5.0-flirt -in "$reversePhaseScanMean" -ref "$processingDirectory"/mainScan_mcMean.nii.gz -omat "$processingDirectory"/singleecho_tmpXfmSecondary.omat
        
        do_cmd fsl5.0-flirt -in "$mainPhaseScan" -ref "$processingDirectory"/mainScan_mcMean.nii.gz -applyxfm -init "$processingDirectory"/singleecho_tmpXfmMain.omat -out "$processingDirectory"/singleecho_mainPhaseAligned.nii.gz
        do_cmd fsl5.0-flirt -in "$reversePhaseScan" -ref "$processingDirectory"/mainScan_mcMean.nii.gz -applyxfm -init "$processingDirectory"/singleecho_tmpXfmSecondary.omat -out "$processingDirectory"/singleecho_secondaryPhaseAligned.nii.gz
        
        do_cmd fsl5.0-fslmaths "$processingDirectory"/singleecho_mainPhaseAligned.nii.gz -Tmean "$processingDirectory"/singleecho_mainPhaseAlignedMean.nii.gz
        do_cmd fsl5.0-fslmaths "$processingDirectory"/singleecho_secondaryPhaseAligned.nii.gz -Tmean "$processingDirectory"/singleecho_secondaryPhaseAlignedMean.nii.gz
        
        # Distortion correction
        printf "0 1 0 $readoutTime \n0 -1 0 $readoutTime" > "$processingDirectory"/singleecho_topupDataIn.txt # Figure out how to set the numbers in this file correctly for any scan. Depends on phase encoding direction!
        do_cmd fsl5.0-fslmerge -t "$processingDirectory"/singleecho_mergeForTopUp.nii.gz "$processingDirectory"/singleecho_mainPhaseAlignedMean.nii.gz "$processingDirectory"/singleecho_secondaryPhaseAlignedMean.nii.gz
        do_cmd fsl5.0-topup --imain="$processingDirectory"/singleecho_mergeForTopUp.nii.gz --datain="$processingDirectory"/singleecho_topupDataIn.txt --config=b02b0.cnf --out="$processingDirectory"/singleecho_topup
        do_cmd fsl5.0-applytopup --imain="$mainScan" --inindex=1 --datain="$processingDirectory"/singleecho_topupDataIn.txt --topup="$processingDirectory"/singleecho_topup --method=jac --out="$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz
    fi
else
    echo "-------------------------------------"
    echo -e "\e[48;5;197m !!!!!  goin str8 to ICA-FIX yo  !!!!! \e[49m"
    echo "-------------------------------------"
    volumetricOutputDirectory="$volumetricOutputDirectoryReal"
fi

# creating mask for ICA-MELODIC then applying band-pass filtering
do_cmd fsl5.0-fslmaths "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz -Tmean "$volumetricOutputDirectory"/tmp.nii.gz
do_cmd fsl5.0-bet "$volumetricOutputDirectory"/tmp.nii.gz "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/func.nii.gz -m -n
do_cmd mv "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/func_mask.nii.gz "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/mask.nii.gz
rm -fv "$volumetricOutputDirectory"/tmp.nii.gz
rm -fv "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz
do_cmd 3dTproject -input "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz -prefix "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz -passband 0.01 666
do_cmd fsl5.0-fslmaths "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz -Tmean "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_mean_orig.nii.gz
do_cmd fsl5.0-fslmaths "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz -Tmean "$processingDirectory"/"$subject"_singleecho_fmrispace_mean.nii.gz

# run MELODIC for ICA-FIX
if [[ ! -f `find "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/filtered_func_data.ica/ -name "melodic_IC.nii.gz"` ]]; then

    do_cmd cp "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz "$ICAOutputDirectory"/filtered_func_data.nii.gz
    do_cmd fsl5.0-melodic --in="$ICAOutputDirectory"/filtered_func_data.nii.gz \
                                    --tr=0.6 \
                                    --nobet \
                                    --mask="$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/mask.nii.gz \
                                    --bgthreshold=3 \
                                    --mmthresh=0.5 \
                                    --report \
                                    --Oall \
                                    --outdir="$ICAOutputDirectory"/filtered_func_data.ica \
                                    --Omean="$ICAOutputDirectory"/mean_func.nii.gz
else
    echo "----------------------------------------"
    echo -e "\e[48;5;197m !!!!  DON'T PANIC; SKIPPIN MELODIC  !!!! \e[49m"
    echo "----------------------------------------"
    volumetricOutputDirectory="$volumetricOutputDirectoryReal"
fi

# Register to Freesurfer space
do_cmd bbregister --s "$subject" --mov "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_mean_orig.nii.gz --reg "$warpOutputDirectoryReal"/"$subject"_singleecho_fmri2fs.lta --o "$warpOutputDirectory"/"$subject"_singleecho_fmri2fs_outbbreg_FIX.nii.gz --init-fsl --bold 
do_cmd lta_convert --inlta "$warpOutputDirectoryReal"/"$subject"_singleecho_fmri2fs.lta --outlta "$warpOutputDirectoryReal"/"$subject"_singleecho_fs2fmri.lta --invert 

mkdir "$ICAOutputDirectory"/mc
mkdir "$ICAOutputDirectory"/reg

do_cmd cp "$volumetricOutputDirectory"/"$subject"_singleecho.1D "$ICAOutputDirectory"/mc/prefiltered_func_data_mcf.par
do_cmd fsl5.0-fslroi "$ICAOutputDirectory"/filtered_func_data.nii.gz "$ICAOutputDirectory"/reg/example_func.nii.gz 397 1
do_cmd cp "$structuralDirectory"/${subject}_t1w_0.8mm_native_brain.nii.gz "$ICAOutputDirectory"/reg/highres.nii.gz
do_cmd cp "$ICAOutputDirectory"/filtered_func_data.ica/mean.nii.gz "$ICAOutputDirectory"/mean_func.nii.gz

do_cmd mri_concatenate_lta "$warpOutputDirectoryReal"/${subject}_t1w_0.8mm_native_brain_t1w2fs.lta \
                           "$warpOutputDirectoryReal"/${subject}_singleecho_fs2fmri.lta \
                           "$warpOutputDirectoryReal"/${subject}_t1w_0.8mm_native_brain_t1w2fmri.lta 
                                            
do_cmd lta_convert --inlta "$warpOutputDirectoryReal"/${subject}_t1w_0.8mm_native_brain_t1w2fmri.lta --outfsl "$ICAOutputDirectory"/reg/highres2example_func.mat

do_cmd mri_vol2vol --mov "$structuralDirectory"/${subject}_t1w_0.8mm_native_brain.nii.gz \
                   --targ "$processingDirectory"/"$subject"_singleecho_fmrispace_mean.nii.gz \
                   --o "$ICAOutputDirectory"/filtered_func_data.ica/t1w2fmri_brain.nii.gz \
                   --lta "$warpOutputDirectoryReal"/${subject}_t1w_0.8mm_native_brain_t1w2fmri.lta
                   

# run ICA-FIX iff melodic has been run
if  [[ -f `find "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/filtered_func_data.ica/ -name "melodic_IC.nii.gz"` ]] && \
    [[ -f `find /host/fladgate/local_raid/MICA-MTL/ICAFIX_training/ -name "MICAMTL_training_15HC_15PX.RData"` ]]; then    
    # if an error occurs during ICA-FIX, try the command line below to install required R packages and re-run
    #Rscript -e "for (i in c('kernlab','ROCR','class','party','e1071','randomForest'))install.packages(i)"
    do_cmd /data_/mica1/01_programs/fix/fix "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/ \
					    /host/fladgate/local_raid/MICA-MTL/ICAFIX_training/MICAMTL_training_15HC_15PX.RData \
					    20 -m -h 100
else
    echo "---------------------------------------------------"
    echo -e "\e[48;5;197m !!!!  NO MELODIC OUTPUT OR NO ICA-FIX TRAINING FILE  !!!! \e[49m"
    echo "---------------------------------------------------"
fi
volumetricOutputDirectory="$volumetricOutputDirectoryReal"

# Change single echo files for clean ones
yes | do_cmd cp -rf "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/filtered_func_data_clean.nii.gz \
                "$processingDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz 

yes | do_cmd cp -rf "$baseDirectory"/proc_rsfmri_FIX/ICA_MELODIC/filtered_func_data_clean.nii.gz \
                "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz 

do_cmd fsl5.0-fslmaths "$processingDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz \
                       -Tmean "$processingDirectory"/"$subject"_singleecho_fmrispace_mean.nii.gz

# Calculate tissue-specific signals
do_cmd mri_concatenate_lta $struct2fs "$warpOutputDirectoryReal"/"$subject"_singleecho_fs2fmri.lta "$processingDirectory"/struct2fmri.lta
for idx in 0 1 2; do 
    [[ $idx == 0 ]] && tissue=CSF
    [[ $idx == 1 ]] && tissue=GM 
    [[ $idx == 2 ]] && tissue=WM
    tissuemap=$(find "$structuralDirectory" -name "*native_brain_pve_${idx}.nii.gz")
    do_cmd mri_vol2vol --mov  "$tissuemap" --targ  "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz --lta "$processingDirectory"/struct2fmri.lta --o "$processingDirectory"/"$subject"_singleecho_"$tissue".nii.gz
    do_cmd fslmeants -i "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz -o "$volumetricOutputDirectory"/"$subject"_singleecho_"$tissue".txt -m "$processingDirectory"/"$subject"_singleecho_"$tissue".nii.gz -w 
done

do_cmd fslmaths "$processingDirectory"/"$subject"_singleecho_WM.nii.gz -add  "$processingDirectory"/"$subject"_singleecho_GM.nii.gz -add  "$processingDirectory"/"$subject"_singleecho_CSF.nii.gz "$processingDirectory"/"$subject"_singleecho_WB.nii.gz
do_cmd fslmeants -i "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz -o "$volumetricOutputDirectory"/"$subject"_singleecho_global.txt -m "$processingDirectory"/"$subject"_singleecho_WB.nii.gz -w 

# Motion confound 
do_cmd fsl5.0-fsl_motion_outliers -i "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz -o "$volumetricOutputDirectory"/"$subject"_singleecho_spikeRegressors_REFRMS.1D -s "$volumetricOutputDirectory"/"$subject"_singleecho_metric_REFRMS.1D --refrms --nomoco 

# register to MNI space
:'  ------------> HEY RAUL!
not sure why I use the mean_orig.nii.gz file below... could be changed 
with "$subject"_singleecho_fmrispace_mean.nii.gz... I think
'
do_cmd fsl5.0-epi_reg --epi="$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_mean_orig.nii.gz \
                      --t1="$structuralDirectory"/${subject}_t1w_0.8mm_native.nii.gz \
                      --t1brain="$structuralDirectory"/${subject}_t1w_0.8mm_native_brain.nii.gz \
                      --out="$processingDirectory"/"$subject"_fmri2t1
                      
do_cmd fsl5.0-applywarp -i "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz \
                        -o "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_fmri2t12mni.nii.gz \
                        -w "$warpOutputDirectoryReal"/"$subject"_native2nlmni_0.8mm.nii.gz \
                        -r ${mnitemplate3mm} \
                        --premat="$processingDirectory"/"$subject"_fmri2t1.mat

# Register to surface
for x in lh rh; do 
    [[ $x == lh ]] && hemisphere=l || hemisphere=r
    
    # Map the non high-passed volumetric timeseries to the surface so we can compute tSNR
    do_cmd mri_vol2surf \
        --mov "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz \
        --reg "$warpOutputDirectoryReal"/"$subject"_singleecho_fmri2fs.lta \
        --projfrac-avg 0.2 0.8 0.1 \
        --trgsubject "$subject" \
        --interp trilinear \
        --hemi "$x" \
        --out "$surfaceOutputDirectory"/"$subject"_singleecho_fmri2fs_"$x"_NoHP.mgh
    
    # Map high-passed timeseries to surface - this is what will be used to generate the connectomes later
    do_cmd mri_vol2surf \
        --mov "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_HP.nii.gz \
        --reg "$warpOutputDirectoryReal"/"$subject"_singleecho_fmri2fs.lta \
        --projfrac-avg 0.2 0.8 0.1 \
        --trgsubject "$subject" \
        --interp trilinear \
        --hemi "$x" \
        --out "$surfaceOutputDirectory"/"$subject"_singleecho_fmri2fs_"$x".mgh
    do_cmd mri_convert "$surfaceOutputDirectory"/"$subject"_singleecho_fmri2fs_"$x".mgh "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x".func.gii
    
    # Register to conte69 and apply smooth
    do_cmd wb_command -metric-resample \
        "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x".func.gii \
        "$conte69Directory"/"$subject"_"$hemisphere"h_sphereReg.surf.gii \
        "$surfaceTemplateDirectory"/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
        ADAP_BARY_AREA \
        "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k.func.gii \
        -area-surfs \
        "$freesurferDirectory"/"$subject"/surf/"$hemisphere"h.midthickness.surf.gii \
        "$conte69Directory"/"$subject"_"$hemisphere"h_midthickness_32k_fs_LR.surf.gii
    do_cmd wb_command -metric-smoothing \
        "$surfaceTemplateDirectory"/fsaverage.${hemisphere^^}.midthickness_orig.32k_fs_LR.surf.gii \
        "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k.func.gii \
        10 \
        "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k_10mm.func.gii
    do_cmd mri_convert "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k.func.gii "$surfaceOutputDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k.mgh
    do_cmd mri_convert "$processingDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k_10mm.func.gii "$surfaceOutputDirectory"/"$subject"_singleecho_fmri2fs_"$x"_c69-32k_10mm.mgh
    
done

:'  +++++++++++++++++ S U B C O R T E X ++++++++++++++++++
---------> CHANGE ME HERE RAUL <--------------'
# Flirt fsfast segmentation (nativepro) to rsfmri space using the inverse of .mat on line 323 (epi_reg)
do_cmd fsl5.0-convert_xfm -omat "$processingDirectory"/"$subject"_t12fmri.mat \
                          -inverse "$processingDirectory"/"$subject"_fmri2t1.mat

do_cmd fsl5.0-flirt -in HC012_t1w_0.8mm_nativepro_all_fast_firstseg.nii.gz \                                   # FSFAST OUTPUT EXAMPLE
                    -ref "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_mean_orig.nii.gz \       # MEAN RSFMRI FROM LINE 320 
                    -applyxfm -init "$processingDirectory"/"$subject"_t12fmri.mat \                            # FROM LINE 378
                    -out HC012_t1w_0.8mm_rsfmri_all_fast_firstseg.nii.gz                                       # OUTPUTS FIRST SEG IN RSFMRI SPACE

# Extract subcortical timeseries (mean, one ts per first-segemented structure, exluding the brainstem (ew))
do_cmd mri_segstats --i "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz \                  # RSFMRI IN ITS OWN SPACE
                    --seg HC012_t1w_0.8mm_rsfmri_all_fast_firstseg.nii.gz \                                    # OUTPUT FROM COMMAND JUST BEFORE
                    --exclude 0 --exclude 16 \                                                                 # EXLUDE EVERYTHING BUT ACC, AMY, CAUD, HIP, PAL, PUT, THAL
                    --avgwf "$volumetricOutputDirectory"/"$subject"_FIRST_rsfmri-timeseries.txt                # ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
:' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


:'  +++++++++++++++++ C E R E B E L L U M +++++++++++++++++
---------> CHANGE ME HERE RAUL <--------------'
do_cmd fsl5.0-flirt -in <cereballar parcellation in nativepro> \                                                    # FSFAST OUTPUT EXAMPLE
                    -ref "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace_mean_orig.nii.gz \            # MEAN RSFMRI FROM LINE 320 
                    -applyxfm -init "$processingDirectory"/"$subject"_t12fmri.mat \                                 # FROM LINE 378
                    -out HC012_t1w_0.8mm_rsfmri_cerebellar.nii.gz                                                   # OUTPUTS FIRST SEG IN RSFMRI SPACE

# Extract subcortical timeseries (mean, one ts per first-segemented structure, exluding the brainstem (ew))
do_cmd mri_segstats --i "$volumetricOutputDirectory"/"$subject"_singleecho_fmrispace.nii.gz \                       # RSFMRI IN ITS OWN SPACE
                    --seg HC012_t1w_0.8mm_rsfmri_cerebellar.nii.gz \                                                # OUTPUT FROM COMMAND JUST BEFORE
                    --exclude 0                                                                                     # EXLUDE EVERYTHING BUT CEREBELLUM
                    --avgwf "$volumetricOutputDirectory"/"$subject"_cerebellar_rsfmri-timeseries.txt                # ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
:' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

# copy stuff back
if [[ $networkProcessing == true ]]; then
    if [[ $directAccess == true ]]; then
        mkdir -p $volumetricOutputDirectoryReal
        mkdir -p $warpOutputDirectoryReal
        mkdir -p $surfaceOutputDirectoryReal
        mkdir -p $ICAOutputDirectoryReal
        yes | do_cmd cp -rv $volumetricOutputDirectory/* $volumetricOutputDirectoryReal
        yes | do_cmd cp -rv $warpOutputDirectory/* $warpOutputDirectoryReal
        yes | do_cmd cp -rv $surfaceOutputDirectory/* $surfaceOutputDirectoryReal
        yes | do_cmd cp -rv $ICAOutputDirectory/* $ICAOutputDirectoryReal
    else
        do_cmd scp -r $volumetricDirectory $USER@$hostname:$warpOutputDirectoryReal
        do_cmd scp -r $surfaceOutputDirectory $USER@$hostname:$surfaceOutputDirectoryReal
        do_cmd scp -r $warpOutputDirectory/* $USER@$hostname:$warpOutputDirectoryReal
        do_cmd scp -r $ICAOutputDirectory/* $USER@$hostname:$ICAOutputDirectoryReal
    fi
fi

# Clear temporary files. Add some checks to be sure there's little chance of unintentional deletes. 
find "$processingDirectory"/network -delete
find "$processingDirectory" -mindepth 1 -maxdepth 1 -type f -regextype posix-extended -regex ".*[.txt|.nii.gz|.omat|_log|.1D|.gii]" -delete 

# Remove processing directory
find "$processingDirectory" -mindepth 0 -maxdepth 0 -type d -empty -delete

# run post-rsfmri | OLD VERSION USING MATLAB
#if [[ $postrsfmri == "MATLAB" ]]; then
		#$MICASOFT_DIR/pipelines/08_micaProcessing/mica_03b_postrsfmriSingleEcho.sh
	#elif [[ $postrsfmri == "python" ]]; then
		#python $MICASOFT_DIR/pipelines/08_micaProcessing/mica_03b_postrsfmriSingleEcho.py "'$outputDirectory'"
	#else 
		#echo "---------------------------------------------------"
		#echo -e "\e[48;5;197m !!!!  NO POST-RSFMRI METHOD SPECIFIED  !!!! \e[49m"
		#echo "---------------------------------------------------"
#done

# run post-rsfmri
labelDirectory="$freesurferDirectory"/"$subject"/label/
python $MICASOFT_DIR/pipelines/08_micaProcessing/mica_03b_postrsfmriSingleEcho.py "$subject" "$outputDirectory" "$labelDirectory"

