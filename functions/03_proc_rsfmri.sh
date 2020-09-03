#!/bin/bash
# Resting state preprocessing
# Written by Casey Paquola and Reinder Vos De Wael (Oct 2018).
# and a tiny bit from Sara (Feb 2019)...
# and a whole lot from Sara (August 2019)
# and incorporation to mica-pipe by Raul (August-September 2020)
# ... YOLO
#
# Resting state fMRI processing with bash:
#
# Preprocessing workflow for rsfmri.
#
# This workflow makes use of AFNI, FSL, ANTs, FIX
#
# Atlas an templates are avaliable from:
#
# https://github.com/MICA-MNI/micaopen/templates
#
#   ARGUMENTS order:
#   $1 : BIDS directory
#   $2 : participant
#   $3 : Out Directory
#   $4 : Temporal directory (default /tmp)
#
# ONLY for scripting and debugging:
# TEST=ON
# source utilities
source $MICAPIPE/functions/utilities.sh

BIDS=$1
id=$2
out=$3
tmp=$4

# Assigns variables names
bids_variables $BIDS $id $out
# print the names on the terminal
bids_print.variables-rsfmri

#------------------------------------------------------------------------------#
Title "Running MICA rsfMRI processing"

# GLOBAL variables for this script
Info "wb_command will use $OMP_NUM_THREADS threads"; fi
Warning "Warning: Topup correction can currently only use AP-PA scans!"

#	Timer
aloita=$(date +%s)
here=`pwd`

# Check tmp dir: temporary directory << processingDirectory
if [ -z ${tmp} ]; then tmp=/tmp; fi
tmp=${tmp}/${RANDOM}_proc-rsfmri_${id}
if [ ! -d $tmp ]; then Do_cmd mkdir -p $tmp; fi

# Set basic parameters.
struct2fs=$(find $dir_warp -name "*t1w2fs.lta")
rsTag="*rsfmri*3mm*bold*AP"
rsTagRegex=".*rsfmri.*3mm.*bold.*AP.*"
mnitemplate3mm=${util_MNIvolumes}/MNI152_T1_3mm.nii.gz

# Check inputs: rsfMRI and phase encoding
if [ ! -f ${bids_fmri} ]; then Error "Subject $id doesn't have acq-AP_bold: \n\t ${subject_bids}/func/"; exit; fi
if [ ! -f ${mainPhaseScan} ]; then Warning "Subject $id doesn't have acq-APse_bold: TOPUP will be skipped"; fi
if [ ! -f ${reversePhaseScan} ]; then Warning "Subject $id doesn't have acq-PAse_bold: TOPUP will be skipped"; fi
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro: run -proc_volumetric"; exit; fi

# Define directories
export SUBJECTS_DIR=$dir_surf

# rsfMRI directories
rsfmri_volum=${proc_rsfmri}/volumetric/   # volumetricOutputDirectory
rsfmri_surf=${proc_rsfmri}/surfaces/      # surfaceOutputDirectory
rsfmri_ICA=$proc_rsfmri/ICA_MELODIC       # ICAOutputDirectory

# Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script).
for x in $rsfmri_surf $rsfmri_volum $rsfmri_ICA; do
    [[ ! -d ${x} ]] && mkdir -p ${x}
done

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# bbreg missing for struct2fs

#------------------------------------------------------------------------------#
# Begining of the REAL processing
# gettin dat readout time
readoutTime=`cat ${mainScan/.nii.gz/}.json | grep "TotalReadoutTime" | grep -Eo [0-9].[0-9]+`

toProcess="$reversePhaseScan $mainScan $mainPhaseScan"

# Loop over all scans for everything before motion correction across scans.
if [[ ! -f ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz ]]; then
    step=-1
    tags=(mainScan mainPhaseScan reversePhaseScan) # this does not match with $toProcess

    for rawNifti in ${toProcess}; do
        # Get basic parameters
        tag=${tags[$step]}
        step=$(echo $step + 1 | bc)

        Info  "TAG: $tag \n\tRAWNIFTI: $rawNifti"

        # Drop first five TRs and reorient.
        if [ "$tag" == "mainScan" ]; then
            Do_cmd nifti_tool -cbl -prefix ${tmp}/${tag}_trDrop.nii.gz -infiles "$rawNifti"'[5..$]'
            Do_cmd 3dresample -orient RPI -prefix ${tmp}/${tag}_reorient.nii.gz -inset ${tmp}/${tag}_trDrop.nii.gz
        else
            Do_cmd 3dresample -orient RPI -prefix ${tmp}/${tag}_reorient.nii.gz -inset $rawNifti
        fi

        # Remove slices to make an even number of slices in all directions (requisite for topup).
        dimensions=`fslhd ${tmp}/${tag}_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
        newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
        Do_cmd fslroi ${tmp}/${tag}_reorient.nii.gz \
                             ${tmp}/${tag}_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` -1 -1

        # Motion correction within scans
        Do_cmd fslmaths ${tmp}/${tag}_sliceCut.nii.gz \
                                -Tmean ${tmp}/${tag}_sliceCutMean.nii.gz
        Do_cmd 3dvolreg -Fourier -twopass -base ${tmp}/${tag}_sliceCutMean.nii.gz \
                        -zpad 4 -prefix ${tmp}/${tag}_mc.nii.gz \
                        -1Dfile ${rsfmri_volum}/${id}_${tag}.1D \
                        ${tmp}/${tag}_sliceCut.nii.gz
        Do_cmd fslmaths ${tmp}/${tag}_mc.nii.gz \
                                -Tmean ${tmp}/${tag}_mcMean.nii.gz
    done

    Do_cmd fsl_motion_outliers -i ${tmp}/mainScan_sliceCut.nii.gz \
                               -o ${rsfmri_volum}/${id}_singleecho_spikeRegressors_FD.1D \
                               -s ${rsfmri_volum}/${id}_singleecho_metric_FD.1D --fd
    mv ${rsfmri_volum}/${id}_mainScan.1D ${rsfmri_volum}/${id}_singleecho.1D

    # Only do distortion correction if field maps were provided, if not then rename the scan to distortionCorrected (just to make the next lines of code easy).
    if [[ $mainPhaseScan == "" ]] || [[ $reversePhaseScan == "" ]]; then
        Warning "skipped TOPUP!!!!!!!"
        mv -v ${tmp}/mainScan_mc.nii.gz ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz
    else
        mainPhaseScanMean=`find ${tmp}    -maxdepth 1 -name "*mainPhaseScan*_mcMean.nii.gz"`
        mainPhaseScan=`find ${tmp}        -maxdepth 1 -name "*mainPhaseScan*_mc.nii.gz"`
        reversePhaseScanMean=`find ${tmp} -maxdepth 1 -name "*reversePhaseScan*_mcMean.nii.gz"`
        reversePhaseScan=`find ${tmp}     -maxdepth 1 -name "*reversePhaseScan*_mc.nii.gz"`
        mainScan=`find ${tmp}             -maxdepth 1 -name "*mainScan*_mc.nii.gz"`

        Do_cmd flirt -in $mainPhaseScanMean -ref ${tmp}/mainScan_mcMean.nii.gz -omat ${tmp}/singleecho_tmpXfmMain.omat
        Do_cmd flirt -in $reversePhaseScanMean -ref ${tmp}/mainScan_mcMean.nii.gz -omat ${tmp}/singleecho_tmpXfmSecondary.omat

        Do_cmd flirt -in $mainPhaseScan -ref ${tmp}/mainScan_mcMean.nii.gz -applyxfm -init ${tmp}/singleecho_tmpXfmMain.omat -out ${tmp}/singleecho_mainPhaseAligned.nii.gz
        Do_cmd flirt -in $reversePhaseScan -ref ${tmp}/mainScan_mcMean.nii.gz -applyxfm -init ${tmp}/singleecho_tmpXfmSecondary.omat -out ${tmp}/singleecho_secondaryPhaseAligned.nii.gz

        Do_cmd fslmaths ${tmp}/singleecho_mainPhaseAligned.nii.gz -Tmean ${tmp}/singleecho_mainPhaseAlignedMean.nii.gz
        Do_cmd fslmaths ${tmp}/singleecho_secondaryPhaseAligned.nii.gz -Tmean ${tmp}/singleecho_secondaryPhaseAlignedMean.nii.gz

        # Distortion correction
        printf "0 1 0 $readoutTime \n0 -1 0 $readoutTime" > ${tmp}/singleecho_topupDataIn.txt # Figure out how to set the numbers in this file correctly for any scan. Depends on phase encoding direction!
        Do_cmd fslmerge -t ${tmp}/singleecho_mergeForTopUp.nii.gz ${tmp}/singleecho_mainPhaseAlignedMean.nii.gz ${tmp}/singleecho_secondaryPhaseAlignedMean.nii.gz
        Do_cmd topup --imain=${tmp}/singleecho_mergeForTopUp.nii.gz --datain=${tmp}/singleecho_topupDataIn.txt --config=b02b0.cnf --out=${tmp}/singleecho_topup
        Do_cmd applytopup --imain="$mainScan" --inindex=1 --datain=${tmp}/singleecho_topupDataIn.txt --topup=${tmp}/singleecho_topup --method=jac --out=${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz
    fi
else
    Warning "!!!!!  goin str8 to ICA-FIX yo  !!!!!"
fi

# creating mask for ICA-MELODIC then applying band-pass filtering
Do_cmd fslmaths ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz -Tmean ${rsfmri_volum}/tmp.nii.gz
Do_cmd fbet ${rsfmri_volum}/tmp.nii.gz ${rsfmri_ICA}/func.nii.gz -m -n
Do_cmd mv ${rsfmri_ICA}/func_mask.nii.gz ${rsfmri_ICA}/mask.nii.gz
rm -fv ${rsfmri_volum}/tmp.nii.gz
rm -fv ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz
Do_cmd 3dTproject -input ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz -prefix ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz -passband 0.01 666
Do_cmd fslmaths ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz -Tmean ${rsfmri_volum}/${id}_singleecho_fmrispace_mean_orig.nii.gz
Do_cmd fslmaths ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz -Tmean ${tmp}/${id}_singleecho_fmrispace_mean.nii.gz

# run MELODIC for ICA-FIX
if [[ ! -f `find ${rsfmri_ICA}/filtered_func_data.ica/ -name "melodic_IC.nii.gz"` ]]; then
    Do_cmd cp ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz ${rsfmri_ICA}/filtered_func_data.nii.gz
    Do_cmd melodic --in=${rsfmri_ICA}/filtered_func_data.nii.gz \
                                    --tr=0.6 \
                                    --nobet \
                                    --mask=${rsfmri_ICA}/mask.nii.gz \
                                    --bgthreshold=3 \
                                    --mmthresh=0.5 \
                                    --report \
                                    --Oall \
                                    --outdir=${rsfmri_ICA}/filtered_func_data.ica \
                                    --Omean=${rsfmri_ICA}/mean_func.nii.gz
else
    Warning "!!!!  DON'T PANIC; SKIPPIN MELODIC  !!!!"
fi

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# WHAT ARE YOU REGISTERING???
# Register to Freesurfer space
Do_cmd bbregister --s ${id} --mov ${rsfmri_volum}/${id}_singleecho_fmrispace_mean_orig.nii.gz --reg ${dir_warp}/${id}_singleecho_fmri2fs.lta --o ${dir_warp}/${id}_singleecho_fmri2fs_outbbreg_FIX.nii.gz --init-fsl --bold
Do_cmd lta_convert --inlta ${dir_warp}/${id}_singleecho_fmri2fs.lta --outlta ${dir_warp}/${id}_singleecho_fs2fmri.lta --invert

mkdir ${rsfmri_ICA}/mc
mkdir ${rsfmri_ICA}/reg

Do_cmd cp ${rsfmri_volum}/${id}_singleecho.1D ${rsfmri_ICA}/mc/prefiltered_func_data_mcf.par
Do_cmd fslroi ${rsfmri_ICA}/filtered_func_data.nii.gz ${rsfmri_ICA}/reg/example_func.nii.gz 397 1
Do_cmd cp $T1nativepro ${rsfmri_ICA}/reg/highres.nii.gz
Do_cmd cp ${rsfmri_ICA}/filtered_func_data.ica/mean.nii.gz ${rsfmri_ICA}/mean_func.nii.gz

Do_cmd mri_concatenate_lta ${dir_warp}/${subject}_t1w_0.8mm_nativepro_t1w2fs.lta \
                           ${dir_warp}/${subject}_singleecho_fs2fmri.lta \
                           ${dir_warp}/${subject}_t1w_0.8mm_nativepro_t1w2fmri.lta

Do_cmd lta_convert --inlta ${dir_warp}/${subject}_t1w_0.8mm_nativepro_t1w2fmri.lta --outfsl ${rsfmri_ICA}/reg/highres2example_func.mat

Do_cmd mri_vol2vol --mov $T1nativepro \ # Is this really nativepro???
                   --targ ${tmp}/${id}_singleecho_fmrispace_mean.nii.gz \
                   --o ${rsfmri_ICA}/filtered_func_data.ica/t1w2fmri_brain.nii.gz \
                   --lta ${dir_warp}/${subject}_t1w_0.8mm_nativepro_t1w2fmri.lta

#------------------------------------------------------------------------------#
# run ICA-FIX if melodic has been run
if  [[ -f `find ${rsfmri_ICA}/filtered_func_data.ica/ -name "melodic_IC.nii.gz"` ]] && \
    [[ -f `find /host/fladgate/local_raid/MICA-MTL/ICAFIX_training/ -name "MICAMTL_training_15HC_15PX.RData"` ]]; then
    # if an error occurs during ICA-FIX, try the command line below to install required R packages and re-run
    # Rscript -e "for (i in c('kernlab','ROCR','class','party','e1071','randomForest'))install.packages(i)"
    Do_cmd /data_/mica1/01_programs/fix/fix ${rsfmri_ICA}/ \ # Need to add this to the main PATH
					    /host/fladgate/local_raid/MICA-MTL/ICAFIX_training/MICAMTL_training_15HC_15PX.RData \  # WHAT is this???
					    20 -m -h 100
else
    Warning "!!!!  NO MELODIC OUTPUT OR NO ICA-FIX TRAINING FILE  !!!!"
fi

#------------------------------------------------------------------------------#
# Change single echo files for clean ones
yes | Do_cmd cp -rf ${rsfmri_ICA}/filtered_func_data_clean.nii.gz ${tmp}/${id}_singleecho_fmrispace_HP.nii.gz
yes | Do_cmd cp -rf ${rsfmri_ICA}/filtered_func_data_clean.nii.gz ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz
Do_cmd fslmaths ${tmp}/${id}_singleecho_fmrispace_HP.nii.gz -Tmean ${tmp}/${id}_singleecho_fmrispace_mean.nii.gz

# Calculate tissue-specific signals
Do_cmd mri_concatenate_lta $struct2fs ${dir_warp}/${id}_singleecho_fs2fmri.lta ${tmp}/struct2fmri.lta
for idx in 0 1 2; do
    [[ $idx == 0 ]] && tissue=CSF
    [[ $idx == 1 ]] && tissue=GM
    [[ $idx == 2 ]] && tissue=WM
    tissuemap=$(find "$structuralDirectory" -name "*nativepro_pve_${idx}.nii.gz")
    Do_cmd mri_vol2vol --mov  "$tissuemap" --targ  ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz --lta ${tmp}/struct2fmri.lta --o ${tmp}/${id}_singleecho_"$tissue".nii.gz
    Do_cmd fslmeants -i ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz -o ${rsfmri_volum}/${id}_singleecho_"$tissue".txt -m ${tmp}/${id}_singleecho_"$tissue".nii.gz -w
done

Do_cmd fslmaths ${tmp}/${id}_singleecho_WM.nii.gz -add  ${tmp}/${id}_singleecho_GM.nii.gz -add  ${tmp}/${id}_singleecho_CSF.nii.gz ${tmp}/${id}_singleecho_WB.nii.gz
Do_cmd fslmeants -i ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz -o ${rsfmri_volum}/${id}_singleecho_global.txt -m ${tmp}/${id}_singleecho_WB.nii.gz -w

# Motion confound
Do_cmd fsl_motion_outliers -i ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz -o ${rsfmri_volum}/${id}_singleecho_spikeRegressors_REFRMS.1D -s ${rsfmri_volum}/${id}_singleecho_metric_REFRMS.1D --refrms --nomoco

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# register to MNI space???????
fmri_mean=${rsfmri_volum}/${id}_singleecho_fmrispace_mean.nii.gz
Do_cmd epi_reg --epi=$fmri_mean \
                      --t1=$T1_nativepro \
                      --t1brain=$T1nativepro_brain \
                      --out=${tmp}/${id}_fmri2t1

Do_cmd applywarp -i ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz \
                        -o ${rsfmri_volum}/${id}_singleecho_fmrispace_fmri2t12mni.nii.gz \
                        -w ${dir_warp}/${id}_native2nlmni_0.8mm.nii.gz \
                        -r ${mnitemplate3mm} \
                        --premat=${tmp}/${id}_fmri2t1.mat
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# Register to surface
for x in lh rh; do
    [[ $x == lh ]] && hemisphere=l || hemisphere=r
    HEMI=`echo $hemisphere | tr [:lower:] [:upper:]`
    # Map the non high-passed volumetric timeseries to the surface so we can compute tSNR
    Do_cmd mri_vol2surf \
        --mov ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz \
        --reg ${dir_warp}/${id}_singleecho_fmri2fs.lta \
        --projfrac-avg 0.2 0.8 0.1 \
        --trgsubject ${id} \
        --interp trilinear \
        --hemi ${x} \
        --out ${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}_NoHP.mgh

    # Map high-passed timeseries to surface - this is what will be used to generate the connectomes later
    Do_cmd mri_vol2surf \
        --mov ${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz \
        --reg ${dir_warp}/${id}_singleecho_fmri2fs.lta \
        --projfrac-avg 0.2 0.8 0.1 \
        --trgsubject ${id} \
        --interp trilinear \
        --hemi ${x} \
        --out ${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}.mgh
    Do_cmd mri_convert ${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}.mgh ${tmp}/${id}_singleecho_fmri2fs_${x}.func.gii

    # Register to conte69 and apply smooth
    Do_cmd wb_command -metric-resample \
        ${tmp}/${id}_singleecho_fmri2fs_${x}.func.gii \
        ${dir_conte69}/${id}_${hemisphere}h_sphereReg.surf.gii \
        ${util_surface}/fs_LR-deformed_to-fsaverage.${HEMI}.sphere.32k_fs_LR.surf.gii \
        ADAP_BARY_AREA \
        ${tmp}/${id}_singleecho_fmri2fs_${x}_c69-32k.func.gii \
        -area-surfs \
        ${dir_surf}/${id}/surf/${hemisphere}h.midthickness.surf.gii \
        ${dir_conte69}/${id}_${hemisphere}h_midthickness_32k_fs_LR.surf.gii
    Do_cmd wb_command -metric-smoothing \
        ${util_surface}/fsaverage.${HEMI}.midthickness_orig.32k_fs_LR.surf.gii \
        ${tmp}/${id}_singleecho_fmri2fs_${x}_c69-32k.func.gii \
        10 \
        ${tmp}/${id}_singleecho_fmri2fs_${x}_c69-32k_10mm.func.gii
    Do_cmd mri_convert ${tmp}/${id}_singleecho_fmri2fs_${x}_c69-32k.func.gii ${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}_c69-32k.mgh
    Do_cmd mri_convert ${tmp}/${id}_singleecho_fmri2fs_${x}_c69-32k_10mm.func.gii ${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}_c69-32k_10mm.mgh

done

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
:'  +++++++++++++++++ S U B C O R T E X ++++++++++++++++++'
# Flirt fsfast segmentation (nativepro) to rsfmri space using the inverse of .mat on line 323 (epi_reg)
Do_cmd convert_xfm -omat ${tmp}/${id}_t12fmri.mat -inverse ${tmp}/${id}_fmri2t1.mat

Do_cmd fsl5.0-flirt -in HC012_t1w_0.8mm_nativepro_all_fast_firstseg.nii.gz \              # FSFAST OUTPUT EXAMPLE
                    -ref ${fmri_mean} \    # MEAN RSFMRI FROM LINE 320
                    -applyxfm -init ${tmp}/${id}_t12fmri.mat \                            # FROM LINE 378
                    -out HC012_t1w_0.8mm_rsfmri_all_fast_firstseg.nii.gz                  # OUTPUTS FIRST SEG IN RSFMRI SPACE

# Extract subcortical timeseries (mean, one ts per first-segemented structure, exluding the brainstem (ew))
Do_cmd mri_segstats --i ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz \               # RSFMRI IN ITS OWN SPACE
                    --seg HC012_t1w_0.8mm_rsfmri_all_fast_firstseg.nii.gz \               # OUTPUT FROM COMMAND JUST BEFORE
                    --exclude 0 --exclude 16 \                                            # EXLUDE EVERYTHING BUT ACC, AMY, CAUD, HIP, PAL, PUT, THAL
                    --avgwf ${rsfmri_volum}/${id}_FIRST_rsfmri-timeseries.txt             # ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
:' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++'


:'  +++++++++++++++++ C E R E B E L L U M +++++++++++++++++'
Do_cmd fsl5.0-flirt -in <cereballar parcellation in nativepro> \                              # FSFAST OUTPUT EXAMPLE
                    -ref ${rsfmri_volum}/${id}_singleecho_fmrispace_mean_orig.nii.gz \        # MEAN RSFMRI FROM LINE 320
                    -applyxfm -init ${tmp}/${id}_t12fmri.mat \                                # FROM LINE 378
                    -out HC012_t1w_0.8mm_rsfmri_cerebellar.nii.gz                             # OUTPUTS FIRST SEG IN RSFMRI SPACE

# Extract subcortical timeseries (mean, one ts per first-segemented structure, exluding the brainstem (ew))
Do_cmd mri_segstats --i ${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz \                   # RSFMRI IN ITS OWN SPACE
                    --seg HC012_t1w_0.8mm_rsfmri_cerebellar.nii.gz \                          # OUTPUT FROM COMMAND JUST BEFORE
                    --exclude 0                                                               # EXLUDE EVERYTHING BUT CEREBELLUM
                    --avgwf ${rsfmri_volum}/${id}_cerebellar_rsfmri-timeseries.txt            # ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
:' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

#------------------------------------------------------------------------------#
# Clean temporary directory
Do_cmd rm -rfv $tmp

#------------------------------------------------------------------------------#
# run post-rsfmri
labelDirectory=${dir_surf}/${id}/label/
python $MICAPIPE/03_proc_rsfmri_post.py ${id} ${proc_rsfmri} ${labelDirectory}

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "rsfMRI processing and post processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\t\t\tlogs:${dir_logs}/post_MPC.txt"
# echo "${id}, proc_rsfmri, RUN, `whoami`, $(date), `printf "%0.3f\n" ${eri}`" >> ${out}/brain-proc.csv
