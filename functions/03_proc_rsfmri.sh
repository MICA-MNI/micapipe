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
# This workflow makes use of AFNI, FSL, ANTs, FIX, enigmatoolbox
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

#------------------------------------------------------------------------------#
Title "Running MICA rsfMRI processing"

# Assigns variables names
bids_variables $BIDS $id $out
# print the names on the terminal
bids_print.variables-rsfmri

# GLOBAL variables for this script
Info "ANTs will use $CORES CORES"
Info "wb_command will use $OMP_NUM_THREADS threads"

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

# Check inputs: rsfMRI and phase encoding
if [ ! -f ${mainScan} ]; then Error "Subject $id doesn't have acq-AP_bold: \n\t ${subject_bids}/func/"; exit; fi
if [ ! -f ${mainPhaseScan} ]; then Warning "Subject $id doesn't have acq-APse_bold: TOPUP will be skipped"; fi
if [ ! -f ${reversePhaseScan} ]; then Warning "Subject $id doesn't have acq-PAse_bold: TOPUP will be skipped"; fi
if [ ! -f ${T1nativepro} ]; then Error "Subject $id doesn't have T1_nativepro: run -proc_volumetric"; exit; fi

# Define directories
export SUBJECTS_DIR=$dir_surf

# rsfMRI directories
rsfmri_volum=${proc_rsfmri}/volumetric   # volumetricOutputDirectory
rsfmri_surf=${proc_rsfmri}/surfaces      # surfaceOutputDirectory
rsfmri_ICA=$proc_rsfmri/ICA_MELODIC      # ICAOutputDirectory

# Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script).
for x in $rsfmri_surf $rsfmri_volum $rsfmri_ICA; do
    [[ ! -d ${x} ]] && mkdir -p ${x}
done

#------------------------------------------------------------------------------#
# Begining of the REAL processing
# gettin dat readout time
readoutTime=`cat ${mainScan/.nii.gz/}.json | grep "TotalReadoutTime" | grep -Eo [0-9].[0-9]+`

# Scans to process
toProcess=($mainScan $mainPhaseScan $reversePhaseScan)
tags=(mainScan mainPhaseScan reversePhaseScan)
singleecho=${rsfmri_volum}/${id}_singleecho_fmrispace.nii.gz

# Loop over all scans for everything before motion correction across scans.
if [[ ! -f ${singleecho} ]]; then
    for i in {0..2}; do
        # Get basic parameters
        rawNifti=${toProcess[$i]}
        tag=${tags[$i]}

        # IF FILE NOT FOUND DON'T RUN
        if [ -f ${rawNifti} ]; then
              Info "Processing TAG-$tag scan, readout time: $readoutTime ms"
              Note "RAWNIFTI:" $rawNifti

              # Drop first five TRs and reorient (same orientation as T1nativepro)
              if [ "$tag" == "mainScan" ]; then
                  Do_cmd nifti_tool -cbl -prefix ${tmp}/${tag}_trDrop.nii.gz -infiles "$rawNifti"'[5..$]'
                  Do_cmd 3dresample -orient LPI -prefix ${tmp}/${tag}_reorient.nii.gz -inset ${tmp}/${tag}_trDrop.nii.gz
              else
                  Do_cmd 3dresample -orient LPI -prefix ${tmp}/${tag}_reorient.nii.gz -inset $rawNifti
              fi


              # Remove slices to make an even number of slices in all directions (requisite for topup).
              dimensions=`fslhd ${tmp}/${tag}_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
              newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
              Do_cmd fslroi ${tmp}/${tag}_reorient.nii.gz ${tmp}/${tag}_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` # I removed -1 -1

              # Motion correction within scans
              Do_cmd fslmaths ${tmp}/${tag}_sliceCut.nii.gz -Tmean ${tmp}/${tag}_sliceCutMean.nii.gz
              Do_cmd 3dvolreg -Fourier -twopass -base ${tmp}/${tag}_sliceCutMean.nii.gz \
                              -zpad 4 -prefix ${tmp}/${tag}_mc.nii.gz \
                              -1Dfile ${rsfmri_volum}/${id}_${tag}.1D \
                              ${tmp}/${tag}_sliceCut.nii.gz
              Do_cmd fslmaths ${tmp}/${tag}_mc.nii.gz -Tmean ${tmp}/${tag}_mcMean.nii.gz
        fi
    done
else
      Info "Subject ${id} has a singleecho_fmrispace processed"
fi

# Calculate slice timing if necessary <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Calculate motion outliers with FSL
if [[ ! -f ${rsfmri_volum}/${id}_singleecho.1D ]]; then
    Do_cmd fsl_motion_outliers -i ${tmp}/mainScan_sliceCut.nii.gz \
                               -o ${rsfmri_volum}/${id}_singleecho_spikeRegressors_FD.1D \
                               -s ${rsfmri_volum}/${id}_singleecho_metric_FD.1D --fd
    Do_cmd mv ${rsfmri_volum}/${id}_mainScan.1D ${rsfmri_volum}/${id}_singleecho.1D
else
    Info "Subject ${id} has a singleecho.1D with motion outliers"
fi

# Only do distortion correction if field maps were provided, if not then rename the scan to distortionCorrected (just to make the next lines of code easy).
if [ ! -f ${mainPhaseScan} ] || [ ! -f ${reversePhaseScan} ]; then
    Warning "No AP or PA acquisition was found, TOPUP will be skip!!!!!!!"
    status="NO-topup"
    Do_cmd mv -v ${tmp}/mainScan_mc.nii.gz ${singleecho}
else
    if [[ ! -f ${rsfmri_volum}/TOPUP.txt ]] && [[ ! -f ${singleecho} ]]; then
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
        Do_cmd applytopup --imain=${mainScan} --inindex=1 --datain=${tmp}/singleecho_topupDataIn.txt --topup=${tmp}/singleecho_topup --method=jac --out=${singleecho}
        # Check if it worked
        if [[ ! -f ${singleecho} ]]; then Error "Something went wrong with TOPUP check ${tmp} and log:\n\t\t${dir_logs}/proc_rsfmri.txt"; exit; fi
        echo "${singleecho}, TOPUP, `whoami`, $(date)" >> ${rsfmri_volum}/TOPUP.txt
        status="TOPUP"
    else
          Info "Subject ${id} has singleecho in fmrispace with TOPUP"
    fi
fi

#------------------------------------------------------------------------------#
Info "!!!!!  goin str8 to ICA-FIX yo  !!!!!"

fmri_mean=${rsfmri_volum}/${id}_singleecho_fmrispace_mean.nii.gz
fmri_mask=${rsfmri_ICA}/mask.nii.gz
fmri_HP=${rsfmri_volum}/${id}_singleecho_fmrispace_HP.nii.gz

# Creates a mask from the motion corrected time series
Do_cmd fslmaths ${singleecho} -Tmean ${rsfmri_volum}/tmp.nii.gz
Do_cmd bet ${rsfmri_volum}/tmp.nii.gz ${rsfmri_ICA}/func.nii.gz -m -n
Do_cmd mv ${rsfmri_ICA}/func_mask.nii.gz ${fmri_mask}
Do_cmd rm -fv ${rsfmri_volum}/tmp.nii.gz

# High-pass filter - Remove all frequencies EXCEPT those in the range
if [[ ! -f ${fmri_HP} ]] ; then
    Do_cmd 3dTproject -input ${singleecho} -prefix $fmri_HP -passband 0.01 666
else
    Info "Subject ${id} has High-pass filter"
fi

# Calculates the mean rsfMRI volume
Do_cmd fslmaths $singleecho -Tmean $fmri_mean

#------------------------------------------------------------------------------#
# run MELODIC for ICA-FIX
melodic_IC=${rsfmri_ICA}/filtered_func_data.ica/melodic_IC.nii.gz
fmri_filtered=${rsfmri_ICA}/filtered_func_data.nii.gz

# melodic will run ONLY if FIX is avaliable
if  [[ -f `which fix` ]]; then
      if [[ ! -f ${melodic_IC} ]]; then
          Do_cmd cp $fmri_HP $fmri_filtered
          Do_cmd melodic --in=${fmri_filtered} \
                          --tr=0.6 \
                          --nobet \
                          --mask=${fmri_mask} \
                          --bgthreshold=3 \
                          --mmthresh=0.5 \
                          --report \
                          --Oall \
                          --outdir=${rsfmri_ICA}/filtered_func_data.ica \
                          --Omean=${rsfmri_ICA}/mean_func.nii.gz
          if [[ -f ${melodic_IC} ]]; status="${status}/melodic"; fi
      else
          Info "Subject ${id} has MELODIC outputs"
      fi
fi

#------------------------------------------------------------------------------#
fmri_in_T1nativepro=${proc_struct}/${id}_singleecho_nativepro_brain.nii.gz
T1nativepro_in_fmri=${rsfmri_ICA}/filtered_func_data.ica/t1w2fmri_brain.nii.gz
str_rsfmri_affine=${dir_warp}/${id}_rsfmri_to_nativepro_
mat_rsfmri_affine=${str_rsfmri_affine}0GenericAffine.mat
fmri_brain=${rsfmri_volum}/${id}_singleecho_fmrispace_brain.nii.gz

# masked mean rsfMRI time series
Do_cmd fslmaths $fmri_mean -mul $fmri_mask $fmri_brain

if [[ ! -f ${fmri_in_T1nativepro} ]] ; then
    Info "Registering fmri space to nativepro"
    Do_cmd antsRegistrationSyN.sh -d 3 -f $T1nativepro_brain -m $fmri_brain -o $str_rsfmri_affine -t a -n $CORES -p d
    Do_cmd antsApplyTransforms -d 3 -i $fmri_brain -r $T1nativepro -t $mat_rsfmri_affine -o $fmri_in_T1nativepro -v -u int

    # t1-nativepro to fmri space
    Do_cmd antsApplyTransforms -d 3 -i $T1nativepro -r $fmri_brain -t [$mat_rsfmri_affine,1] -o ${T1nativepro_in_fmri} -v -u int

else
    Info "Subject ${id} has a rsfMRI in T1nativepro space"
fi

#------------------------------------------------------------------------------#
# Register rsfMRI to Freesurfer space with Freesurfer
fmri2fs_lta=${dir_warp}/${id}_singleecho_fmri2fs.lta
if [[ ! -f ${fmri2fs_lta} ]] ; then
  Info "Registering fmri to FreeSurfer space"
    Do_cmd bbregister --s $id --mov $fmri_mean --reg ${fmri2fs_lta} --o ${dir_warp}/${id}_singleecho_fmri2fs_outbbreg_FIX.nii.gz --init-fsl --bold
else
    Info "Subject ${id} has a lta transformation matrix from fmri to Freesurfer space"
fi

#------------------------------------------------------------------------------#
# run ICA-FIX IF melodic ran succesfully
fix_output=${rsfmri_ICA}/filtered_func_data_clean.nii.gz
fmri_processed=${rsfmri_volum}/${id}_singleecho_fmrispace_clean.nii.gz

if  [[ -f ${melodic_IC} ]] && [[ -f `which fix` ]]; then
    if [[ ! -f ${fix_output} ]] ; then
          Info "Getting ICA-FIX requirements"
          # FIX requirements
          if [ ! -d ${rsfmri_ICA}/mc ]; then mkdir ${rsfmri_ICA}/mc; fi
          if [ ! -d ${rsfmri_ICA}/reg ]; then mkdir ${rsfmri_ICA}/reg; fi

          Do_cmd cp ${rsfmri_volum}/${id}_singleecho.1D ${rsfmri_ICA}/mc/prefiltered_func_data_mcf.par
          Do_cmd fslroi ${fmri_filtered} ${rsfmri_ICA}/reg/example_func.nii.gz 397 1
          Do_cmd cp $T1nativepro ${rsfmri_ICA}/reg/highres.nii.gz
          Do_cmd cp ${rsfmri_ICA}/filtered_func_data.ica/mean.nii.gz ${rsfmri_ICA}/mean_func.nii.gz

          # REQUIRED by FIX ${rsfmri_ICA}/reg/highres2example_func.mat
          # Get transformation matrix T1native to rsfMRI space (ICA-FIX requirement)
          Do_cmd antsApplyTransforms  -v 1 -o Linear[${tmp}/highres2example_func.mat,0] -t [$mat_rsfmri_affine,1]

          # Transform matrix: ANTs (itk binary) to text
          Do_cmd ConvertTransformFile 3 ${tmp}/highres2example_func.mat ${tmp}/highres2example_func.txt

        # Fixing the transformations incompatiility between ANTS and FSL <<<<<<<<<<
          tmp_ants2fsl_mat=${tmp}/itk2fsl_highres2example_func.mat
          # Transform matrix: ITK text to matrix (FSL format)
          Do_cmd lta_convert --initk ${tmp}/highres2example_func.txt --outfsl $tmp_ants2fsl_mat --src $T1nativepro --trg $fmri_brain
          # apply transformation with FSL
          Do_cmd flirt -in $T1nativepro -out ${tmp}/t1w2fmri_brain_ants2fsl.nii.gz  -ref $fmri_brain -applyxfm -init $tmp_ants2fsl_mat
          # correct transformation matrix
          Do_cmd flirt -v -in ${tmp}/t1w2fmri_brain_ants2fsl.nii.gz -ref $T1nativepro_in_fmri -omat ${tmp}/ants2fsl_fixed.omat -cost mutualinfo -searchcost mutualinfo -dof 6
          # concatenate the matrices to fix the transformation matrix
          Do_cmd convert_xfm -concat ${tmp}/ants2fsl_fixed.omat -omat ${rsfmri_ICA}/reg/highres2example_func.mat $tmp_ants2fsl_mat


          Info "Running ICA-FIX"
          Do_cmd fix ${rsfmri_ICA}/ ${MICAPIPE}/functions/MICAMTL_training_15HC_15PX.RData 20 -m -h 100

          # Replace file if melodic ran correctly - Change single-echo files for clean ones
          if [[ ! -f ${fix_output} ]]; then
              yes | Do_cmd cp -rf $fix_output $fmri_processed
              status="${status}/FIX"
          else
              Error "melodic ran but ICA-FIX failed check log file:\n\t${dir_logs}/proc_rsfmri.txt"; exit
          fi
    else
          Info "Subject ${id} has filtered_func_data_clean from ICA-FIX already"
    fi
else
    Warning "!!!!  Melodic Failed and/or ICA-FIX was not found, check the software installation !!!!
                   If you've installed FIX try to install required R packages and re-run:
                   'kernlab','ROCR','class','party','e1071','randomForest'"
    Do_cmd cp -rf $fmri_HP $rsfmri_processed # OR cp -rf  $singleecho $rsfmri_processed <<<<<<<<<<<<<<<<<<<<< NOT SURE YET
    status="${status}/NO-fix"
fi

#------------------------------------------------------------------------------#
Info "Calculating tissue-specific and global signals changes"
global_signal=${rsfmri_volum}/${id}_singleecho_global.txt
if [[ ! -f ${global_signal} ]] ; then
      tissues=(CSF GM WM)
      for idx in {0..2}; do
           tissue=${tissues[$idx]}
           tissuemap=${proc_struct}/${id}_t1w_${res}mm_nativepro_brain_pve_${idx}.nii.gz
           tissue_series=${rsfmri_volum}/${id}_singleecho_${tissue}.txt
           if [[ ! -f ${tissue_series} ]] ; then
           Do_cmd antsApplyTransforms -d 3 -i $tissuemap -r $fmri_mean -t [$mat_rsfmri_affine,1] -o ${tmp}/${id}_singleecho_${tissue}.nii.gz -v -u int
           Do_cmd fslmeants -i $fmri_processed -o ${tissue_series} -m ${tmp}/${id}_singleecho_${tissue}.nii.gz -w
         else
             Info "Subject ${id} has $tissue time-series"
         fi
      done
      Do_cmd fslmaths ${tmp}/${id}_singleecho_WM.nii.gz -add  ${tmp}/${id}_singleecho_GM.nii.gz -add  ${tmp}/${id}_singleecho_CSF.nii.gz ${tmp}/${id}_singleecho_WB.nii.gz
      Do_cmd fslmeants -i $fmri_processed -o ${global_signal} -m ${tmp}/${id}_singleecho_WB.nii.gz -w
else
      Info "Subject ${id} has Global time-series"
fi

# Motion confound
spikeRegressors=${rsfmri_volum}/${id}_singleecho_spikeRegressors_REFRMS.1D
if [[ ! -f ${spikeRegressors} ]] ; then
    Do_cmd fsl_motion_outliers -i $fmri_processed -o ${spikeRegressors} -s ${rsfmri_volum}/${id}_singleecho_metric_REFRMS.1D --refmse --nomoco
else
    Info "Subject ${id} has a spike Regressors from fsl_motion_outliers"
fi

#------------------------------------------------------------------------------#
# Register to surface
for x in lh rh; do
    [[ $x == lh ]] && hemisphere=l || hemisphere=r
    HEMI=`echo $hemisphere | tr [:lower:] [:upper:]`
    out_surf=${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}_c69-32k_10mm.mgh
    if [[ ! -f ${out_surf} ]] ; then
          # Map the non high-passed volumetric timeseries to the surface so we can compute tSNR
          Do_cmd mri_vol2surf \
              --mov ${singleecho} \
              --reg ${fmri2fs_lta} \
              --projfrac-avg 0.2 0.8 0.1 \
              --trgsubject ${id} \
              --interp trilinear \
              --hemi ${x} \
              --out ${rsfmri_surf}/${id}_singleecho_fmri2fs_${x}_NoHP.mgh

          # Map high-passed timeseries to surface - this is what will be used to generate the connectomes later
          Do_cmd mri_vol2surf \
              --mov $fmri_processed \
              --reg ${fmri2fs_lta} \
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
          Do_cmd mri_convert ${tmp}/${id}_singleecho_fmri2fs_${x}_c69-32k_10mm.func.gii ${out_surf}
    else
          Info "Subject ${id} has a singleecho fmri2fs ${x}_conte69-32 10mm surface"
    fi
done

#------------------------------------------------------------------------------#
#                           S U B C O R T E X
# Subcortical segmentation (nativepro) to rsfmri space
T1_seg_subcortex=${dir_volum}/${id}_t1w_${res}mm_nativepro_subcortical.nii.gz
rsfmri_subcortex=${rsfmri_volum}/${id}_singleecho_fmrispace_subcortical.nii.gz
timese_subcortex=${rsfmri_volum}/${id}_singleecho_timeseries_subcortical.txt

if [[ ! -f ${timese_subcortex} ]] ; then
      Do_cmd antsApplyTransforms -d 3 -i $T1_seg_subcortex -r $fmri_mean -n GenericLabel -t [$mat_rsfmri_affine,1] -o $rsfmri_subcortex -v -u int
      # Extract subcortical timeseries
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i $fmri_processed --seg $rsfmri_subcortex --exclude 0 --exclude 16 --avgwf $timese_subcortex
else
      Info "Subject ${id} has rsfmri subcortical time-series"
fi

#------------------------------------------------------------------------------#
#                           C E R E B E L L U M
Info "Getting cerebellar timeseries"
T1_seg_cerebellum=${dir_volum}/${id}_t1w_${res}mm_nativepro_cerebellum.nii.gz
rsfmri_cerebellum=${rsfmri_volum}/${id}_singleecho_fmrispace_cerebellum.nii.gz
timese_cerebellum=${rsfmri_volum}/${id}_singleecho_timeseries_cerebellum.txt

if [[ ! -f ${timese_cerebellum} ]] ; then
      Do_cmd antsApplyTransforms -d 3 -i $T1_seg_cerebellum -r $fmri_mean -n GenericLabel -t [$mat_rsfmri_affine,1] -o $rsfmri_cerebellum -v -u int
      # Extract subcortical timeseries (mean, one ts per first-segemented structure, exluding the brainstem (ew))
      # Output: ascii text file with number of rows equal to the number of frames and number of columns equal to the number of segmentations reported
      Do_cmd mri_segstats --i $fmri_processed --seg $rsfmri_cerebellum --exclude 0 --avgwf ${timese_cerebellum}
else
      Info "Subject ${id} has rsfmri cerebellar time-series"
fi

#------------------------------------------------------------------------------#
# run post-rsfmri
Info "Running rsfMRI post processing"
labelDirectory=${dir_surf}/${id}/label/
python $MICAPIPE/functions/03_proc_rsfmri_post.py ${id} ${proc_rsfmri} ${labelDirectory}

#------------------------------------------------------------------------------#
# Clean temporary directory
Do_cmd rm -rfv $tmp

#------------------------------------------------------------------------------#
# QC notification of completition
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=`echo print $eri/60 | perl`

# Notification of completition
Title "rsfMRI processing and post processing ended in \033[38;5;220m `printf "%0.3f\n" ${eri}` minutes \033[38;5;141m:\n\t\tlogs:${dir_logs}/proc_rsfmri.txt"
echo "${id}, proc_rsfmri, ${status}, `whoami`, $(date), `printf "%0.3f\n" ${eri}`" >> ${out}/brain-proc.csv
