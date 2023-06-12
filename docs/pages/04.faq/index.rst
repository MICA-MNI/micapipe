.. _faqA:

.. title:: Frequent Asked Questions

******************
FAQ
******************

.. contents:: Table of Contents

Registration issues
================================================

Q: How do I register a volume from rsfMRI or DWI space to MNI152 and vice-versa?
    **A:** although micapipe calculates all the necessary transformation matrices and warp fields to register from any native space to MNI152,
    this procedure is not implemented at any point in the pipeline.

    Using the following information you should be able to register between any spaces using ANTs.
    The paths are relative to the subject directory (e.g. ``out/micapipe/sub-01/ses-01``)

    .. code-block:: bash
      :linenos:
      :caption: Map from MNI152 0.8mm space to DWI native space

      # Subjects identification
      subjectID=sub-01_ses-01

      # Map from MNI152 0.8mm to DWI space
      Input=file_in_space-MNI152_0.8mm.nii.gz
      Output=file_from_MNI152_0.8mm_in_space-dwi.nii.gz

      antsApplyTransforms -d 3 -i ${Input} \
                -r dwi/${subjectID}_space-dwi_desc-b0.nii.gz \
                -n GenericLabel \   # Interpolation!
                -t xfm/${subjectID}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_1Warp.nii.gz \
                -t xfm/${subjectID}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_0GenericAffine.mat \
                -t [xfm/${subjectID}_space-dwi_from-dwi_to-nativepro_mode-image_desc-0GenericAffine.mat,1] \
                -t [xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat,1] \
                -t xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1InverseWarp.nii.gz \
                -o ${Output} -v

    .. code-block:: bash
      :linenos:
      :caption: Map from DWI native space to MNI152 0.8mm

      # Map from DWI space to MNI152 0.8mm
      Input=file_in_space-dwi.nii.gz
      Output=file_from_dwi_in_space-MNI152_0.8mm.nii.gz

      antsApplyTransforms -d 3 -i ${Input} \
                -r ${MICAPIPE_DIR}/MNI152Volumes/MNI152_T1_0.8mm_brain.nii.gz \
                -n GenericLabel \
                -t xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz \
                -t xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat \
                -t xfm/${subjectID}_space-dwi_from-dwi_to-nativepro_mode-image_desc-0GenericAffine.mat \
                -t [xfm/${subjectID}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_0GenericAffine.mat,1] \
                -t xfm/${subjectID}_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_1Warp.nii.gz \
                -o ${Output} -v
      }

    .. code-block:: bash
      :linenos:
      :caption: From MNI152 to rsfMRI space

      # Map from MNI152 0.8mm to rsfMRI space
      Input=file_in_space-MNI152_0.8mm.nii.gz
      Output=file_from_MNI152_0.8mm_in_space-rsfMRI.nii.gz

      antsApplyTransforms -d 3 \
          -i ${Input} \
          -r func/volumetric/${subjectID}_space-rsfmri_desc-singleecho_brain.nii.gz \
          -t xfm/${subjectID}_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1Warp.nii.gz \
          -t xfm/${subjectID}_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat \
          -t [xfm/${subjectID}_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat,1] \
          -t [xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat,1] \
          -t xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1InverseWarp.nii.gz \
          -o ${Output} -v

    .. code-block:: bash
      :linenos:
      :caption: From rsfMRI to MNI152 space

      # Map from rsfMRI space to MNI152 0.8mm
      Input=file_in_space-rsfMRI.nii.gz
      Output=file_from_rsfMRI_in_space-MNI152_0.8mm.nii.gz

      antsApplyTransforms -d 3 \
          -i func/volumetric/${subjectID}_space-rsfmri_desc-singleecho_brain.nii.gz \
          -r ${MICAPIPE_DIR}/MNI152Volumes/MNI152_T1_0.8mm_brain.nii.gz \
          -t xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz \
          -t xfm/${subjectID}_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat \
          -t xfm/${subjectID}_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat \
          -t [xfm/${subjectID}_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat,1] \
          -t xfm/${subjectID}_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1InverseWarp.nii.gz \
          -o ${Output} -v -u int

    .. code-block:: bash
      :linenos:
      :caption: From rsfMRI to nativepro space

      # Map from rsfMRI space to nativepro 2mm
      Input=file_in_space-rsfMRI.nii.gz
      Output=file_from_rsfMRI_in_space-nativepro_2mm.nii.gz

      antsApplyTransforms -d 3 \
          -i func/volumetric/${subjectID}_space-rsfmri_desc-singleecho_brain.nii.gz \
          -r ${MICAPIPE_DIR}/MNI152Volumes/MNI152_T1_2mm_brain.nii.gz \
          -t xfm/${subjectID}_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat \
          -t [xfm/${subjectID}_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat,1] \
          -t xfm/${subjectID}_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1InverseWarp.nii.gz \
          -o ${Output} -v -u int
          
Surface issues
================================================

Q: What if my database already has quality checked *Freesurfer* outputs?
    **A:** If you have an already processed and QC-ed freesurfer directory remember to use the flag ``-freesurfer_dir`` when you run the pipeline!!
    It will make a copy of your data with a compatible naming inside ``out/freesurfer/``. It's up to you to keep the original *freesurfer* directory or erase it.

    .. code-block:: bash
       :caption:  Usage example of -freesurfer_dir flag
       :linenos:

        # Run micapipe
        mica-pipe -bids rawdata -out derivatives -sub 01 \
              -proc_freesurfer -freesurfer_dir <path_to_my_subject_freesurfer_dir> \


Q: How do I modify the smoothing on the surfaces?
    **A:** By default micapipe ``-Morphology`` only applies a smoothing of 10mm over the surfaces.
    If you need a different FWHM you can use either: ``mri_surf2surf`` from *freesurfer* or ``wb_command -metric-smoothing`` from *workbench*.

    In the next examples we'll use the subject ``sub-01`` session ``ses-01``:

    .. code-block:: bash
       :caption:  20mm smoothing of the left hemisphere using *freesurfer* tools
       :linenos:

        # OutDir is the directory with -Morphology outpus
        outDir=out/micapipe/sub-01/ses01/anat/surfaces/morphology

        # Declare the micapipe's freesurfer directory variable
        export SUBJECTS_DIR=out/freesurfer

        mri_surf2surf --hemi lh \
            --fwhm-trg 20 \
            --srcsubject sub-01_ses01 \
            --srcsurfval "${outDir}/sub-01_ses01_space-fsnative_desc-lh_thickness.mgh" \
            --trgsubject fsaverage5 \
            --trgsurfval "${outDir}/sub-01_ses01_space-fsaverage5_desc-lh_thickness_20mm.mgh"
            "${outDir}/lh_curv_20mm_c69-32k.func.gii"

    .. code-block:: bash
       :caption:  20mm smoothing of the left hemisphere using *WorkBench* tools
       :linenos:

       MICAPIPE_DIR=<path to the micapipe repository>

       # For WorkBench the first step is to convent the mgh surface file to GIFTI
        mri_convert "${outDir}/sub-01_ses01_space-conte69-32k_desc-lh_thickness.mgh" "/tmp/lh_curv_c69-32k_thickness.func.gii"

        wb_command -metric-smoothing \
            "${MICAPIPE_DIR}/surfaces/fsaverage.L.midthickness_orig.32k_fs_LR.surf.gii" \
            "/tmp/lh_curv_c69-32k_thickness.func.gii" \
            20 \
            "/tmp/lh_curv_20mm_c69-32k.func.gii"     # This is the 20 mm FWHM surface

        # Convert from GIFTI back to MGH
        mri_convert "/tmp/lh_curv_10mm_c69-32k.func.gii" "${outDir}/sub-01_ses01_space-conte69-32k_desc-lh_thickness_20mm.mgh"


Resting state issues
================================================

Q: How do I process multiple rsfMRI If I have different runs in the same session?
    **A:** Right now the pipeline does not manage multiple runs on the same session.
    We recommend to concatenate all the runs into one file and process it.
    Using ``mrcat`` from MRtrix3 for example:

    .. code-block:: bash
       :caption:  Concatenates all rsfMRI runs and process the output
       :linenos:

        # Inside the func directory contatenate all the runs
        mrcat sub-01_task-rest_run-1_bold.nii.gz sub-01_task-rest_run-2_bold.nii.gz sub-01_task-rest_desc-cat_bold.nii.gz

        # Copy the json file from run-1
        cp sub-01_task-rest_run-1_bold.json sub-01_task-rest_desc-cat_bold.json

        # Run micapipe and specify the name of the concatenated rsfmri
        mica-pipe -bids rawdata -out derivatives -sub 01 \
              -proc_rsfmri -mainScanStr task-rest_desc-cat_bold \


Q: Can I process multi-echo rsfMRI acquisition with micapipe?
    **A:** Although is planned to be included in a future release, right now `micapipe` cannot handle in any way multi-echo acquisitions.

Q: How do I train an ICA-FIX ``RData`` file?
    **A:** The default processing of rsfMRI was optimized with the MICs dataset. Thus, we generated a custom training file in order to use FIX with our dataset.
    If you want to use FIX to clean the noisy components from your rsfMRI time series, we recommend to train your own weights file using `FSL-FIX instructions <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide#Training_datasets>`_

Q: What if I already have rsfMRI processed data? (eg. from fmriPREP)
    **A:** Unfortunately, due to the multiple nested steps which take place inside ``proc_rsfmri`` is not possible to use an already processed rsfMRI volume.
    If you are thrilled to develop its implementation, please help us!

Q: Why does `micapipe` not apply slice timing correction in rsfMRI processing?
    **A:** `micapipe` was optimized for multi-band acquisitions, which do not require slice timing correction.

Q: If my database has the field maps instead of the reverse phase encoding acquisition of rsfMRI and DWI, can micapipe use them for the geometric distortion correction?
    **A:** The usage of field maps is not implemented in *micapipe*.
    However you could apply the field map based correction to their correspondent images and then run micapipe, choosing the fieldmap-based corrected images to process.


DWI issues
================================================

Q: What if I already have DWI processed data?
    **A:** If you have your DWI already processed and you don't want to run this step again. you can use the flag ``-dwi_processed`` to the processed DWI,
    it must be in ``mif`` format with *bvals*, *bvecs*, *PhaseEncodingDirection* and *ReadoutTime* encoded.

    .. code-block:: bash
       :caption:  Usage example of -dwi_processed flag
       :linenos:

        # Run micapipe
        mica-pipe -bids rawdata -out derivatives -sub 01 \
              -proc_dwi -dwi_processed <path>/dwi_preprocessed.mif \

Q: Can I save the tractogram (tck file) generated in ``-SC`` ?
    **A:** Yes, you need to use the flag ``-keep_tck`` to keep the connectome generated by the pipeline, which is erased by default.

Q: Does ``micapipe`` compress the tractogram (tck file) generated in ``-SC`` when they are saved with ``-keep_tck``?
    **A:** No, by default, our pipeline does not save the tractograms or compress them after connectome calculations (``-SC``). The main reason is storage and computational time, which is a significant issue even when compression is used (for 10 million streamlines, uncompressed/compressed tractograms are ~14G/3G per case). When running the pipeline, it is possible to keep tractograms (optional argument ``-keep_tck``), but it’s up to the final user to compress them.
    We included a short tutorial regarding this issue in the section `How to downsample a tractogram <../04.tckdownsample/index.html>`_

Q: My dataset contains multiple and separated DWI shells, can I process them individually or should I merge them before ``-proc_dwi``?
    **A:** The default behavior of micapipe is to use all the dwi shells provided and perform a rigid registration before concatenating them. While we found this step to increase robustness of our processing for some cases in the MICs dataset, applying a rotation prior to invoking ``dwipreproc`` (which runs eddy/topup and thus also corrects for motion within a given ``mif`` file) may increase data interpolation between shells.
    If the user wants to avoid this rigid registration between shells is possible to concatenate the shells prior processing and call the resulting file with their corresponding updated json, bvecs and bval files. We recommend using ``mrcat`` and ``mrconvert`` to perform this step.

    If the user wants to process only an individual shell from a dataset acquisition is also possible using the flag ``-dwi_main`` to select only one DWI nifti file.

    **Processing multiple shells with default parameters**

    .. code-block:: bash
       :linenos:

       # Path relative to main MICs rawdata directory
       sub_id="sub-HC001/ses-01/dwi/sub-HC001_ses-01"

       # Processing command
       micapipe -sub HC001 -ses 01 -proc_dwi \
         -bids . \
         -out ../derivatives \
         -dwi_main ${sub_id}_acq-b2000-91_dir-AP_dwi.nii.gz,${sub_id}_acq-b700-41_dir-AP_dwi.nii.gz \
         -dwi_rpe ${sub_id}_dir-PA_dwi.nii.gz \
         -dwi_acq multi-shell


    **Processing a single shell**

    .. code-block:: bash
       :linenos:

       # Processing command
       micapipe -sub HC001 -ses 01 -proc_dwi \
         -bids . \
         -out ../derivatives \
         -dwi_main ${sub_id}_acq-b700-41_dir-AP_dwi.nii.gz \
         -dwi_rpe ${sub_id}_dir-PA_dwi.nii.gz \
         -dwi_acq shell-b700

    **Concatenating shells before processing with micapipe**

    .. code-block:: bash
        :linenos:

        # From the subject's DWI directory
        cd sub-HC001/ses-01/dwi/
        sub="sub-HC001_ses-01"

        # Convert each shell to mif to store the associated bvecs, bvals and json files
        # Shell b2000
        mrconvert -json_import ${sub}_acq-b2000-91_dir-AP_dwi.json \
        -fslgrad ${sub}_acq-b2000-91_dir-AP_dwi.bvec ${sub}_acq-b2000-91_dir-AP_dwi.bval \
        ${sub}_acq-b2000-91_dir-AP_dwi.nii.gz \
        ${sub}_acq-b2000-91_dir-AP_dwi.mif

        # Shell b700
        mrconvert -json_import ${sub}_acq-b700-41_dir-AP_dwi.json \
        -fslgrad ${sub}_acq-b700-41_dir-AP_dwi.bvec ${sub}_acq-b700-41_dir-AP_dwi.bval \
        ${sub}_acq-b700-41_dir-AP_dwi.nii.gz \
        ${sub}_acq-b700-41_dir-AP_dwi.mif

        # Concatenate shells and export bvec, bval and json files
        mrcat ${sub}_acq-b2000-91_dir-AP_dwi.mif ${sub}_acq-b700-41_dir-AP_dwi.mif ${sub}_acq-MultiShell_dir-AP_dwi.mif

        # Convert mif to nifti
        mrconvert -export_grad_fsl ${sub}_acq-MultiShell_dir-AP_dwi.bvec ${sub}_acq-MultiShell_dir-AP_dwi.bval \
        -json_export ${sub}_acq-MultiShell_dir-AP_dwi.json \
        ${sub}_acq-MultiShell_dir-AP_dwi.mif \
        ${sub}_acq-MultiShell_dir-AP_dwi.nii.gz

        # Processing command
        micapipe -sub HC001 -ses 01 -proc_dwi \
        -bids \
        -out \
        -dwi_main ${sub}_acq-MultiShell_dir-AP_dwi.nii.gz \
        -dwi_rpe ${sub}_dir-PA_dwi.nii.gz
        -dwi_acq MultiShell


Parcellation issues
================================================

Q: Can I use a different cortical / subcortical / cerebellar atlas not included in the micapipe?
    **A:** At the present moment this feature is not included. If you wan't to help us implementing this new feature you are very welcome.
