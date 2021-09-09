.. _faqA:

.. title:: Frequent Asked Questions

FAQ
================================================

Q: How do I register a volume from rsfMRI or DWI space to MNI152?
    **A:** although micapipe calculates all the necessary transformation matrices and warp fields to register from any native space to MNI152, this procedure is not implemented at any point in the pipeline.
    Using the following information you should be able to register between any spaces using ANTs.

    .. code-block:: bash
      :linenos:
      :caption: Map from MNI152 0.8mm space to DWI nativespace

      # map_mni2dwi
      function map_mni2dwi() {
      # map_dwi ROIinMNI.nii.gz BIDSPATH CASE OUTPATH
      outname=`basename ${1}`
      antsApplyTransforms -d 3 -i ${1} \
                -r ${2}/${3}/ses-01/dwi/${3}_ses-01_space-dwi_desc-b0.nii.gz \
                -n GenericLabel \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_1Warp.nii.gz \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_0GenericAffine.mat \
                -t [${2}/${3}/ses-01/xfm/${3}_ses-01_space-dwi_from-dwi_to-nativepro_mode-image_desc-0GenericAffine.mat,1] \
                -t [${2}/${3}/ses-01/xfm/${3}_ses-01_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat,1] \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1InverseWarp.nii.gz \
                -o ${4}/${3}_space-dwi_t1w_${outname} -v
      }

    .. code-block:: bash
      :linenos:
      :caption: Map from DWI native space to MNI152 0.8mm

      # map_dwi2mni
      function map_dwi2mni() {
      # map_dwi ROIinMNI.nii.gz BIDSPATH CASE OUTPATH
      outname=`basename ${1}`
      do_cmd antsApplyTransforms -d 3 -i ${1} \
                -r /data/mica3/boristest/ROIs/MNI152_T1_1mm.nii.gz \
                -n GenericLabel \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_space-dwi_from-dwi_to-nativepro_mode-image_desc-0GenericAffine.mat \
                -t [${2}/${3}/ses-01/xfm/${3}_ses-01_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_0GenericAffine.mat,1] \
                -t ${2}/${3}/ses-01/xfm/${3}_ses-01_space-dwi_from-dwi_to-dwi_mode-image_desc-SyN_1Warp.nii.gz \
                -o ${4}/${3}_mnispace_${outname} -v
      }

    .. code-block:: bash
      :linenos:
      :caption: From T1nativepro to rsfMRI space

      antsApplyTransforms -d 3 \
          -i anat/sub-CC110045_ses-pre_space-nativepro_t1w_brain.nii.gz \
          -r func/volumetric/sub-CC110045_ses-pre_space-rsfmri_desc-singleecho_brain.nii.gz \
          -t xfm/sub-CC110045_ses-pre_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1Warp.nii.gz \
          -t xfm/sub-CC110045_ses-pre_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat \
          -t [xfm/sub-CC110045_ses-pre_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat,1] \
          -o func/volumetric/sub-CC110045_ses-pre_space-rsfmri_t1w.nii.gz -v -u int

    .. code-block:: bash
      :linenos:
      :caption: From rsfMRI to T1nativepro space

      antsApplyTransforms -d 3 \
          -i func/volumetric/sub-CC110045_ses-pre_space-rsfmri_desc-singleecho_brain.nii.gz \
          -r anat/sub-CC110045_ses-pre_space-nativepro_desc-t1wbold.nii.gz \
          -t xfm/sub-CC110045_ses-pre_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat \
          -t [xfm/sub-CC110045_ses-pre_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat,1] \
          -t xfm/sub-CC110045_ses-pre_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1InverseWarp.nii.gz \
          -o anat/sub-CC110045_ses-pre_space-nativepro_desc-rsfmri_bold.nii.gz -v -u int


Q: How do I modify the smoothing on the surfaces?
    A: By default micapipe -Morphology only applies a smoothing of 10mm over the surfaces

Q: How do I process the rsfMRI If I have different runs in the same session?
    A: Right now the pipeline does not manages multiple runs on the same session. We recomend to concatenate all the runs and process the output file.

Q: Can I process multi-echo rsfMRI acquisition with micapipe?
    A: Right now micapipe can not handle in any way multi-echo acquisitions. Although is planned to be included in a future release.

Q: What if my database already has Freesurfer outputs with Quality Check?
    A:

Q: What if I have more than one set of DWI acquisitions?
    A:

Q: What if I already have DWI processed data?
    A:

Q: Why does micapipe not apply slice timing correction in rsfMRI processing?
    A: Micapipe was optimized for multi-band acquisitions which do not require slice timing correction.

Q: How do I train an ICA-FIX RData file?
    A:

Q: What if I already have rsfMRI processed data? (eg. from fmriPREP)
    A:

Q: Can I use a different Cortical / subcortical / cerebellar atlas  not included in the micapipe?
    A: At the present time this feature is not included in the latest release.

Q: If my database has the field maps instead of the reverse phase encoding acquisition of rsfMRI and DWI, can micapipe use them for the geometric distrotion correction?
    A: The usage of field maps ins not implemented in micapipe. However you could apply the field map based correction to the correspondent images and then run micapipe.
