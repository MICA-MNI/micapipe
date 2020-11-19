.. _volumetric:

.. title:: Volumetric processing

Volumetric processing
======================================

Native structural images were anonymized and de-identified using AFNI’s reface plus tool. This tool aligns an input volume, in this case each participant’s first T1w scan, to a reference dataset, on which an average-face mask is applied. All other structural scans included in this release (e.g., additional T1w and qT1 images, and pre-processed T1w images) were non-linearly registered to the refaced T1w scan, and the resulting warp-field was subsequently applied to all structural images. This procedure allowed warping of the face mask to each structural scan while maintaining original image intensities and positioning of the brain in native space. The additional use of a brain mask applied to each structural scan during this procedure ensured that brain tissues were unaffected by the warp-field.

Each T1w scan was deobliqued and reoriented. Both scans were then linearly co-registered and averaged, automatically corrected for intensity nonuniformity, and intensity normalized. Resulting images were skull-stripped, and subcortical structures were segmented using FSL FIRST. Different tissue types (cortical and subcortical grey matter, white matter, cerebrospinal fluid) were segmented to perform anatomically constrained tractography. Pre-processed T1w images as well as segmented tissue types were co-registered to MNI152 using symmetric image normalization (SyN) implemented in Advanced Neuroimaging Tools (ANTs).


Proc-structural processing
--------------------------

.. tabs::

    .. tab:: What is it?
    
        **Proc-structural** processing description :

            - Create the t1w_nativepro for structural processing (Reorient to LPI with AFNI, If multiple T1w were provided --> Register and average to the first T1w, Coregistrer) 
            - Intensity non-uniform correction (N4)
            - Rescale intensity [100,0]
            - Brain mask (Brain extraction)
            - FSL first on the t1w_nativepro --> Subcortical 
            - FSL fast on the t1w_natipro --> Tissue types (Gray, White, CSF)
            - Registration to MNI space (No linear, 0.8 and 2mm)
            - Generate a five-tissue-type image for anatomically constrained tractography (5 Tissues Type (5TT))

    .. tab:: How to use it? 

        **Proc-structural** command line :

        .. parsed-literal:: 
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-proc_structural**    
    
    .. tab:: What are the outcomes?

        **Proc-structural** processing outputs :

        .. parsed-literal:: 

         **proc-structural directories**
            - /proc_struct
            - /proc_struct/first --> subcortical structure segmentations (FSL FIRST 51)
            - /xfms
         **proc-structural output files**
            - This directory contains two data spaces, **nativepro** and **MNI152**, indicated accordingly in all filenames. 
                - T1nativepro
                - T1nativepro_brain
                - T1nativepro_first
                - T1nativepro_5TT
                - T1_MNI152
                - T1_MNI152_brain
                - T1_nativepro_brain_to_0.8mm_MNI152
                - T1_nativepro_brain_to_2mm_MNI152


Post-structural processing
--------------------------
.. admonition:: Prerequisites 🖐🏼

     You need to run the *proc-structural* and *proc-freesurfer* before this stage

.. tabs::

    .. tab:: What is it?
    
        **Post-structural** processing description :

            - Registration (Compute affine matrix from Freesurfer space to nativepro)
            - Cerebellar parcellation from MNI152_0.8mm to T1_Nativepro (inverse transformation)
            - Create parcellation on nativepro space + Register the annot surface parcellation to the T1-freesurfer volume + Register parcellation (native freesurfer) to nativepro
            - Compute warp of native structural to Freesurfer and apply to 5TT and first (build the conte69-32k sphere + resample white and pial surfaces to conte69-32k)

    .. tab:: How to use it? 

        **Post-structural** command line :

        .. parsed-literal:: 
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-post_structural**    
    
    .. tab:: What are the outcomes?

        **Post-structural** processing outputs :

        .. parsed-literal:: 

         **post-structural directories**
            - /proc_struct
            - /proc_struct/volumetric --> parcellation schemes in volumetric nativepro space
         **post-structural output files**
            - This directory contains two data spaces, **nativepro** and **MNI152**, indicated accordingly in all filenames. 
                - T1_nativepro_glasser-360
                - T1_nativepro_schaefer
                - T1_nativepro_subcortical
                - T1_nativepro_vosdewael
            - Skull-stripped equivalent ?
            - Brain mask ?
            
