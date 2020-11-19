.. _dwi:

.. title:: DWI processing

Diffusion weighted images processing
======================================

Informations about the DWI processing...
DWI-Diffusion weighted images processing with MRtrix3

Proc-dwi processing
--------------------------

.. tabs::

    .. tab:: What is it?
    
        **Proc-dwi** processing description :

        - DWI denoise, bias filed correction and concatenation
        - DWI FSL preprocessing 
        - TOPUP ?
        - Registering corrected DWI-b0 to T1nativepro (transformation and inverse transformation)
        - Creating a DWI binary mask of processed volumes
        - Registering 5TT file to DWI-b0 space
        - Calculating basic DTI metrics
        - Calculating multi-shell multi-tissue, response function and fiber orientation distribution
        - QC of the tractography

    .. tab:: How to use it? 

        **Proc-dwi** command line :

         .. parsed-literal:: 
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-proc_dwi**    
    
    .. tab:: What are the outcomes?

        **Proc-dwi** processing outputs :

        .. parsed-literal:: 

         **proc-dwi directories**
            - /proc_dwi
         **proc-dwi output files**
            - dwi_dnsN4.mif
            - dwi_concatenate.mif
            - dwi_concatenate_denoised.mif
            - dwi_corr.mif

            - dwi_mask_nii.gz
            - dwi_5tt.nii.gz
            - dwi_b0.nii.gz
            - dwi_nativepro.nii.gz
            - t1w_dwi_nii.gz
            - dwi_to_nativepro_


Post-dwi processing
--------------------------
.. admonition:: Prerequisites 🖐🏼

     You need to run the *proc-dwi* before this stage

.. tabs::

    .. tab:: What is it?
    
        **Post-dwi** processing description :

        - Post tractography and connectome generation
        - Registering Cerebellua parcellation to DWI-b0 space
        - Registering Subcortical parcellation to DWI-b0 space
        - Building the tracts streamlines connectomes
        - Building the cortical-subcortical connectomes
        - Building the cortical connectomes 
        - Building the cortical-subcortical-cerebellar connectomes
        - Calculate the edge lenghts
        - QC of the tractogram 

    .. tab:: How to use it? 

        **Proc-dwi** command line :

         .. parsed-literal:: 
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-post_dwi**    
    
    .. tab:: What are the outcomes?

        **Post-dwi** processing outputs :

        .. parsed-literal:: 

         **post-dwi directories**
            - 
         **post-dwi output files**


