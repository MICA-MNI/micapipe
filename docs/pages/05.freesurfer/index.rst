.. _freesurfer:

.. title:: Freesurfer processing

Freesurfer processing
======================================

The **Freesurfer** will allow you to generate the cortical surface segmentation of your data. Cortical surface segmentations were generated from native T1w scans using Freesurfer 6.0. A series of equivolumetric surfaces were also constructed for each participant between the pial and white matter boundaries. These surfaces were used for systematic sampling of qT1 image intensities, to compute individual microstructural profile similarity matrices (see next section). Here, qT1 images were co-registered to native FreeSurfer space of each participant using boundary-based registration. 

Freesurfer processing
--------------------------

.. tabs::

    .. tab:: What is it?
    
        **Freesurfer** processing description :

            - Recon-all processing
            - Segmentation of skull stripping, B1 bias field correction and gray-white matter
            - Reconstruction of cortical surface models
            - Nonlinear registration of the cortical surface 
            - Statistical analysis of group morphometry differences 
        
        For more informations about Freesurfer, see the `Freesurfer web site <https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki>`_
        

    .. tab:: How to use it? 

        **Freesurfer** command line : 

        .. parsed-literal:: 

            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-proc_freesurfer**
    
    .. tab:: What are the outcomes?

        **Freesurfer** processing outputs :

        .. parsed-literal:: 

         **proc-freesurfer directories**
            - /proc_struct/surfaces --> individual subject surfaces or surface-mapped structural features
            - /proc_struct/surfaces/sub-<subject_id> --> native Freesurfer space derivatives 
         **proc-freesurfer output files**
            - native surfaces
            - subject-specific parcellation annotation files
            - native midsurface for each hemispheres 

        All files in the /surfaces subdirectory consist of individual subject surfaces or surface-mapped structural features. These features are organized in four distinct subdirectories
        First, subject-specific outputs from Freesurfer (/surfaces/sub-HC#) include native Freesurfer space derivatives such as native surfaces and subject-specific parcellation annotation files, as well as native midsurface files for each hemispheres generated using workbench command tools. 
        
        Second, we provide pial, white matter and midsurfaces mapped to the Conte69 template space (/surfaces/conte69). This directory includes each subject’s native surface models mapped to the 32,000-vertex Conte69 surface template deformed to fsaverage space (indicated by “conte69-32k_fs_LR”), as well as the corresponding subject sphere resampled to the Conte69 template. All Conte69 template-mapped subject surfaces are made available as surface-type gifti files (.surf.gii).

   

Quality control
--------------------------
After the freesurfer processing, you will need to complete the quality control (QC) to corrected for any segmentation errors. 

    .. parsed-literal:: 
        
        **What you need for QC?**
        - Freesurfer 
        - Files path
             /proc_struct/surfaces/sub-<subject_id>/mri --> **brainmask.mgz**
            - /proc_struct/surfaces/sub-<subject_id>/surf --> **lh.pial lh.white rh.pial rh.white**
        - Add control point + manuel edits

