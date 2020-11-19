.. _restingstate:

.. title:: Resting state processing

Resting state fMRI processing
======================================

Informations about the resting state processing...

Proc-rsfmri processing
--------------------------
.. admonition:: Prerequisites 🖐🏼

     You need to run the *proc-structural* and *proc-freesurfer* before this stage

.. tabs::

    .. tab:: What is it?
    
        **Proc-rsfmri** processing description :

            - Processing single --> Drop first five TRs and reorient (same orientation as T1nativepro) + Motion correction within scans
            - Calculate motion outliers with FSL 
            - Distorsion correction with TOPUP
            - ICA-FIX preparation --> calculate the mean rsfmri volume, create a mask, masked mean rsfmri time series
            - High-pass filter - Remove all frequencies except those in the range
            - Run MELODIC for ICA-FIX --> generate the independents components (ICA)
            - Registration to nativepro + inverse transformation (nativepro to fmri space)
            - Registration rsfmri to Freesurfer space (line 292)

    .. tab:: How to use it? 

        **Proc-rsfmri** command line :

        .. parsed-literal:: 
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-proc_rsfmri**    
    
    .. tab:: What are the outcomes?

        **Proc-rsfmri** processing outputs :
        
        .. parsed-literal:: 
            **rsfmri directories**
                - /proc_rsfmri/volumetric
                - /proc_rsfmri/surfaces
                - /proc_rsfmri/ICA_MELODIC
            **rsfmri output files**
            - 
