.. _restingstateproc:

.. title:: rs-fMRI

Resting-state fMRI processing
============================================================

This module performs all pre-processing of subject resting-state functional MRI (rsfMRI) scans, in preparation for the construction of functional connectomes. This pipeline is optimized to leverage spin-echo images with reverse phase encoding for distortion correction. For increased functionality, our pipeline can also handle protocols in which a single fieldmap is acquired.

.. image:: sankey_rsfmri.png
   :align: center

-proc_rsfmri
--------------------------------------------------------

.. admonition:: Prerequisites 🖐🏼

     You need to run ``-proc_structural``, ``-proc_freesurfer`` and ``-post_structural`` before this stage

.. image:: proc_rsfmri.png
   :scale: 70 %
   :align: center

.. tabs::

    .. tab:: Processing steps

            - Remove first five TRs and reorient to LPI
            - Reorients to standard
            - Perform motion correction within rsfMRI scan and provided fieldmaps by registering each volume to the scan's own average
            - Calculate motion outliers
            - Apply distortion correction on motion-corrected images
            - Calculates binary mask from motion and distortion corrected volume
            - High-pass filtering of rsfMRI timeseries to remove frequencies below 0.01Hz
            - Run Multivariate Exploratory Linear Optimized Decomposition into Independent Components (MELODIC) on filtered timeseries
            - Compute linear and non-linear registrations between fMRI and T1-nativepro space, as well as boundary-based registration between fMRI and native Freesurfer space
            - Run FMRIB's ICA-based Xnoiseifier (ICA-FIX) using specified training file. Note that if ICA-FIX is not found on the user's system, or if MELODIC failed, ICA-FIX will be skipped and further processing will be performed using high-pass filtered timeseries
            - Calculate global and tissue-specific signal from processed timeseries (can be used for global-signal regression)
            - Compute motion confounds matrix from processed timeseries, to remove effects of frames with large motion in the timeseries
            - Register processed timeseries to native cortical surface. Minimially pre-processed (i.e. motion and distortion corrected) timeseries are also registered to the native cortical surface to compute statistics such as temporal signal-to-noise
            - Surface-template registration of native surface timeseries (fsaverage5, conte69)
            - Native surface, fsaverage5, and conte69-mapped timeseries are each smoothed with a 10mm Gaussian kernel
            - Use previously computed registrations to align cerebellar and subcortical parcellation to fMRI space
            - Regress motion spikes from cerebellar, subcortical, and cortical timeseries in linear model
            - Concatenate cerebellar, subcortical, and parcellated cortical timeseries and compute correlation matrix



    .. tab:: Usage

        **Terminal:**

        .. parsed-literal::
            $ mica-pipe **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-proc_rsfmri**

        **Docker command:**

        .. parsed-literal::
            $ docker -proc_rsfmri

        **Optional arguments:**

        -mainScanStr          ``<str>`` String to manually identify the main scan for rsfMRI processing (eg. *func/sub-001_<mainScanStr>.nii.gz*).
                              Default string is *task-rest_acq-AP_bold*.
        -fmri_pe              ``<path>`` Full path to the main phase encoding scan for rsfMRI processing.
                              Default string is *task-rest_acq-APse_bold*.
        -fmri_rpe             ``<path>`` Full path to the reverse phase encoding scan for rsfMRI processing. If the path is wrong topup will be skipped!.
                              Default string is *task-rest_acq-PAse_bold*.
        -mainScanRun          ``<num>`` If multiple resting-state scan runs exist in the BIDS rawdata,
                              you must specify which scan to process using this flag (e.g. '1').
        -phaseReversalRun     ``<num>`` If multiple phase encoding runs exist in the BIDS directory (only main phase is checked),
                              you must specify which file to process using this flag (e.g. '1').
        -topupConfig          ``<path>`` Specify path to config file that should be used for distortion correction using topup.
                              Default is *${FSLDIR}/etc/flirtsch/b02b0_1.cnf*.
        -smoothWithWB         Specify this option to use workbench tools for surface-based smoothing (more memory intensive), The
                              default smoothing is performed with freesurfer tools: *mri_surf2surf*.
        -regress_WM_CSF       Specify this option to perform white matter and CSF signal regression of timeseries.
                              Default = no white matter and CSF signal regression.
        -GSR                  Specify this option to perform global signal regression of timeseries.
                              Default = no global regression.
        -noFIX                Specify this option to skip ICA-FIX processing.
                              Default = FIX runs with the default training file.
        -icafixTraining       ``<path>`` Path to specified ICA-FIX training file for nuisance signal regression
                              (file.RData). Default is *${MICAPIPE}/functions/MICAMTL_training_15HC_15PX.RData*.
        -regAffine            Specify this option to perform an Affine registration ONLY from rsfMRI to T1w.
                              Default is non linear registration using ANTs-SyN.
                              We recommend this option for rsfMRI acquisitions with low resolution and/or low SNR.
        -sesAnat              ``<str>`` If longitudinal data is provided, this flag allows to register the current *functional* session to the desired *anatomical* session
                              Default processing is independent.

        .. admonition:: Distortion correction ⚠️

                 If the script can't find the *fmri_rpe* (reverse phase encoding), or *fmri_pe* (phase encoding) images,
                 the distortion correction will be skipped. If you provided the path to the *fmri_pe* and *fmri_rpe* images,
                 ensured that the paths are correct!!! On the following table you can see the possible scenarios:

                 =========  ========  ======================
                       Inputs                 Output
                 -------------------  ----------------------
                 fmri_rpe   fmri_pe           topup
                 =========  ========  ======================
                    Yes        Yes    runs using pe and rpe
                    Yes        No     runs using main as pe
                    No         No     skipped
                 =========  ========  ======================

        .. admonition:: WARNING: ⚠️ Melodic and FIX ⚠️

                FIX and Melodic are used by default to remove nuisance variable signal. However our default parameters might not suit all databases.
                Our default training file used for FIX was trained in-house, on a subset of 30 participants.
                Scans were acquire on a 3T Siemens Magnetom Prisma-Fit equipped with a 64-channel head coil.
                rs-fMRI scans of 7 minutes were acquired using multiband accelerated 2D-BOLD echo-planar imaging
                (3mm isotropic voxels, TR=600ms, TE=30ms, flip angle=52°, FOV=240×240mm2, slice thickness=3mm, mb factor=6, echo spacing=0.54ms).
                If your acquisition parameters are similar, feel free to use the defaults options in ``-proc_rsfmri``.

                Otherwise we recommend you to `train your own dataset for FIX <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide#Training_datasets>`_,
                or do not use fix and try a different method for nuisance variable signal removal, e.g.:

                .. code-block:: bash
                   :caption: On the next example FIX and Melodic will be skipped, but global signal, white matter and CSF regressions will be applied:
                   :linenos:

                   mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory> \
                             -proc_rsfmri -noFIX -regress_WM_CSF -GSR


        .. admonition:: Longitudinal acquisitions ⚠️

                 If your database contains multiple sessions (different days) and you wish to register the rsfMRI to the same
                 anatomical session, you should use the ``-sesAnat`` flag. For example if you database looks like:

                 .. parsed-literal::
                     sub-01
                     ├── ses-func01
                     │   └── func
                     ├── ses-func02
                     │   └── func
                     └── ses-struct01
                         └── anat

                 You should specify the ``-sesAnat struct01`` to register each session to the same anatomical volume.

                 .. code-block:: bash
                    :caption: On the next example, func01 and func02 will be registered to the anatomical image in struct01:
                    :linenos:

                     mica-pipe -sub 01 -ses func01 -proc_rsfmri -sesAnat struct01 -bids <bidsDir> -out <outDir>
                     mica-pipe -sub 01 -ses func02 -proc_rsfmri -sesAnat struct01 -bids <bidsDir> -out <outDir>

    .. tab:: Outputs

        Directories created or populated by **-proc_rsfmri**:

        .. parsed-literal::

            - <outputDirectory>/micapipe/func
            - <outputDirectory>/micapipe/func/ICA_MELODIC
            - <outputDirectory>/micapipe/func/surfaces
            - <outputDirectory>/micapipe/func/volumetric
            - <outputDirectory>/micapipe/xfms

        Files generated by **-proc_rsfmri**:

        .. parsed-literal::
            - All outputs generated by MELODIC, or necessary for ICA-FiX, are stored in:
                *<outputDirectory>/micapipe/<sub>/func/ICA_MELODIC*

            - All volumetric processing outputs are stored in
                *<outputDirectory>/micapipe/<sub>/func/volumetric*

                - rsfMRI processing json card:
                    *<sub>_space-rsfmri_desc-singleecho_clean.json*

                - Motion confounds processing (<tag> = reversePhaseScan, mainPhaseScan):
                    *<sub>_space-rsfmri_<tag>.1D*

                - Motion outliers and metric values:
                    *<sub>_space-rsfmri_spikeRegressors_FD.1D*
                    *<sub>_space-rsfmri_metric_FD.1D*

                - Motion and distortion corrected image:
                    *<sub>_space-rsfmri_desc-singleecho.nii.gz*

                - Mean motion and distortion corrected image:
                    *<sub>_space-rsfmri_desc-singleecho_mean.nii.gz*

                - Skull-stipped mean motion and distortion corrected image:
                    *<sub>_space-rsfmri_desc-singleecho_brain.nii.gz*

                - High-passed, motion and distortion corrected image:
                    *<sub>_space-rsfmri_desc-singleecho_HP.nii.gz*

                - Nuisance-signal regressed timeseries (i.e. output of ICA-FIX):
                    *<sub>_space-rsfmri_desc-singleecho_clean.nii.gz*

                - Tissue-specific mean signal (<tissue> = CSF, GM, or WM):
                    *<sub>_space-rsfmri_pve_<tissue>.txt*

                - Global mean signal:
                    *<sub>_space-rsfmri_global.txt*

                - Subcortical segmentation in fMRI space:
                    *<sub>_space-rsfmri_desc-singleecho_subcortical.nii.gz*

                - Mean signal in each subcortical parcel:
                    *<sub>_space-rsfmri_desc-singleecho_timeseries_subcortical.txt*

                - Cerebellar segmentation in fMRI space:
                    *<sub>_space-rsfmri_desc-singleecho_cerebellum.nii.gz*

                - Mean signal in each cerebellar parcel:
                    *<sub>_space-rsfmri_desc-singleecho_timeseries_cerebellum.txt*

                - Parcel statistics for cerebellum, to screen for any missing parcels:
                    *<sub>_space-rsfmri_desc-singleecho_cerebellum_roi_stats.txt*


            - Vertexwise cortical timeseries (<hemi> = rh, lh)
                stored in *<outputDirectory>/micapipe/func/surfaces*:

                - Motion and distortion corrected timeseries mapped to native cortical surface:
                    *<sub>_rsfmri_space-fsnative_<hemi>_NoHP.mgh*

                - Fully pre-processed timeseries mapped to native cortical surface:
                    *<sub>_rsfmri_space-fsnative_<hemi>.mgh*
                    *<sub>_rsfmri_space-fsnative_<hemi>_10mm.mgh*

                - Timeseries mapped to fsaverage5 template:
                    *<sub>_rsfmri_space-fsaverage5_<hemi>.mgh*
                    *<sub>_rsfmri_space-fsaverage5_<hemi>_10mm.mgh*

                - Timeseries mapped to conte69 template:
                    *<sub>_rsfmri_space-conte69-32k_<hemi>.mgh*
                    *<sub>_rsfmri_space-conte69-32k_<hemi>_10mm.mgh*

                - Vertexwise and smoothed timeseries on conte69 template, following regression of motion spikes:
                    *<sub>_rsfmri_space-conte69-32k_desc-timeseries_clean.txt*

            - Temporal signal-to-noise ratio computed on native cortical surface from motion and distortion correction timesries:
                *<sub>_rsfmri_desc-tSNR.txt*

            - Functional connectome matrices (r values) generated from smoothed, parcellated timeseries sampled in subcortex, cerebellum, and cortical surface
               <parc> = up to 18 parcellations

                - Conte69 cortical surface:
                    *<sub>_rsfmri_space-conte69-32k_atlas-<parc>_desc-FC.txt*

                - Native cortical surface:
                    *<sub>_rsfmri_space-fsnative_atlas-<parc>_desc-FC.txt*

                - Contatenated timeseries sampled in subcortex, cerebellum, and parcellated native cortical surface models:
                    *<sub>_rsfmri_space-fsnative_atlas-<parc>_desc-timeseries.txt*

            - rsfMRI registration files are found in *<outputDirectory>/micapipe/<sub>/xfms*

                - Boundary based registration from rsfMRI space to native freesurfer space:
                    *<sub>_from-rsfmri_to-fsnative_bbr_outbbreg_FIX.nii.gz*
                    *<sub>_from-rsfmri_to-fsnative_bbr.dat*
                    *<sub>_from-rsfmri_to-fsnative_bbr.dat.log*
                    *<sub>_from-rsfmri_to-fsnative_bbr.dat.mincost*
                    *<sub>_from-rsfmri_to-fsnative_bbr.dat.param*
                    *<sub>_from-rsfmri_to-fsnative_bbr.dat.sum*

                - Affine registration between T1w nativepro and rsfmri space:
                    *<sub>_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat*
                    *<sub>_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_InverseWarped.nii.gz*
                    *<sub>_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_Warped.nii.gz*

                - Non-linear registrations between T1w in dwi space to wmNorm in dwi space:
                    *<sub>_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat*
                    *<sub>_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1InverseWarp.nii.gz*
                    *<sub>_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1Warp.nii.gz*
                    *<sub>_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_InverseWarped.nii.gz*
                    *<sub>_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_Warped.nii.gz*
