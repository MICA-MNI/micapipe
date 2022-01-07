.. _databases:

.. title:: Processing databases

Processing databases with ``micapipe``
============================================================

.. contents:: Table of Contents

In this section, you will find examples of how ``mica-pipe`` optional arguments where used to process a variety of databases with different acquisitions parameters.
Under the tab ``database structure`` the file organization of each dataset is listed. However you'll only find those files relevant for the pipeline, not the full structure.
Further information about the datasets such as the *source*, and *references* can be found here as well.

Microstructure-Informed Connectomics (MICs)
--------------------------------------------------------

.. figure:: mics_logo.png
    :align: left
    :scale: 20 %

Reference
   Royer, J., Rodriguez-Cruces, R., Tavakol, S., Lariviere, S., Herholz, P., Li, Q, Vos de Wael, R., Paquola, C., Benkarim, O., Park, B., Lowe, A. J., Margulies, D., Smallwood, J., Bernasconi, A., Bernasconi, N., Frauscher, B., Bernhardt, B. C.
   *An Open MRI Dataset for Multiscale Neuroscience*. bioRxiv 2021.08.04.454795

+--------+----------------------------------------------------------------------------------+
| Source | `CONP Portal: MICA-MICs <https://portal.conp.ca/dataset?id=projects/mica-mics>`_ |
+--------+----------------------------------------------------------------------------------+
| doi    | https://doi.org/10.1101/2021.08.04.454795                                        |
+--------+----------------------------------------------------------------------------------+

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            sub-HC001
            └── ses-01
                ├── anat
                │   ├── sub-HC001_ses-01_acq-inv1_T1map.json
                │   ├── sub-HC001_ses-01_acq-inv1_T1map.nii.gz
                │   ├── sub-HC001_ses-01_acq-mp2rage_T1map.json
                │   ├── sub-HC001_ses-01_acq-mp2rage_T1map.nii.gz
                │   ├── sub-HC001_ses-01_T1w.json
                │   └── sub-HC001_ses-01_T1w.nii.gz
                ├── dwi
                │   ├── sub-HC001_ses-01_acq-b2000-91_dir-AP_dwi.bval
                │   ├── sub-HC001_ses-01_acq-b2000-91_dir-AP_dwi.bvec
                │   ├── sub-HC001_ses-01_acq-b2000-91_dir-AP_dwi.json
                │   ├── sub-HC001_ses-01_acq-b2000-91_dir-AP_dwi.nii.gz
                │   ├── sub-HC001_ses-01_acq-b300-11_dir-AP_dwi.bval
                │   ├── sub-HC001_ses-01_acq-b300-11_dir-AP_dwi.bvec
                │   ├── sub-HC001_ses-01_acq-b300-11_dir-AP_dwi.json
                │   ├── sub-HC001_ses-01_acq-b300-11_dir-AP_dwi.nii.gz
                │   ├── sub-HC001_ses-01_acq-b700-41_dir-AP_dwi.bval
                │   ├── sub-HC001_ses-01_acq-b700-41_dir-AP_dwi.bvec
                │   ├── sub-HC001_ses-01_acq-b700-41_dir-AP_dwi.json
                │   ├── sub-HC001_ses-01_acq-b700-41_dir-AP_dwi.nii.gz
                │   ├── sub-HC001_ses-01_dir-PA_dwi.bval
                │   ├── sub-HC001_ses-01_dir-PA_dwi.bvec
                │   ├── sub-HC001_ses-01_dir-PA_dwi.json
                │   └── sub-HC001_ses-01_dir-PA_dwi.nii.gz
                └── func
                    ├── sub-HC001_ses-01_task-rest_acq-AP_bold.json
                    ├── sub-HC001_ses-01_task-rest_acq-AP_bold.nii.gz
                    ├── sub-HC001_ses-01_task-rest_acq-APse_bold.json
                    ├── sub-HC001_ses-01_task-rest_acq-APse_bold.nii.gz
                    ├── sub-HC001_ses-01_task-rest_acq-PAse_bold.json
                    └── sub-HC001_ses-01_task-rest_acq-PAse_bold.nii.gz


    .. tab:: Usage

        MICs and micapipe were developed practically at the same time, therefore they are fully compatible and the default arguments are enough to process it.
        Modalities include high-resolution anatomical (T1-weighted), microstructurally-sensitive (quantitative T1), diffusion-weighted (three shells, 140 directions in total) and resting-state functional imaging.
        DWI and rsfMRI include a reverse phase encoding acquisition for geometric distortion correction.

        .. code-block:: bash
           :linenos:

            mica-pipe -bids rawdata -out derivatives -sub HC001 -ses 01 -all

        .. code-block:: bash
           :caption:  Using the -all flag is the same as using all the processing flags!
           :linenos:

            mica-pipe -bids rawdata -out derivatives -sub HC001 -ses 01 \
                      -proc_structural \
                      -proc_freesurfer \
                      -post_structural \
                      -proc_dwi \
                      -SC \
                      -proc_rsfmri \
                      -GD \
                      -Morphology \
                      -MPC \
                      -QC_subj

        proc_structural
           This section was processed with the default parameters (T1w = ``sub-HC001_ses-01_T1w.nii.gz``).

        proc_freesurfer
           This section was processed with the default parameters (T1w = ``sub-HC001_ses-01_T1w.nii.gz``).

        post_structural
           This section was processed with the default parameters.

        Morphology
           This section was processed with the default parameters.

        GD
           This section was processed with the default parameters.

        proc_dwi
           This section was processed with the default parameters.

           DWI shells to process: ``acq-b2000-91_dir-AP_dwi``, ``acq-b300-11_dir-AP_dwi`` and ``acq-b700-41_dir-AP_dwi``.

           Reverse phase encoding image = ``sub-HC001_ses-01_dir-PA_dwi.nii.gz``.

        SC
           This section was processed with the default parameters.

        proc_rsfmri
           This section was processed with the default parameters. Using melodic and FIX for nuisance regression and non linear registration (SyN) to T1w-nativepro space.

           Main scan = ``task-rest_acq-AP_bold``

           Main phase encoding = ``task-rest_acq-APse``

           Reverse phase encoding = ``task-rest_acq-PAse``

        MPC
           This section was processed with the default parameters.

           Microstructural image = ``acq-mp2rage_T1map``

           Image for registration = ``acq-inv1_T1map``.

Epilepsy and Cognition (EpiC-UNAM)
--------------------------------------------------------

Reference
   Rodríguez-Cruces R., Bernhardt B. C., & Concha L. (2020). *Multidimensional associations between cognition and connectome organization in temporal lobe epilepsy.* NeuroImage, Volume 213, June 2020, 116706.

+--------+--------------------------------------------------+
| Source | Unpublished                                      |
+--------+--------------------------------------------------+
| doi    | https://doi.org/10.1016/j.neuroimage.2020.116706 |
+--------+--------------------------------------------------+

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            sub-00367
            └── ses-01
                ├── anat
                │   ├── sub-00367_ses-01_T1w.json
                │   └── sub-00367_ses-01_T1w.nii.gz
                ├── dwi
                │   ├── sub-00367_ses-01_acq-b2000_dir-AP_dwi.bval
                │   ├── sub-00367_ses-01_acq-b2000_dir-AP_dwi.bvec
                │   ├── sub-00367_ses-01_acq-b2000_dir-AP_dwi.json
                │   ├── sub-00367_ses-01_acq-b2000_dir-AP_dwi.nii.gz
                │   ├── sub-00367_ses-01_dir-PA_dwi.bval
                │   ├── sub-00367_ses-01_dir-PA_dwi.bvec
                │   ├── sub-00367_ses-01_dir-PA_dwi.json
                │   └── sub-00367_ses-01_dir-PA_dwi.nii.gz
                └── func
                    ├── sub-00367_ses-01_task-rest_bold.json
                    └── sub-00367_ses-01_task-rest_bold.nii.gz


    .. tab:: Usage

        This dataset contains T1 weighted images, a DWI single shell (b=2000, 60 directions) with a reverse phase encoding acquisition (``dir-PA_dwi``),
        and a 5 minutes length resting state fMRI.

        .. code-block:: bash
           :linenos:

            mica-pipe -bids rawdata -out derivatives -sub 00367 -ses 01
                      -proc_structural \
                      -proc_freesurfer \
                        -freesurfer_dir ./freesurfer_processed_00367 \
                      -post_structural \
                      -proc_dwi \
                      -SC \
                      -proc_rsfmri \
                        -mainScanStr task-rest_bold \
                        -regress_WM_CSF \
                        -noFIX \
                      -GD \
                      -Morphology \
                      -QC_subj

        proc_structural
           This section was processed with the default parameters (T1w = ``sub-00367_ses-01_T1w.nii.gz``).

        proc_freesurfer
           This section was already processed, therefore we used the flag ``-freesurfer_dir`` to copy the already processed freesurfer files inside the ``derivatives/freesurfer/sub-00367_ses-01`` directory

        post_structural
           This section was processed with the default parameters.

        Morphology
           This section was processed with the default parameters.

        GD
           This section was processed with the default parameters.

        proc_dwi
           This section was processed with the default parameters.

           DWI shells to process: ``acq-b2000_dir-AP_dwi``

           Reverse phase encoding image = ``sub-00367_ses-01_dir-PA_dwi.nii.gz``

        SC
           This section was processed with the default parameters.

        proc_rsfmri
           White matter and CSF signal was regressed from the time-series for nuisance regression, instead of Melodic/FIX. Non linear registration (SyN) was used between the rsfMRI and T1-nativepro (default).
           No reverse phase encoding image was used.

           Main scan = ``task-rest_bold`` (default)

Cambridge Centre for Ageing and Neuroscience (Cam-CAN)
--------------------------------------------------------

Reference
   Shafto M.A., Tyler L.K., Dixon M., Taylor J.R., Rowe J.B., Cusack R., Calder A.J., Marslen-Wilson W.D., Duncan J., Dalgleish T., Henson R.N., Brayne C., Cam-CAN, & Matthews F.E. (2014).
   *The Cambridge Centre for Ageing and Neuroscience (Cam-CAN) study protocol: a cross-sectional, lifespan, multidisciplinary examination of healthy cognitive ageing*. BMC Neurology, 14 (204)

+--------+--------------------------------------------------------------------------------+
| Source | `Cam-CAN Data repository <https://www.cam-can.org/index.php?content=dataset>`_ |
+--------+--------------------------------------------------------------------------------+
| doi    | https://doi.org/10.1186/s12883-014-0204-1                                      |
+--------+--------------------------------------------------------------------------------+

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            sub-CC110045
            └── ses-pre
                ├── anat
                │   ├── mti_MTR.nii
                │   ├── sub-CC110045_T1w.json
                │   └── sub-CC110045_T1w.nii.gz
                ├── dwi
                │   ├── sub-CC110045_dwi.bval
                │   ├── sub-CC110045_dwi.bvec
                │   ├── sub-CC110045_dwi.json
                │   └── sub-CC110045_dwi.nii.gz
                └── func
                    ├── sub-CC110045_ses-pre_ses-pre_task-Rest_bold.json
                    └── sub-CC110045_ses-pre_ses-pre_task-Rest_bold.nii.gz

    .. tab:: Usage

        .. code-block:: bash
           :linenos:

            mica-pipe -bids rawdata -out derivatives -sub CC110045 -ses pre \
                      -proc_structural \
                      -proc_freesurfer \
                        -freesurfer_dir ./freesurfer_processed_CC110045 \
                      -post_structural \
                      -proc_dwi \
                      -SC \
                      -proc_rsfmri \
                        -mainScanStr task-Rest_bold \
                        -regress_WM_CSF \
                        -noFIX \
                      -GD \
                      -Morphology \
                      -MPC \
                        -microstructural_img sub-CC110045/ses-pre/mti_MTR.nii \
                      -QC_subj \

        proc_structural
           This section was processed with the default parameters (T1w = ``sub-CC110045_T1w.nii.gz``).

        proc_freesurfer
           This section was already processed, therefore we used the flag ``-freesurfer_dir`` to copy the already processed freesurfer files inside the ``derivatives/freesurfer/sub-CC110045_ses-pre`` directory

        post_structural
           This section was processed with the default parameters.

        Morphology
           This section was processed with the default parameters.

        GD
           This section was processed with the default parameters.

        proc_dwi
           This section was processed with the default parameters. No reverse phase encoding image was used.

           DWI to process: ``sub-CC110045_dwi.nii.gz``

        SC
           This section was processed with the default parameters.

        proc_rsfmri
           White matter and CSF signal was regressed from the time-series for nuisance regression, instead of Melodic/FIX. Non linear registration (SyN) was used between the rsfMRI and T1-nativepro (default).
           No reverse phase encoding image was used.
           Main scan = ``task-Rest_bold``

SUDMEX_CONN
--------------------------------------------------------
Reference
   Angeles-Valdez, D., Rasgado-Toledo, J., Issa-Garcia, V., Balducci, T., Villicana, V., Valencia, A., Gonzalez-Olvera, J. J., Reyes-Zamorano, E., Garza-Villarreal, E. A. 2021. *SUDMEX CONN: The Mexican MRI dataset of patients with cocaine use disorder*. medRxiv 2021.09.03.21263048

+--------+--------------------------------------------------------------------------------------------------------------------------------------------+
| Source | `Openneuro: SUDMEX_CONN: The Mexican dataset of cocaine use disorder patients. <https://openneuro.org/datasets/ds003346/versions/1.1.1>`_  |
+--------+--------------------------------------------------------------------------------------------------------------------------------------------+
| doi    | https://doi.org/10.1101/2021.09.03.21263048                                                                                                |
+--------+--------------------------------------------------------------------------------------------------------------------------------------------+

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            sub-010
            ├── anat
            │   ├── sub-010_T1w.json
            │   └── sub-010_T1w.nii.gz
            ├── dwi
            │   ├── sub-010_dwi.bval
            │   ├── sub-010_dwi.bvec
            │   ├── sub-010_dwi.json
            │   └── sub-010_dwi.nii.gz
            ├── fmap
            │   ├── sub-010_dir-PA_run-01_epi.json
            │   ├── sub-010_dir-PA_run-01_epi.nii.gz
            │   ├── sub-010_dir-PA_run-02_epi.json
            │   └── sub-010_dir-PA_run-02_epi.nii.gz
            └── func
                ├── sub-010_task-rest_bold.json
                └── sub-010_task-rest_bold.nii.gz

    .. tab:: Usage

        This dataset contains T1 weighted images, a HARDI-DWI single shell (b=3000, 120 directions) with a reverse phase encoding acquisition (``dir-PA_run-01_epi``),
        and a resting state fMRI of 10 minutes with its respective reverse phase encoding acquisition (``dir-PA_run-02_epi``). Both phase encoding images are in the ``fmap`` directory.

        .. code-block:: bash
           :linenos:

            subjectDir=ds003346/sub-010

            mica-pipe -bids ds003346 -out derivatives -sub 010 \
                      -proc_structural \
                      -proc_freesurfer \
                        -freesurfer_dir ./freesurfer_processed_010 \
                      -post_structural \
                      -proc_dwi \
                        -dwi_rpe ${subjectDir}/fmap/sub-110_dir-PA_run-01_epi.nii.gz \
                      -SC \
                      -proc_rsfmri \
                        -mainScanStr task-Rest_bold \
                        -fmri_rpe ${subjectDir}/fmap/sub-110_dir-PA_run-02_epi.nii.gz
                        -regress_WM_CSF \
                        - GSR \
                        -noFIX \
                      -GD \
                      -Morphology \
                      -QC_subj \

        proc_structural
           This section was processed with the default parameters (T1w = ``sub-010_T1w.nii.gz``).

        proc_freesurfer
           This section was already processed, therefore we used the flag ``-freesurfer_dir`` to copy the already processed freesurfer files inside the ``derivatives/freesurfer/sub-010`` directory

        post_structural
           This section was processed with the default parameters.

        Morphology
           This section was processed with the default parameters.

        GD
           This section was processed with the default parameters.

        proc_dwi
           This section was processed with the default parameters.

           DWI to process: ``sub-010_dwi.nii.gz``

           Reverse phase encoding image = ``fmap/sub-110_dir-PA_run-01_epi.nii.gz``

        SC
           This section was processed with the default parameters.

        proc_rsfmri
           White matter and CSF signal was regressed from the time-series for nuisance regression, instead of Melodic/FIX. Non linear registration (SyN) was used between the rsfMRI and T1-nativepro (default). Global signal regression was applied.

           Main scan = ``task-rest_bold`` (default)

           Reverse phase encoding =  ``map/sub-110_dir-PA_run-02_epi.nii.gz``

Midnight Scan Club MSC
--------------------------------------------------------

Reference
   Gordon E. M., Laumann T. O., Gilmore A. W., Newbold D. J., Greene D. J., Berg J. J., Ortega M., Hoyt-Drazen C., Gratton C., Sun H., Hampton J. M., oalson R. S., Nguyen A. L., McDermott K. B., Shimony J. S., Snyder A. Z. Schlaggar B. L., Petersen S. E., Nelson S. M., Dosenbach N. U. (2017).
   *Precision functional mapping of individual human brains.* Neuron, 95(4), 791-807.

+--------+--------------------------------------------------------------------------------------------------------------+
| Source | `Openneuro: The Midnight Scan Club (MSC) dataset <https://openneuro.org/datasets/ds000224/versions/1.0.3>`_  |
+--------+--------------------------------------------------------------------------------------------------------------+
| doi    | https://doi.org/10.1016/j.neuron.2017.07.011                                                                 |
+--------+--------------------------------------------------------------------------------------------------------------+

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            sub-MSC01
            ├── ses-func01
            │   └── func
            │       └── sub-MSC01_ses-func01_task-rest_bold.nii.gz
            ├── ses-func02
            │   └── func
            │       └── sub-MSC01_ses-func02_task-rest_bold.nii.gz
            ├── ses-func03
            │   └── func
            │   │   └── sub-MSC01_ses-func03_task-rest_bold.nii.gz
            ├── ses-func04
            │   └── func
            │       └── sub-MSC01_ses-func04_task-rest_bold.nii.gz
            ├── ses-func05
            │   └── func
            │       └── sub-MSC01_ses-func05_task-rest_bold.nii.gz
            ├── ses-func06
            │   └── func
            │       └── sub-MSC01_ses-func06_task-rest_bold.nii.gz
            ├── ses-func07
            │   └── func
            │       └── sub-MSC01_ses-func07_task-rest_bold.nii.gz
            ├── ses-func08
            │   └── func
            │       └── sub-MSC01_ses-func08_task-rest_bold.nii.gz
            ├── ses-func09
            │   └── func
            │       └── sub-MSC01_ses-func09_task-rest_bold.nii.gz
            ├── ses-func10
            │   └── func
            │       └── sub-MSC01_ses-func10_task-rest_bold.nii.gz
            └── ses-struct01
                └── anat
                    ├── sub-MSC01_ses-struct01_run-01_T1w.nii.gz
                    ├── sub-MSC01_ses-struct01_run-02_T1w.nii.gz
                    └── sub-MSC01_ses-struct01_run-01_T2w.nii.gz


    .. tab:: Usage

        This dataset contains two sessions of structural imaging, (``struct01`` and ``struct02``). However we only used the first one for processing. It does not contain DWI images, thus this step was omitted.
        Resting state functional acquisitions were acquired in different sessions (``func01`` to ``func10``), hence each session was register to the first anatomical image (``struct01``).

            .. code-block:: bash
               :linenos:

               # Pipeline for the anatomical images struct01
               mica-pipe -bids ds000224/ -out derivatives/ -sub MSC01 -ses struct01 \
                          -proc_structural -t1wStr run-01_T1w,run-02_T1w \
                          -proc_freesurfer \
                          -post_structural \
                          -proc_rsfmri \
                            -regAffine -noFIX -regress_WM_CSF \
                            -mainScanStr task-rest_bold \
                          -GD \
                          -Morphology \
                          -QC_subj

               # Pipeline for the functional images
               for N in {01..10}; do
               mica-pipe -bids ds000224/ -out derivatives/ -sub MSC01 -ses func${N} \
                          -proc_rsfmri \
                          -mainScanStr task-rest_bold \
                          -regress_WM_CSF \
                          -noFIX \
                          -sesAnat struct01 \
                          -QC_subj \
               done

            proc_structural
               Structural processing used the files ``sub-MSC01_ses-struct01_run-01_T1w.nii.gz`` and ``sub-MSC01_ses-struct01_run-01_T2w.nii.gz`` to generate the nativepro file

            proc_freesurfer
               Freesurfer used ``sub-MSC01_ses-struct01_run-01_T1w.nii.gz`` and ``sub-MSC01_ses-struct01_run-01_T2w.nii.gz`` to generate the native surfaces.

            post_structural
               This section was processed with the default parameters.

            Morphology
               This section was processed with the default parameters.

            GD
               This section was processed with the default parameters.

            proc_rsfmri
               Each functional session was processed individually and registered to the anatomical session ``struct01``. In this example we use a ``for`` loop to iterate over each session.
               White matter and CSF signal was regressed from the time-series for nuisance regression, instead of Melodic/FIX. Non linear registration (SyN) was used between the rsfMRI and T1-nativepro (default).
               No reverse phase encoding image was used.

               Main scan = ``task-rest_bold``

Auditory localization with 7T fMRI (Audiopath)
--------------------------------------------------------

Reference
   Sitek, K. R., Gulban, O. F., Calabrese, E., Johnson, G. A., Lage-Castellanos, A., Moerel, M., Ghosh S. S. & De Martino, F. (2019).
   *Mapping the human subcortical auditory system using histology, postmortem MRI and in vivo MRI at 7T*. Elife, 8, e48932. 10.7554/eLife.48932.

+--------+--------------------------------------------------------------------------------------------------------------+
| Source | `Openneuro: Auditory localization with 7T fMRI <https://openneuro.org/datasets/ds001942/versions/1.2.0>`_    |
+--------+--------------------------------------------------------------------------------------------------------------+
| doi    | https://doi.org/10.7554/eLife.48932                                                                          |
+--------+--------------------------------------------------------------------------------------------------------------+

`Audiopath <https://openneuro.org/datasets/ds001942/versions/1.2.0>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            sub-S01
            └── ses-SES01
                ├── anat
                │   ├── sub-S01_ses-SES01_acq-ShortInv_run-01_T1w.json
                │   ├── sub-S01_ses-SES01_acq-ShortInv_run-01_T1w.nii.gz
                │   ├── sub-S01_ses-SES01_run-01_T1w.json
                │   ├── sub-S01_ses-SES01_run-01_T1w.nii.gz
                │   ├── sub-S01_ses-SES01_run-02_T1w.json
                │   └── sub-S01_ses-SES01_run-02_T1w.nii.gz
                ├── dwi
                │   ├── sub-S01_ses-SES01_dir-AP_run-01_dwi.bval
                │   ├── sub-S01_ses-SES01_dir-AP_run-01_dwi.bvec
                │   ├── sub-S01_ses-SES01_dir-AP_run-01_dwi.json
                │   ├── sub-S01_ses-SES01_dir-AP_run-01_dwi.nii.gz
                │   ├── sub-S01_ses-SES01_dir-AP_run-02_dwi.bval
                │   ├── sub-S01_ses-SES01_dir-AP_run-02_dwi.bvec
                │   ├── sub-S01_ses-SES01_dir-AP_run-02_dwi.json
                │   ├── sub-S01_ses-SES01_dir-AP_run-02_dwi.nii.gz
                │   ├── sub-S01_ses-SES01_dir-AP_run-03_dwi.bval
                │   ├── sub-S01_ses-SES01_dir-AP_run-03_dwi.bvec
                │   ├── sub-S01_ses-SES01_dir-AP_run-03_dwi.json
                │   ├── sub-S01_ses-SES01_dir-AP_run-03_dwi.nii.gz
                │   ├── sub-S01_ses-SES01_dir-PA_run-01_dwi.bval
                │   ├── sub-S01_ses-SES01_dir-PA_run-01_dwi.bvec
                │   ├── sub-S01_ses-SES01_dir-PA_run-01_dwi.json
                │   ├── sub-S01_ses-SES01_dir-PA_run-01_dwi.nii.gz
                │   ├── sub-S01_ses-SES01_dir-PA_run-02_dwi.bval
                │   ├── sub-S01_ses-SES01_dir-PA_run-02_dwi.bvec
                │   ├── sub-S01_ses-SES01_dir-PA_run-02_dwi.json
                │   ├── sub-S01_ses-SES01_dir-PA_run-02_dwi.nii.gz
                │   ├── sub-S01_ses-SES01_dir-PA_run-03_dwi.bval
                │   ├── sub-S01_ses-SES01_dir-PA_run-03_dwi.bvec
                │   ├── sub-S01_ses-SES01_dir-PA_run-03_dwi.json
                │   └── sub-S01_ses-SES01_dir-PA_run-03_dwi.nii.gz
                ├── fmap
                │   ├── sub-S01_ses-SES01_acq-gre_dir-AP_run-01_epi.json
                │   ├── sub-S01_ses-SES01_acq-gre_dir-AP_run-01_epi.nii.gz
                │   ├── sub-S01_ses-SES01_acq-gre_dir-PA_run-01_epi.json
                │   └── sub-S01_ses-SES01_acq-gre_dir-PA_run-01_epi.nii.gz
                └── func
                    ├── sub-S01_ses-SES01_task-rest_bold.json
                    └── sub-S01_ses-SES01_task-rest_bold.nii.gz

    .. tab:: Usage

        This 7T in vivo dataset contains high resolution anatomical MRI, diffusion-weighted MRI and resting state fMRI. DWI images were acquired twice with the phase encoding direction of the second acquisition reversed with respect to the first.
        No real quantitative image was acquired, instead we used the T1w ShortInv for the MPC even though is not ideal.
        Additionally, it contains rsfMRI reverse phase encoding acquisitions for geometric distortion corrections.

            .. code-block:: bash
               :linenos:

                subjectDir=./ds001942/sub-S01/ses-SES01
                dwiDir=./ds001942/sub-S01/ses-SES01/dwi
                sub=sub-S02_ses-SES01

                mica-pipe -bids ds001942/ -out . -sub S01 -ses SES01 \
                          -proc_structural -t1wStr run-01_T1w,run-02_T1w \
                          -proc_freesurfer -hires \
                          -post_structural \
                          -proc_dwi \
                            -dwi_main ${dwiDir}/${sub}_dir-AP_run-01_dwi.nii.gz,${dwiDir}/${sub}_dir-AP_run-02_dwi.nii.gz,${dwiDir}/${sub}_dir-AP_run-03_dwi.nii.gz \
                            -dwi_rpe ${dwiDir}/${sub}_dir-PA_run-01_dwi.nii.gz,${dwiDir}/${sub}_dir-PA_run-02_dwi.nii.gz,${dwiDir}/${sub}_dir-PA_run-03_dwi.nii.gz \
                            -rpe_all \
                          -SC \
                          -proc_rsfmri \
                            -regAffine -noFIX -regress_WM_CSF \
                            -mainScanStr task-rest_bold \
                            -fmri_pe  ${subjectDir}/fmap/${sub}_acq-gre_dir-AP_run-01_epi.nii.gz \
                            -fmri_rpe ${subjectDir}/fmap/${sub}_acq-gre_dir-PA_run-01_epi.nii.gz \
                          -GD \
                          -Morphology \
                          -MPC \
                            -microstructural_img ${subjectDir}/anat/${sub}_acq-ShortInv_run-01_T1w.nii.gz \
                            -microstructural_reg ${subjectDir}/anat/${sub}_acq-ShortInv_run-01_T1w.nii.gz \
                          -QC_subj

            proc_structural
               Structural processing used the files ``sub-S01_ses-SES01_run-01_T1w.nii.gz`` and ``sub-S01_ses-SES01_run-02_T1w.nii.gz`` to generate the nativepro file

            proc_freesurfer
               Freesurfer used ``sub-S01_ses-SES01_run-02_T1w.nii.gz`` to generate the native surfaces. We used the ``hires`` option as is recommended for 7T isometric acquisition with high resolution or below, in this case the voxel size is 0.3x0.3x0.3 mm3.

            post_structural
               This section was processed with the default parameters.

            Morphology
               This section was processed with the default parameters.

            GD
               This section was processed with the default parameters.

            proc_dwi
               DWI directions were acquired twice with opposite phase encoding directions, thus we used the option ``rpe_all`` for the processing.

                  Main DWI shells to process: ``dir-AP_run-01``, ``dir-AP_run-02`` and ``dir-AP_run-03``

                  Reverse phase encoding DWI = ``dir-PA_run-01``, ``dir-PA_run-02`` and ``dir-PA_run-03``

            SC
               This section was processed with the default parameters.

            proc_rsfmri
               White matter and CSF signal was regressed from the time-series for nuisance regression, instead of Melodic/FIX. Only affine registration between the rsfMRI and T1-nativepro was performed.

                  Main scan = ``task-rest_bold``

                  Main phase encoding = ``fmap/${sub}_acq-gre_dir-AP_run-01_epi``

                  Reverse phase encoding = ``fmap/${sub}_acq-gre_dir-PA_run-01_epi``

            MPC
               The microstructural image as the image to register to freesurfer space were the same = ``acq-ShortInv_run-01_T1w``.
