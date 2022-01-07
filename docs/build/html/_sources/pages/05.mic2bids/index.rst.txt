.. _mic2bids:

.. title:: Mic2bids

From DICOMS to BIDS: ``mic2bids``
================================================

This section describes how to organize a database from sorted *DICOMs* files to `BIDS format <https://bids.neuroimaging.io>`_. As the naming and sorting of DICOMs can be quite unique to each imaging protocol,
 MRI sequence names and outputs must be adapted to be compatible with each dataset.

First, you will have to look at the organization of the sorted DICOM directory. This directory should contain subdirectories in which the DICOM images are sorted according to each sequence completed by the participant. You must identify the sequences in the protocol that you want to transform to NIfTI: These sequences will then populate the BIDS directory.

Let's suppose that your sorted DICOM directory looks like the following:

.. parsed-literal::
    S1_02_AAHead_Scout_64ch-head-coil_MPR_sag
    S1_03_AAHead_Scout_64ch-head-coil_MPR_cor
    S1_04_AAHead_Scout_64ch-head-coil_MPR_tra
    S1_08_mprage-0.8iso_ORIG
    S1_10_dwi_b2000_90
    S1_11_dwi_b700_40
    S1_12_dwi_b300_10
    S1_13_dwi_b0_5PA
    S1_16_rsfmri-3mm-bold_AP
    S1_17_rsfmri-3mm_se_AP
    S1_18_rsfmri-3mm_se_PA
    S1_20_mprage-0.8iso_ORIG
    S1_29_mp2rage-0.8iso-sess2_INV1
    S1_30_mp2rage-0.8iso-sess2_T1_Images
    S1_31_mp2rage-0.8iso-sess2_UNI_Images
    S1_32_mp2rage-0.8iso-sess2_INV2

In this example, we don't want images in the DICOM directories containing the string `AAHead_Scout_64ch-head-coil_MPR` to be included in the BIDS directory, because those acquisitions are for calibration purposes only.

Once you've determined which directories you want to transform to NIfTI and organize into BIDS, you should identify their corresponding MRI sequence. You should determine the BIDS naming equivalent of each unique sequence in your dataset. In other words, you should match the string of each DICOM directory to a BIDS string that will identify the corresponding sequence.

.. admonition:: BIDS naming 🗄

    Take a look at the `BIDS specifications <https://bids-specification.readthedocs.io/en/stable/>`_ if you have any doubts about the BIDS string that best matches your MRI acquisition.

The next step is to transform your DICOM directory to a NIfTI file with the matching BIDS naming, and place it under its corresponding BIDS directory. Different online tools can transform DICOMs to NIfTI. For example, we use `dcm2niix <https://github.com/rordenlab/dcm2niix>`_ with the flag ``-b``, which creates a json card for each NIfTI file.

We can refer to to the following table to see the DICOM and BIDS naming of each sequence, and which BIDS subdirectory it will store each file. This example uses ``sub-01`` as subject's identification and
``ses-01`` as the session name. As we acquired two distinct T1w images, using the string ``-run`` will differentiate them.

.. list-table:: DICOM to BIDS
    :widths: 30 20 10 5 20
    :header-rows: 1
    :class: tight-table

    * - **Sequence**
      - **DICOM directory name**
      - **BIDS string**
      - **BIDS directory**
      - **Example subject naming**
    * - T1 weighted (T1w)
      - S1_08_mprage-0.8iso_ORIG
      - ``_T1w``
      - /anat
      - ``sub-01_ses-01_run-1_T1w.nii.gz``
    * - T1 weighted (T1w)
      - S1_20_mprage-0.8iso_ORIG
      - ``_T1w``
      - /anat
      - ``sub-01_ses-01_run-2_T1w.nii.gz``      
    * - Quantitative T1 (qT1), first inversion time
      - S1_29_mp2rage-0.8iso-sess2_INV1
      - ``_acq-inv1_T1map``
      - /anat
      - ``sub-01_ses-01_acq-inv1_T1map.nii.gz``
    * - Quantitative T1 (qT1), second inversion time
      - S1_32_mp2rage-0.8iso-sess2_INV2
      - ``_acq-inv2_T1map``
      - /anat
      - ``sub-01_ses-01_acq-inv2_T1map.nii.gz``
    * - Quantitative T1 (qT1), T1 map
      - S1_30_mp2rage-0.8iso-sess2_T1_Images
      - ``_acq-mp2rage_T1map``
      - /anat
      - ``sub-01_ses-01_acq-mp2rage_T1map.nii.gz``
    * - Quantitative T1 (qT1), synthetic T1
      - S1_31_mp2rage-0.8iso-sess2_UNI_Images
      - ``_acq-uni_T1map``
      - /anat
      - ``sub-01_ses-01_acq-uni_T1map.nii.gz``
    * - Diffusion weighted Imaging (DWI), b=2000
      - S1_10_dwi_b2000_90
      - ``_acq-b2000_dir-AP_dwi``
      - /dwi
      - ``sub-01_ses-01_acq-b2000_dir-AP_dwi.nii.gz``
    * - Diffusion weighted Imaging (DWI), b=700
      - S1_11_dwi_b700_40
      - ``_acq-b700_dir-AP_dwi``
      - /dwi
      - ``sub-01_ses-01_acq-b700_dir-AP_dwi.nii.gz``
    * - Diffusion weighted Imaging (DWI), b=300
      - S1_12_dwi_b300_10
      - ``_acq-b300_dir-AP_dwi``
      - /dwi
      - ``sub-01_ses-01_acq-b300_dir-AP_dwi.nii.gz``
    * - Diffusion weighted Imaging (DWI), b=0
      - S1_13_dwi_b0_5PA
      - ``_acq-b0_dir-PA_dwi``
      - /dwi
      - ``sub-01_ses-01_acq-b0_dir-PA_dwi.nii.gz``
    * - Blood-oxygen-level dependent imaging (BOLD)
      - S1_16_rsfmri-3mm-bold_AP
      - ``_task-rest_acq-AP_bold``
      - /func
      - ``sub-01_ses-01_task-rest_acq-AP_bold.nii.gz``
    * - Fieldmap for BOLD (AP)
      - S1_17_rsfmri-3mm_se_AP
      - ``_task-rest_acq-APse_bold``
      - /func
      - ``sub-01_ses-01_task-rest_acq-APse_bold.nii.gz``
    * - Fieldmap for BOLD (PA)
      - S1_18_rsfmri-3mm_se_PA
      - ``_task-rest_acq-PAse_bold``
      - /func
      - ``sub-01_ses-01_task-rest_acq-PAse_bold.nii.gz``


Let's suppose we are organizing the directories of ``sub-01_ses-01``, the main BIDS structure should be as follows:

.. parsed-literal::
    sub-01
    └── ses-01
        ├── anat
        ├── dwi
        └── func

The final BIDS structure should look like:

.. parsed-literal::
    sub-HC001
    └── ses-01
        ├── anat
        │   ├── sub-01_ses-01_acq-inv1_T1map.json
        │   ├── sub-01_ses-01_acq-inv1_T1map.nii.gz
        │   ├── sub-01_ses-01_acq-inv2_T1map.json
        │   ├── sub-01_ses-01_acq-inv2_T1map.nii.gz
        │   ├── sub-01_ses-01_acq-mp2rage_T1map.json
        │   ├── sub-01_ses-01_acq-mp2rage_T1map.nii.gz
        │   ├── sub-01_ses-01_acq-uni_T1map.json
        │   ├── sub-01_ses-01_acq-uni_T1map.nii.gz
        │   ├── sub-01_ses-01_T1w.json
        │   └── sub-01_ses-01_T1w.nii.gz
        ├── dwi
        │   ├── sub-01_ses-01_acq-b2000_dir-AP_dwi.bval
        │   ├── sub-01_ses-01_acq-b2000_dir-AP_dwi.bvec
        │   ├── sub-01_ses-01_acq-b2000_dir-AP_dwi.json
        │   ├── sub-01_ses-01_acq-b2000_dir-AP_dwi.nii.gz
        │   ├── sub-01_ses-01_acq-b300_dir-AP_dwi.bval
        │   ├── sub-01_ses-01_acq-b300_dir-AP_dwi.bvec
        │   ├── sub-01_ses-01_acq-b300_dir-AP_dwi.json
        │   ├── sub-01_ses-01_acq-b300_dir-AP_dwi.nii.gz
        │   ├── sub-01_ses-01_acq-b700_dir-AP_dwi.bval
        │   ├── sub-01_ses-01_acq-b700_dir-AP_dwi.bvec
        │   ├── sub-01_ses-01_acq-b700_dir-AP_dwi.json
        │   ├── sub-01_ses-01_acq-b700_dir-AP_dwi.nii.gz
        │   ├── sub-01_ses-01_dir-PA_dwi.bval
        │   ├── sub-01_ses-01_dir-PA_dwi.bvec
        │   ├── sub-01_ses-01_dir-PA_dwi.json
        │   └── sub-01_ses-01_dir-PA_dwi.nii.gz
        └── func
            ├── sub-01_ses-01_task-rest_acq-AP_bold.json
            ├── sub-01_ses-01_task-rest_acq-AP_bold.nii.gz
            ├── sub-01_ses-01_task-rest_acq-APse_bold.json
            ├── sub-01_ses-01_task-rest_acq-APse_bold.nii.gz
            ├── sub-01_ses-01_task-rest_acq-PAse_bold.json
            └── sub-01_ses-01_task-rest_acq-PAse_bold.nii.gz

It is the same procedure for each subject in your dataset.

The bash script ``mic2bids`` completes all the previous steps automatically using ``dcm2niix``. However it was written to match the name of the DICOMs in the `MICA-MICs dataset <https://doi.org/10.1101/2021.08.04.454795>`_ with their corresponding BIDS name. Feel free to adapt it to your own needs, by modifying the DICOM strings and BIDS names in the lines 138-146 of ``mic2bids``, or use it as a guideline to build your own script.

.. parsed-literal::
    # -----------------------------------------------------------------------------------------------
    # CHANGE THIS regex (regular expressions) ACCORDING TO YOUR DICOMS NAMING
    orig=("*mprage-0.8iso*" "*_INV1" "*_INV2" "*_T1_Images" "*UNI_Images" "*FLAIR*" "*rsfmri-3mm-bold_AP" "*rsfmri-3mm_se_AP" "*rsfmri-3mm_se_PA")
    origDWI=("*dwi_b700_40" "*dwi_b300_10" "*dwi_b0_5PA" "*_dwi_b2000*")

    # New BIDS-naming, follow the BIDS specification:
    # https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/01-magnetic-resonance-imaging-data.html
    bids=(T1w acq-inv1_T1map acq-inv2_T1map acq-mp2rage_T1map acq-uni_T1map FLAIR task-rest_acq-AP_bold task-rest_acq-APse_bold task-rest_acq-PAse_bold)
    bidsDWI=(acq-b700-NUM_dir-AP_dwi acq-b300-NUM_dir-AP_dwi dir-PA_dwi acq-b2000-NUM_dir-AP_dwi)

Once you've done all the necessary modifications you can run the scrip:

.. parsed-literal::
   mic2bids  -in <DICOMS_directory> -bids <BIDS directory path> -id <subject> -ses <session>

.. list-table::
    :widths: 75 750
    :header-rows: 1
    :class: tight-table

    * - **Argument**
      - **Description**
    * - ``-in`` ``<path>``
      - Input directory with the subject's DICOMS directories (FULL PATH)
    * - ``-id`` ``string``
      - Subject identification for the new BIDS directory ``-id`` CAN be different than ``-in`` (DICOMS directory name)
    * - ``-ses`` ``string``
      - Flag to specify the session name (Default is ``ses-pre``)
    * - ``-bids`` ``<path>``
      - Path to BIDS directory ( . or FULL PATH)
    * - ``-force``
      - This flag will overwrite the existing BIDS directory for that subject

Once you have ordered all your subjects, add the rest of the mandatory files to your BIDS directory:

  - CHANGES
  - dataset_description.json
  - participants.json
  - participants.tsv
  - README

Finally, remember to validate your dataset with the `BIDS validator <https://bids-standard.github.io/bids-validator/>`_ tool prior to running micapipe!
