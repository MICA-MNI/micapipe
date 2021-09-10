.. _mic2bids:

.. title:: Mic2bids

From DICOMS to BIDS: ``mic2bids``
================================================

This section was meant to be an example of how to organize a database from sorted *DICOMS* to `BIDS <https://bids.neuroimaging.io>`_. As dicom naming and sorting can be quite unique to each imaging protocol,
the MRI sequence names and outputs must be adapted to be compatible with each dataset.

The first step starts by looking at the organization of the sorted *DICOMS* directory. You must identify the protocol names that you want to transform to *NIFTI*, and which will populate the BIDS directory.
Let's suppose that your sorted *DICOMS* directory looks like the following:

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

From this example you won't include those dicom directories with the string `AAHead_Scout_64ch-head-coil_MPR` in the BIDS directory, because those acquisitions are for calibration purposes only.


Once you've determined the directories you want to transform to *NIFTI*, you should identify their correspondence MRI sequence, for example:

+---------------------------------------------+------------------------------------------+
|                   Sequence                  |     DICOM directory name                 |
+=============================================+==========================================+
| T1 weighted (T1w)                           | - S1_08_mprage-0.8iso_ORIG               |
|                                             | - S1_20_mprage-0.8iso_ORIG               |
+---------------------------------------------+------------------------------------------+
| Quantitative T1 weighted (qT1w)             | - S1_29_mp2rage-0.8iso-sess2_INV1        |
|                                             | - S1_30_mp2rage-0.8iso-sess2_T1_Images   |
|                                             | - S1_31_mp2rage-0.8iso-sess2_UNI_Images  |
|                                             | - S1_32_mp2rage-0.8iso-sess2_INV2        |
+---------------------------------------------+------------------------------------------+
| Diffusion weighted Imaging (DWI)            | - S1_10_dwi_b2000_90                     |
|                                             | - S1_11_dwi_b700_40                      |
|                                             | - S1_12_dwi_b300_10                      |
|                                             | - S1_13_dwi_b0_5PA                       |
+---------------------------------------------+------------------------------------------+
| Blood-oxygen-level dependent imaging (bold) | - S1_16_rsfmri-3mm-bold_AP               |
|                                             | - S1_17_rsfmri-3mm_se_AP                 |
|                                             | - S1_18_rsfmri-3mm_se_PA                 |
+---------------------------------------------+------------------------------------------+

Next, you should determine the BIDS naming equivalent for each of the unique sequences on your dataset. In other words you should match the string of the DICOM directory to a BIDS string
that will identify the correspondent files, in this example it would be as follows:

+----------+---------------------------------+--------------------------+
| Sequence | DICOM string                    | BIDS string              |
+==========+=================================+==========================+
| T1w      | mprage-0.8iso_ORIG              | _T1w                     |
+----------+---------------------------------+--------------------------+
| qT1w     | mp2rage-0.8iso-sess2_INV1       | _acq-inv1_T1map          |
+----------+---------------------------------+--------------------------+
| qT1w     | mp2rage-0.8iso-sess2_T1_Images  | _acq-mp2rage_T1map       |
+----------+---------------------------------+--------------------------+
| qT1w     | mp2rage-0.8iso-sess2_UNI_Images | _acq-uni_T1map           |
+----------+---------------------------------+--------------------------+
| qT1w     | mp2rage-0.8iso-sess2_INV2       | _acq-inv2_T1map          |
+----------+---------------------------------+--------------------------+
| DWI      | dwi_b2000_90                    | _acq-b2000_dir-AP_dwi    |
+----------+---------------------------------+--------------------------+
| DWI      | dwi_b700_40                     | _acq-b700_dir-AP_dwi     |
+----------+---------------------------------+--------------------------+
| DWI      | dwi_b300_10                     | _acq-b300_dir-AP_dwi     |
+----------+---------------------------------+--------------------------+
| DWI      | dwi_b0_5PA                      | _acq-b0_dir-PA_dwi       |
+----------+---------------------------------+--------------------------+
| bold     | rsfmri-3mm-bold_AP              | _task-rest_acq-AP_bold   |
+----------+---------------------------------+--------------------------+
| bold     | rsfmri-3mm_se_AP                | _task-rest_acq-APse_bold |
+----------+---------------------------------+--------------------------+
| bold     | rsfmri-3mm_se_PA                | _task-rest_acq-PAse_bold |
+----------+---------------------------------+--------------------------+

Take a look at the `BIDS specifications <https://bids-specification.readthedocs.io/en/stable/>`_, if you have any doubts about what is the correct BIDS string, which best matches to your MRI acquisition name.

The next step is to transform your *DICOMS* directory to a *NIFTI* file with the matching BIDS naming, and under the correspondent BIDS directory.
There are different tools you can use to transform *DICOMS* to *NIFTI*, for example we used `dcm2niix` with the flag ``-b``, which creates a json card for each NIFTI file.

Let's suppose we are organizing the directories of ``sub-01_ses-01``, the main BIDS structure should be as follows:

.. parsed-literal::
    sub-01
    └── ses-01
        ├── anat
        ├── dwi
        └── func

On the next table you can find the corresponding BIDS names for each DICOM directory. This example uses ``sub-01`` as subject's identification and
``ses-01`` as the session name. As we had two T1w images, using the string ``-run`` will identify each one of them.

+----------+---------------------------------------+----------------+----------------------------------------------+
| Sequence | From DICOM directory                  | BIDS directory | To BIDS name                                 |
+==========+=======================================+================+==============================================+
| T1w      | S1_08_mprage-0.8iso_ORIG              | anat           | sub-01_ses-01_run-1_T1w.nii.gz               |
+----------+---------------------------------------+----------------+----------------------------------------------+
| T1w      | S1_20_mprage-0.8iso_ORIG              | anat           | sub-01_ses-01_run-2_T1w.nii.gz               |
+----------+---------------------------------------+----------------+----------------------------------------------+
| DWI      | S1_10_dwi_b2000_90                    | dwi            | sub-01_ses-01_acq-b2000_dir-AP_dwi.nii.gz    |
+----------+---------------------------------------+----------------+----------------------------------------------+
| DWI      | S1_11_dwi_b700_40                     | dwi            | sub-01_ses-01_acq-b700_dir-AP_dwi.nii.gz     |
+----------+---------------------------------------+----------------+----------------------------------------------+
| DWI      | S1_12_dwi_b300_10                     | dwi            | sub-01_ses-01_acq-b300_dir-AP_dwi.nii.gz     |
+----------+---------------------------------------+----------------+----------------------------------------------+
| DWI      | S1_13_dwi_b0_5PA                      | dwi            | sub-01_ses-01_acq-b0_dir-PA_dwi.nii.gz       |
+----------+---------------------------------------+----------------+----------------------------------------------+
| bold     | S1_16_rsfmri-3mm-bold_AP              | func           | sub-01_ses-01_task-rest_acq-AP_bold.nii.gz   |
+----------+---------------------------------------+----------------+----------------------------------------------+
| bold     | S1_17_rsfmri-3mm_se_AP                | func           | sub-01_ses-01_task-rest_acq-APse_bold.nii.gz |
+----------+---------------------------------------+----------------+----------------------------------------------+
| bold     | S1_18_rsfmri-3mm_se_PA                | func           | sub-01_ses-01_task-rest_acq-PAse_bold.nii.gz |
+----------+---------------------------------------+----------------+----------------------------------------------+
| qT1w     | S1_29_mp2rage-0.8iso-sess2_INV1       | anat           | sub-01_ses-01_acq-inv1_T1map.nii.gz          |
+----------+---------------------------------------+----------------+----------------------------------------------+
| qT1w     | S1_30_mp2rage-0.8iso-sess2_T1_Images  | anat           | sub-01_ses-01_acq-mp2rage_T1map.nii.gz       |
+----------+---------------------------------------+----------------+----------------------------------------------+
| qT1w     | S1_31_mp2rage-0.8iso-sess2_UNI_Images | anat           | sub-01_ses-01_acq-uni_T1map.nii.gz           |
+----------+---------------------------------------+----------------+----------------------------------------------+
| qT1w     | S1_32_mp2rage-0.8iso-sess2_INV2       | anat           | sub-01_ses-01_acq-inv2_T1map.nii.gz          |
+----------+---------------------------------------+----------------+----------------------------------------------+

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
The bash script ``mic2bids`` does all the previous steps automatically using ``dcm2niix``, however it was written to match the naming from the `MICs dataset <https://doi.org/10.1101/2021.08.04.454795>`_ DICOMS to its corresponding BIDS conform.
Feel free to adapt it to your own necessities, by modifying the DICOMS string and BIDS names in the lines 138-146.

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

**Arguments:**

		-in 	  ``path`` Input directory with the subject's DICOMS directories (FULL PATH)
		-id 	  ``string`` Subject identification for the new BIDS directory
			      *-id* CAN be different than *-in* (DICOMS directory name)
		-ses 	  ``string`` flag to specify the session name (Default is 'ses-pre')
		-bids   ``path`` Path to BIDS directory ( . or FULL PATH)
		-force 	flag that will overwrite the directory

Once you have ordered all your subjects, add the rest of the mandatory files to your BIDS directory:

  - CHANGES
  - dataset_description.json
  - participants.json
  - participants.tsv
  - README

Finally, remember to validate your dataset with the `BIDS validator <https://bids-standard.github.io/bids-validator/>`_ tool!
