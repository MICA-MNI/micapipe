.. _micapipe_anonymize:

.. title:: Micapipe Anonymize

``micapipe_anonymize``
================================================

This function anonymizes the anatomical images from the BIDS directory for data sharing. Three different methods are available for defacing/refacing. This tool uses a custom template and a set of ROIs specifically developed to identify the face and skull.

The full head template was created using the T1w images (resolution of 0.8x0.8x0.8mm) of 60 randomly selected healthy individuals from the MICA-MICs dataset (*MICs60_T1_0.8mm.nii.gz*). An inter-subject non-linear registration was performed without any mask, then the template was built using the mean of the normalized images. Three masks were generated: an ROI that covers the face, a brain mask, and a brain and neck mask. 

Unlike other algorithms, ``micapipe_anonymize`` supports different anatomical modalities, and has shown good performance across datasets and modalities.

.. image:: anonymize.png
   :align: center
   :scale: 70 %

.. tabs::


      .. tab:: Processing steps

         - Identify the T1w image for registration (run-1 in case of multiple acquisitions)
         - Affine registration from the T1w image to the face atlas
         - Anonymize images according to selected method(s):

            - Deface: Apply the affine transformation to the face mask and crop the native T1w by applying the mask
            - Reface: Apply the affine transformation to the brain-neck mask, use its inverse to crop the native T1w face, and replace with the atlas face.
            - Warpface: Create a non-linear warpfield from the T1w native face to the refaced atlas and apply the warpfield to the face.
         - An affine registration from each anatomical acquisition (e.g. quantitative T1, T2w...) to the main T1w is calculated. Then, the atlas face, the face mask and the warp field are registered to the native anatomical volume and applied.


      .. tab:: Usage

         The basic usage of ``micapipe_anonymize`` requires the same inputs as ``mica-pipe``:

         .. parsed-literal::
            $ micapipe_anonymize **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-method** **<options>**

         .. list-table::
            :widths: 10 1000
            :header-rows: 1

            * - **Options**
              - **Description**
            * - ``-sub``
              - Corresponds to subject ID. Even if your data is in BIDS, we exclude the ``sub-`` substring from the ID code (e.g. to process data for ``sub-HC001``, you would specify ``-sub HC001``). However if you forget the ``sub-`` micapipe-anonymize will manage it.
            * - ``-out``
              - Output directory path. Inside this directory the pipeline will create a new folder for each processed subject, containing the anonymized outputs.
            * - ``-bids``
              - Path to **rawdata** BIDS directory.
            * - method
              - Specify one or more methods to anonymize the structural images: ``-warpface``, ``-reface``, ``-deface``, or ``-all``

         This tool can also handle a number of **optional arguments**:

         .. list-table::
            :widths: 10 1000
            :header-rows: 1

            * - **Optional argument**
              - **Description**
            * - ``-ses`` ``string``
              - This optional flag allows the user to specify a session name (e.g. 01, 02, pre, post...). If omitted, all processing will be managed as a single session.
            * - ``-dilate`` ``num``
              - Dilation of the refaced mask (default is 6, set higher if the brain gets cropped during processing)
            * - ``-robust``
              - If reface-warped yields poor results, try running with this flag. This option to run a ROBUST registration, at the expense of more computation time
            * - ``-nocleanup``
              - Include this flag to keep (i.e. not delete) the temporary directory at script completion.
            * - ``-threads`` ``num``
              - Number of threads to use during processing (Default is 6)

         **Example:**

         Suppose you have a BIDS directory called ``dataset`` alongside an empty directory called ``dataset_anonymized``. Under the ``dataset`` directory, we see a number of anatomical images from subject 01 ``sub-01``, session 01 ``ses-01``. In this case, the participant completed three types of structural acquisitions: T1w, FLAIR and T2map.

         .. parsed-literal::
            dataset/sub-01/ses-01/anat/
            ├── sub-01_ses-01_FLAIR.json
            ├── sub-01_ses-01_FLAIR.nii.gz
            ├── sub-01_ses-01_T1w.json
            ├── sub-01_ses-01_T1w.nii.gz
            ├── sub-01_ses-01_T2map.json
            └── sub-01_ses-01_T2map.nii.gz
            dataset_anonymized

         We can run ``micapipe_anonymize`` (all methods) on these acquisitions using the following command:

         .. code-block:: bash
            :linenos:
            :caption: Example
            :emphasize-lines: 2

            micapipe_anonymize -sub 01 -ses 01 -out dataset_anonymized -bids dataset
            -all


      .. tab:: Outputs

         Directories created by this script will be contained in ``<out>/<sub>``. In our example, the output directory would then consist of ``dataset_anonymized/sub-01/``. Each output file contains a string that identifies it with the method used to anonymize it.

        .. parsed-literal::
              dataset/sub-01/ses-01/anat/
              ├── sub-01_ses-01_FLAIR.json
              ├── sub-01_ses-01_FLAIR.nii.gz
              ├── sub-01_ses-01_T1w.json
              ├── sub-01_ses-01_T1w.nii.gz
              ├── sub-01_ses-01_T2map.json
              └── sub-01_ses-01_T2map.nii.gz
              dataset_anonymized/sub-01/
              ├── sub-01_ses-01_FLAIR_defaced.nii.gz
              ├── sub-01_ses-01_FLAIR_refaced.nii.gz
              ├── sub-01_ses-01_FLAIR_warpfaced.nii.gz
              ├── sub-01_ses-01_T1w_defaced.nii.gz
              ├── sub-01_ses-01_T1w_refaced.nii.gz
              ├── sub-01_ses-01_T1w_warpfaced.nii.gz
              ├── sub-01_ses-01_T2map_defaced.nii.gz
              ├── sub-01_ses-01_T2map_refaced.nii.gz
              └── sub-01_ses-01_T2map_warpfaced.nii.gz

.. image:: sagittal.gif
   :align: center
   :scale: 50 %
