.. _micapipe_anonymize:

.. title:: Micapipe Anonymize

``micapipe_anonymize``
================================================

This function anonymizes the anatomical images from the BIDS directory for data sharing, using three different methods. It uses a custom template and a set of ROI specially developed to identify the face and skull.
The full head template was created using 60 randomly selected healthy subjects T1w images with a resolution of 0.8x0.8x0.8mm (*MICs60_T1_0.8mm.nii.gz*). An inter-subject non-linear registration was performed without any mask, then the template was built using the mean of the normalized images.
Three masks were generated: a ROI that covers the face, a brain mask, and a brain and neck mask. Derived from the inverted brain-neck mask an image of the face and skull was generated, which is used to replace the native face.
Unlike other algorithms *micapipe_anonymize* supports different anatomical modalities, and it has shown good performance.

.. image:: anonymize.png
   :align: center
   :scale: 70 %

.. tabs::
       .. tab:: Processing steps

        - Identify the T1w (run-1 in case of multiple acquisitions).
        - Affine registration from the T1w to the face atlas.
        - Deface:
           - Apply the affine transformation to the face mask and crop the native T1w with it.
        - Reface:
           - Apply the affine transformation to the brainNeck mask, use its inverse to crop the native T1w face, and replace with the atlas face.
        - Warpface:
           - Create a non-linear warpfield from the T1w native face to the refaced atlas.
           - Apply the warpfield to the face.
        - An affine registration from each anatomical acquisition to the main T1w is calculated. Then the atlas face, the face mask and the warp field are register to the native anatomical volume and applied.

       .. tab:: Usage

        The basic usage requires the same inputs as in mica-pipe:

        -sub            Corresponds to subject ID. Even if your data is in BIDS, we exclude the "sub-" substring from the ID code (e.g. sub-HC001 is entered as -sub HC001). However if you forget the ``sub-`` micapipe will manage it.
        -out            Output directory path. Following BIDS, this corresponds to the "derivatives" directory associated with your dataset. Inside this directory the pipeline will create a new folder called ``micapipe``, containing all the derivatives.
        -bids           Path to rawdata BIDS directory.
        -ses            This optional flag allows the user to specify a session name (e.g. 01, 02, pre, post...). If omitted, all processing will be managed as a single session.

        And one or more methods:

        -warpface 	     Method to anonymize the anat images by creating a new Warped-face
        -reface 	       Method to anonymize the anat images by creating a new face
        -deface 	       Method to anonymize the anat images by erasing the face
        -all             Runs all the previous algorithms

        **Example:**

        Let's suppose you have  a BIDS directory called ``dataset`` and an empty directory called ``dataset_anonymized``.
        Under the ``dataset`` directory we'll use the anatomical images from subject 01 ``sub-01`` session 01 ``ses-01``.
        In this case it contains three types of acquisitions: T1w, FLAIR and T2map.

        .. parsed-literal::
            dataset/sub-01/ses-01/anat/
            ├── sub-01_ses-01_FLAIR.json
            ├── sub-01_ses-01_FLAIR.nii.gz
            ├── sub-01_ses-01_T1w.json
            ├── sub-01_ses-01_T1w.nii.gz
            ├── sub-01_ses-01_T2map.json
            └── sub-01_ses-01_T2map.nii.gz
            dataset_anonymized

        .. code-block:: bash
          :linenos:
          :caption: In this example the script will run the three methods for anonymize
          :emphasize-lines: 2

          micapipe_anonymize -sub 01 -ses 01 -out dataset_anonymized -bids dataset
          -all

        **Optional arguments:**

        -ses 	             ``str`` Optional flag that indicates the session name (if omitted will manage as SINGLE session)
        -dilate 	         ``num`` Dilation of the refaced mask (default is 6, set higher if the brain is cropped)
        -robust 	         If reface-warped isn't great, TRY this option to run a ROBUST registration (More computation time)
        -nocleanup 	       Do NOT DELETE temporal directory at script completion.
        -threads           ``num`` Number of threads (Default is 6)

       .. tab:: Outputs

        Directories created by this script will be ``<out>/<sub>``. In our example the output directory is ``dataset_anonymized/sub-01/``.

        Each output file contains a string that identifies it with the method used to anonymized.

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
