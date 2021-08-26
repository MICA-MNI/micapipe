.. _autotract:

.. title:: Automatic Bundle Segmentation

Automatic Bundle Segmentation
================================================

Perform automatic virtual dissection of a full-brain tractogram

.. image:: autotract.png
   :align: center

https://github.com/lconcha/auto_tracto/blob/master/Readme.md

The adaptation is only partial, as in AutoPtx and XTRACT the seeding of Streamlines is performed for each bundle. Here, a full tractogram is provided, where the user had the option to seed with whatever strategy was preferred. It is assumed, however, that a full-brain seeding approach was used (white matter mask, GM/WM border, etc.)
Another difference is the use of STOP ROIs, which only make sense if seeding is performed per bundle, and not if filtering a full tractogram. STOP ROIs are used, nonetheless, but rather as termination criteria that will truncate the streamlines. Thus, STOP ROIs should be much larger than usual, to avoid the appearance of multiple short-length truncated streamlines.
The quality of the full brain tractogram will determine the quality of bundle separation. It is highly recommended to provide a tractogram with more than one million streamlines, and one that has been checked for errors. Strategies such as anatomically-contstrained tractography (ACT) and spherical deconvolution informed filtering of tractograms (SIFT), both available in MRTrix3 should aid in obtaining such high-quality tractograms.

https://github.com/lconcha/auto_tracto/blob/master/tract_definitions.md

.. tabs::
       .. tab:: Processing steps

          - A.
          - B.

       .. tab:: Usage

          The basic usage requires the same inputs as in mica-pipe:

          .. code-block:: bash
            :linenos:
            :caption: Stand alone script
            :emphasize-lines: 2

            03_auto_tracts.sh -sub 01 -ses 01 -out dataset_anonymized -bids dataset
            -all

          .. code-block:: bash
            :linenos:
            :caption: Optional argument included in ``-SC``
            :emphasize-lines: 2

            mica-pipe -sub <subject> -ses <session> -out <outputDirectory> -bids <BIDS>
            -SC -autoTract

       .. tab:: Outputs

          Directories created by this script will be ...

          .. parsed-literal::
                <subject>/<dwi>/autotract
