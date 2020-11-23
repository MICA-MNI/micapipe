.. _what_need:

.. title:: How to get started?

Usages Notes
======================================

Data format
-----------
To use the micapipe, we recommand to used Brain Imaging Data Structure (BIDS) format. Before the execution of the micapipe,
you will need to convert your dataset in valid BIDS format. You can find the information about the BIDS specification `here <https://bids-specification.readthedocs.io/en/stable/>`_.
We recommend you to validate your dataset after the conversion online, `BIDS Validator <https://bids-standard.github.io/bids-validator/>`_.

.. admonition:: Essentials to have ☝🏼

     You will need to download `dcm2niix <https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage>`_ to complete the conversion of your dataset.

Set the environment
--------------------


Software versions requirements
-------------------------------
micapipe written using...

   micapipe requires software toolx before to start:

    - **Freesurfer**  6.0     (https://surfer.nmr.mgh.harvard.edu/)
    - **FSL**         6.0     (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
    - **AFNI**        20.2.06 (https://afni.nimh.nih.gov/download)
    - **MRtrix3**     3.0.0   (https://www.mrtrix.org)
    - **ANTs**        2.3.3   (https://github.com/ANTsX/ANTs)
    - **workbench**   1.3.2
    - **ENIGMA toolbox** (https://github.com/MICA-MNI/ENIGMA)

   Optional software:
    - **FIX** (FMRIB's ICA-based Xnoiseifier) v1.06
