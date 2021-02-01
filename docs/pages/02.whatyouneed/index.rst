.. _what_need:

.. title:: Getting started

Usages Notes
============================================================

Data format
--------------------------------------------------------
Micapipe requires that your data be formatted in accordance with the Brain Imaging Data Structure (BIDS) format. You can find the information about the BIDS specification `here <https://bids-specification.readthedocs.io/en/stable/>`_. We strongly recommend that you validate your data structure after the conversion, notably using the `BIDS Validator <https://bids-standard.github.io/bids-validator/>`_.

.. admonition:: DICOM to BIDS: need some help? 🤓

     You can find the function used for BIDS conversion of the MICs dataset on our repository `here <https://github.com/MICA-LAB/micapipe/blob/master/functions/mic2bids>`_. As dicom naming and sorting can be quite unique to each imaging protocol, you may have to adapt this script to be compatible with your own dataset. First, your dicoms should be sorted into unique directories for each sequence before running this script. You can then modify the filename strings listed in the script (line 150-160) to correspond to the specific naming scheme of your dicom directories and their associated BIDS naming convention. Note that you will need `dcm2niix <https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage>`_ to convert sorted dicoms into NifTI files.

Set the environment
--------------------------------------------------------
If you are running a bare-metal installation of micapipe, you will need to set up your environment accordingly.

First, add micapipe to your path: ::

     export MICAPIPE=/Path/To/Cloned/Micapipe/Repo
     PATH=${PATH}:${MICAPIPE}:${MICAPIPE}/functions
     export PATH

To check if this set correctly, try displaying the help menu by running the following command from the terminal. You should see a colorful list of arguments and flags for customized runs of micapipe: ::

     mica-pipe -help

Then, you will need to also add the all dependencies (see next section for a complete list) to your PATH. For example, to add ANTs to your path: ::

     export ANTSDIR="/Path/To/ANTs"
     PATH=${PATH}:${ANTSDIR}
     export PATH

You can define distinct DIR variables for each dependency, and add them to the PATH.

.. admonition:: Why we love containers 😍

     No need to make changes to your local environment if you are going for a Docker or Singularity installation! This is all handled within the container.

Software version requirements
--------------------------------------------------------
Micapipe relies on several software dependencies. If you are opting for a bare-metal installation, you will need to set up these dependencies for all micapipe modules to run smoothly.

     - **Freesurfer**  6.0     (https://surfer.nmr.mgh.harvard.edu/)
     - **FSL**         6.0     (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
     - **AFNI**        20.2.06 (https://afni.nimh.nih.gov/download)
     - **MRtrix3**     3.0.0   (https://www.mrtrix.org)
     - **ANTs**        2.3.3   (https://github.com/ANTsX/ANTs)
     - **workbench**   1.3.2   (https://www.humanconnectome.org/software/connectome-workbench)
     - **FIX**         1.06    (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX) *optional*
     - **R**           3.6.3   (https://www.r-project.org)
     - **python**      3.7.6   (https://www.python.org/downloads/)

.. admonition:: Notes on FIX 🧐

     `FIX <https://www.sciencedirect.com/science/article/abs/pii/S1053811913011956?via%3Dihub>`_ (FMRIB’s ICA-based Xnoiseifier) is used in micapipe for removal of nuisance variable signal in resting-state fMRI data. For bare-metal installations, this portion of the functional processing will only run if FIX is found on the user's system. Note that FIX has several dependencies, specifically FSL, R and one of the following: MATLAB Runtime Component (MCR), full MATLAB or Octave. Version 1.06 of FIX relies on MATLAB 2017b/MCR v93. Additionally, it requires the following R libraries: 'kernlab','ROCR','class','party','e1071','randomForest'. 