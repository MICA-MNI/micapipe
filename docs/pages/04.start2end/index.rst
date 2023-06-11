.. _startend:

.. title:: Start to finish with `micapipe`

Processing step by step: start to finish with ``micapipe``
============================================================

In this section, you will find two examples with the necessary steps of processing a dataset with micapipe.

1. Download an open access dataset
--------------------------------------------------------

The first step is to identify a dataset and to ensure that you have the computational resources, including storage and computational processing power. In this this example, we will focus on two datasets: Human Connectome Project (HCP) and (MICA-MICs). Other widely used repositories for BIDS-compliant datasets may be found on e.g., OpenNeuro (https://openneuro.org) or the Canadian Open Neuroscience Platform (CONP; http://portal.conp.ca).

.. tabs::

    .. tab:: HCP

        .. figure:: hcp_logo.png
            :align: left
            :scale: 40 %

        To download the dataset, it is necessary to create an account first and accept the open access terms. Additionally, HCP uses a third-party software called Aspera Connect to boost data transfer speed. Further information about this software and its installation can be found on the HCP website. Once all the requirements are fulfilled, login to  https://db.humanconnectome.org, and select the data you would like to download, in this case the WU-Minn HCP Retest Data with the processing filter of “Unprocessed”.
        We will use micapipe here to process the structural, resting state fMRI run-1 and Diffusion data (dir97_dir-RL).
        This dataset requires about 350 GB of storage.

    .. tab:: MICs

        .. figure:: ../04.databases/mics_logo.png
            :align: left
            :scale: 17 %

        There are two options available to download the MICs dataset, a direct download or using Datalad.
        This dataset requires approximately 70 GB of storage. The dataset and additional documentation can be found on the CONP website:
        https://portal.conp.ca/dataset?id=projects/mica-mics.

2. Converting to BIDS
--------------------------------------------------------

.. tabs::

    .. tab:: HCP

        The HCP dataset was created before the rise of BIDS. Here, we provide a custom-build script that will transform an HCP directory into BIDS format, using the metadata provided in the HCP S1200 release reference manual.

        1. Clone the micapipe-supplementary github directory:

        .. code-block:: bash
           :linenos:

           git clone https://github.com/MICA-MNI/micapipe-supplementary.git

        2. Change into the project directory:

        .. code-block:: bash
           :linenos:

           cd micapipe-supplementary /functions

        3. Run the code specifying the full path to the HCP data and the output HCP directory in BIDS.

        .. code-block:: bash
           :linenos:

           ./hcp2bids -in <full_path_to>/HCP_data -out <full_path_to>/HCP_bids

    .. tab:: MICs

        This dataset is BIDS-compliant, and no further conversion is needed.

3. Validating BIDS
--------------------------------------------------------

At this point, MICs and HCP are BIDS conform. However, any new dataset that has been recently acquired and that you wish to make BIDS-compliant (see specifications at http://bids-specification.readthedocs.io), should be validated with tools provided by BIDS, such as the BIDS-validator (https://bids-standard.github.io/bids-validator/). Another example can be found on the tutorial `From Dicoms to BIDS: mic2bids <../05.mic2bids/index.html>`_

4. Running micapipe
--------------------------------------------------------

Once micapipe has been installed (see Installation), one can run the pipeline. From the main directory of the dataset the command would be:

.. tabs::

    .. tab:: HCP

        Running HCP for subject 250932

        .. code-block:: bash
           :linenos:

           micapipe -bids HCP_bids -out derivatives -sub 250932 \
                  -proc_structural \
                  -proc_freesurfer \
                  -post_structural \
                  -proc_dwi -dwi_acq dir97 \
        	       -dwi_main sub-250932dwi/sub-250932_acq-dir97_dir-LR_dwi.nii.gz \
                  -dwi_rpe sub-250932/dwi/sub-250932_acq-dir97_dir-RL_sbref.nii.gz  \
                  -SC -tracts 20M \
                  -proc_rsfmri \
                  -mainScanStr task-rest_dir-LR_run-2_bold \
                  -fmri_rpe sub-250932/func/sub-250932_task-rest_dir-RL_run-1_bold.nii.gz\
                  -regress_WM_CSF -noFIX -regAffine \
                  -GD \
                  -Morphology \
                  -QC_subj

    .. tab:: MICs

        Running MICs for subject HC001 session 01:

        .. code-block:: bash
           :linenos:

           micapipe -bids rawdata -out derivatives -sub HC001 -ses 01 \
                  -proc_structural \
                  -proc_freesurfer \
                  -post_structural \
                  -proc_dwi \
                  -SC -tracts 20M \
                  -proc_rsfmri \
                  -GD \
                  -Morphology \
                  -MPC \
                  -QC_subj

5. Visualize the QC report
--------------------------------------------------------

The individual QC tool is a html report with detailed information of each processing step, which can be used for rapid visualization of processing status, core registrations, and data matrices by parcellation scheme and module. The files can be found under each subject’s directory QC and opened with any browser:

.. code-block:: bash

   <outputDirectory>/micapipe/<sub>/QC/<sub>_micapipe_qc.html

The group level QC generates a report with all completed and processed modules by subject. The report consists of a color-coded table with rows as subjects and columns as the pipeline modules. The file is located under the micapipe directory as shown below:

.. code-block:: bash

   <outputDirectory>/micapipe/micapipe_progress.html
