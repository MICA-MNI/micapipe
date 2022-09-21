.. _download_page:

.. title:: Install micapipe

Installation
============================================================

There are two distinct ways to install and use micapipe: Either through virtualization/container technology, that is `Docker`_ or `Singularity`_, or in a `"Bare-metal" installation`_. Using a container-based method is highly recommended as they entail entire operating systems through kernel-level virtualization and thus include all software necessary to run micapipe, while at the same time presenting a lightweight alternative to virtual machines. Once you are ready to run micapipe, please see the `Micapipe usage overview <https://micapipe.readthedocs.io/en/latest/pages/01.execution/index.html>`_ for details.

Docker
--------------------------------------------------------

In order to run micapipe in a Docker container, Docker must be `installed <https://docs.docker.com/engine/installation/>`_ on your system. Once Docker is installed, you can get micapipe by running the following command in the terminal of your choice:

.. code-block:: bash

    docker pull micalab/micapipe:version

Where ``version`` is the specific version of ``micapipe`` you would like to use. For example, if you want to employ the latest/most up to date ``version`` you can either run 

.. code-block:: bash

    docker pull micalab/micapipe:latest

or the same command without the ``:latest`` tag, as Docker searches for the ``latest`` tag by default. However, as the ``latest`` version is subject to changes and not necessarily in sync with the most recent ``numbered version``, it is recommend to utilize the latter to ensure reproducibility. For example, if you want to employ ``micapipe v0.1.2`` the command would look as follows:

.. code-block:: bash

    docker pull micalab/micapipe:v0.1.2

After the command finished (it may take a while depending on your internet connection), you can run micapipe like this:

.. code-block:: bash

    $ docker run -ti --rm \
        -v path/to/your/bids_dataset:/bids_dataset:ro \
        -v path/to/your/bids_dataset/derivatives:/output_directory \
        -v path/to/your/working_directory:/tmp \
        -v path/to/your/freesurfer_license_file.txt:/opt/freesurfer-6.0.0/freesurfer_license_file.txt \
        micalab/micapipe:latest \
        -bids /bids_dataset \
        -out /output_directory \
        -all -ses 01 \
        -threads 10 -tracts 10M

In the above example, the flag ``-all`` would run all micapipe modules for the specified dataset. Please have a look at the examples under `Micapipe usage overview <https://micapipe.readthedocs.io/en/latest/pages/01.execution/index.html>`_ to get more information and familiarize yourself with micapipe's functionality!


Singularity
--------------------------------------------------------

For security reasons, many HPCs (e.g., TACC) do not allow Docker containers, but support `Singularity <https://github.com/singularityware/singularity>`_ containers. Depending on the Singularity version available to you, there are two options to get micapipe as a ``Singularity image``.

Preparing a Singularity image (Singularity version >= 2.5)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the version of Singularity on your HPC is modern enough you can create a ``Singularity image`` directly on the HCP. This is as simple as: 

.. code-block:: bash

    $ singularity build /my_images/micapipe-<version>.simg docker://micalab/micapipe:<version>

Where ``<version>`` should be replaced with the desired version of micapipe that you want to download. For example, if you want to use ``micapipe v0.1.2``, the command would look as follows.

.. code-block:: bash

    $ singularity build /my_images/micapipe-v0.1.2.simg docker://micalab/micapipe:v0.1.2


Preparing a Singularity image (Singularity version < 2.5)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this case, start with a machine (e.g., your personal computer) with ``Docker`` installed and the use `docker2singularity <https://github.com/singularityware/docker2singularity>`_ to create a ``Singularity image``. You will need an active internet connection and some time. 

.. code-block:: bash

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /absolute/path/to/output/folder:/output \
        singularityware/docker2singularity \
        micalab/micapipe:<version>

Where ``<version>`` should be replaced with the desired version of micapipe that you want to download and ``/absolute/path/to/output/folder`` with the absolute path where the created ``Singularity image`` should be stored. Sticking with the example of ``micapipe v0.1.2``, this would look as follows:

.. code-block:: bash

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /absolute/path/to/output/folder:/output \
        singularityware/docker2singularity \
        micalab/micapipe:v0.1.2

Beware of the back slashes, expected for Windows systems. The above command would translate to Windows systems as follows:

.. code-block:: bash

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v D:\host\path\where\to\output\singularity\image:/output \
        singularityware/docker2singularity \
        micalab/micapipe:<version>

You can then transfer the resulting ``Singularity image`` to the HPC, for example, using ``scp``:

.. code-block:: bash

    $ scp micalab_micapipe<version>.simg <user>@<hcpserver.edu>:/my_images

Where ``<version>`` should be replaced with the version of micapipe that you used to create the ``Singularity image``, ``<user>`` with your ``user name`` on the HPC and ``<hcpserver.edu>`` with the address of the HPC.  

Running a Singularity Image
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the data to be preprocessed is also on the HPC, you are ready to run micapipe:

.. code-block:: bash

    $ singularity run --cleanenv \
        -B path/to/your/bids_dataset:/bids_dataset:ro \
        -B path/to/your/bids_dataset/derivatives:/output_directory \
        -B path/to/your/working_directory:/tmp \
        -B path/to/your/freesurfer_license_file.txt:/opt/freesurfer-6.0.0/license.txt \
        /path/to/container/micapipe.simg \
        -bids /bids_dataset \
        -out /output_directory \
        -all -ses 01 \
        -threads 10 -tracts 10M

.. admonition:: Some things to consider with Singularity 🙆‍♀️

    - Make sure to check the name of the created ``Singularity image``, as that might diverge based on the method you used. Here and going forward it is assumed that you used ``Singularity >= 2.5`` and thus ``micapipe-<version>.simg`` instead of ``micalab_micapipe<version>.simg``.
    - Singularity by default `exposes all environment variables from the host inside the container <https://github.com/singularityware/singularity/issues/445>`_. Because of this your host libraries could be accidentally used instead of the ones inside the container. To avoid such situation we recommend using the ``--cleanenv`` singularity flag.
    - Depending on how Singularity is configured on your cluster, it might or might not automatically ``bind`` (``mount`` or ``expose``) ``host folders`` to the container. If this is not done automatically you will need to ``bind`` the necessary folders using the ``-B <host_folder>:<container_folder>`` Singularity argument.


"Bare-metal" installation
--------------------------------------------------------

.. warning::

   This method is not recommended! Using a `Docker`_ or a `Singularity`_ might avoid a lot of headaches...

For this route, you will need to make sure all of micapipe's `External Dependencies`_ are installed. These tools must be installed and their binaries available in the system's ``$PATH``. A relatively interpretable description of how your environment can be set-up is found in the `Dockerfile <https://github.com/MICA-MNI/micapipe/blob/master/Dockerfile>`_ as well as in the `init.sh <https://github.com/MICA-MNI/micapipe/blob/master/functions/init.sh>`_ script provided in the micapipe repository. 

Micapipe can be directly downloaded from Github as follows:

.. code-block:: bash

    $ git clone https://github.com/MICA-LAB/micapipe.git

Paths to all dependencies will need to be changed manually to `Set the environment`_.

Set the environment
^^^^^^^^^^^^^^^^^^^
If you are running a bare-metal installation of micapipe, you will need to set up your environment accordingly.

First, add micapipe to your ``$PATH``:

.. code-block:: bash

     $ export MICAPIPE=/Path/To/Cloned/Micapipe/Repo
     $ PATH=${PATH}:${MICAPIPE}:${MICAPIPE}/functions
     $ export PATH

To check if this set correctly, try displaying the help menu by running the following command from the terminal. You should see a colorful list of arguments and flags for customized runs of micapipe:

.. code-block:: bash

     $ mica-pipe -help

Then, you will need to also add the all dependencies (see next section for a complete list) to your ``$PATH``. For example, to add ANTs to your ``$PATH``:

.. code-block:: bash

     $ export ANTSDIR="/Path/To/ANTs"
     $ PATH=${PATH}:${ANTSDIR}
     $ export PATH

You can define distinct DIR variables for each dependency, and add them to the ``$PATH``.

.. admonition:: Why we love containers 😍

     No need to make changes to your local environment if you are going for a Docker or Singularity installation! This is all handled within the container.


External Dependencies
^^^^^^^^^^^^^^^^^^^^^
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

.. admonition:: Notes on ``FIX`` 🧐

     `FIX <https://www.sciencedirect.com/science/article/abs/pii/S1053811913011956?via%3Dihub>`_ (FMRIB’s ICA-based Xnoiseifier) is used in micapipe for removal of nuisance variable signal in resting-state fMRI data. For bare-metal installations, this portion of the functional processing will only run if FIX is found on the user's system. Note that FIX has several dependencies, specifically FSL, R and one of the following: MATLAB Runtime Component (MCR), full MATLAB or Octave. Version 1.06 of FIX relies on MATLAB 2017b/MCR v93. Additionally, it requires the following R libraries: 'kernlab','ROCR','class','party','e1071','randomForest'.

