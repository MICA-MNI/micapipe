.. _download_page:

.. title:: Install micapipe

Installation
============================================================

There are two ways to install and use micapipe: Either through virtualization/container technology, that is `Docker`_ or `Singularity`_, or in a `"Bare-metal" installation`_. Using a container-based method is highly recommended as they entail entire operating systems through kernel-level virtualization and thus include all software necessary to run micapipe, while at the same time presenting a lightweight alternative to virtual machines. Once you are ready to run micapipe, please see the `Micapipe usage overview <https://micapipe.readthedocs.io/en/latest/pages/01.whatyouneed/index.html>`_ for details.

Docker
--------------------------------------------------------

In order to run micapipe in a Docker container, Docker must be `installed <https://docs.docker.com/engine/installation/>`_ on your system. Once Docker is installed, you can get micapipe by running the following command in the terminal of your choice:

.. image:: https://img.shields.io/docker/v/micalab/micapipe?color=blue&label=docker%20version
  :target: https://hub.docker.com/r/micalab/micapipe
  :alt: Docker Image Version (latest by date)

Where ``v0.2.3`` is the latest version of ``micapipe``. We highly recommend to use the latest version because we are constantly integrating request and fixing issues. For older version check our `dockerhub <https://hub.docker.com/r/micalab/micapipe/tags>`_.

.. code-block:: bash

    docker pull micalab/micapipe:v0.2.3


After the command finished (it may take a while depending on your internet connection), you can run micapipe, in the following example only `proc_structural` is being called:

.. code-block:: bash
   :linenos:

   bids=<path to bids dataset>
   out=<path to outputs>
   fs_lic=<path to freesurfer licence>
   docker run -ti --rm \
       -v ${bids}:/bids \
       -v ${out}:/out \
       -v ${tmp}:/tmp \
       -v ${fs_lic}:/opt/licence.txt \
       micalab/micapipe:v0.2.3 \
           -bids /bids -out /out -fs_licence /opt/licence.txt -threads 10 \
           -sub $sub -ses ${ses} \
           -proc_structural


.. admonition:: Notes on ``Docker`` 🧐

      In order to Docker to run you must be sure that your bids and out directories have the correct permissions for `others`:

      bids: must have permissions to read and execute
      > `chmod go+xXr -R ${bids}`

      out:  must have permissions to read, write and execute
      > `chmod go+xXrw -R ${bids}`

      tmp:  must have permissions to read, write and execute
      > `chmod go+xXrw -R ${tmp}`

      fs_dir: must have permissions to read and execute
      > `chmod go+xXr -R ${fs_dir}`

.. warning::

   Docker and singularity have been tested only on Linux and macOS systems. Unfortunately currently we don't have support for Windows users.



Singularity
--------------------------------------------------------

For security reasons, many HPCs (e.g., TACC) do not allow Docker containers, but support `Singularity <https://github.com/singularityware/singularity>`_ containers. Depending on the Singularity v0.2.0 available to you, there are two options to get micapipe as a ``Singularity image``.

Preparing a Singularity image (Singularity version >= 2.5)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the version of Singularity on your HPC is modern enough you can create a ``Singularity image`` directly on the HCP. This is as simple as:

.. code-block:: bash

    $ singularity build /my_images/micapipe-<version>.simg docker://micalab/micapipe:<version>

Where ``<version>`` should be replaced with the desired v0.2.3 of micapipe that you want to download. For example, if you want to use ``micapipe v0.2.3``, the command would look as follows.

.. code-block:: bash

    $ singularity build /my_images/micapipe-v0.2.3.simg docker://micalab/micapipe:v0.2.3


Preparing a Singularity image (Singularity version < 2.5)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this case, start with a machine (e.g., your personal computer) with ``Docker`` installed and the use `docker2singularity <https://github.com/singularityware/docker2singularity>`_ to create a ``Singularity image``. You will need an active internet connection and some time.

.. code-block:: bash

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /absolute/path/to/output/folder:/output \
        singularityware/docker2singularity \
        micalab/micapipe:<version>

Where ``<version>`` should be replaced with the desired version of micapipe that you want to download and ``/absolute/path/to/output/folder`` with the absolute path where the created ``Singularity image`` should be stored. Sticking with the example of ``micapipe v0.2.3``, this would look as follows:

.. code-block:: bash

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v /absolute/path/to/output/folder:/output \
        singularityware/docker2singularity \
        micalab/micapipe:v0.2.3

Beware of the back slashes, expected for Windows systems. The above command would translate to Windows systems as follows:

.. code-block:: bash

    $ docker run --privileged -t --rm \
        -v /var/run/docker.sock:/var/run/docker.sock \
        -v D:\host\path\where\to\output\singularity\image:/output \
        singularityware/docker2singularity \
        micalab/micapipe:v0.2.3

You can then transfer the resulting ``Singularity image`` to the HPC, for example, using ``scp``:

.. code-block:: bash

    $ scp micalab_micapipe_v0.2.3.simg <user>@<hcpserver.edu>:/my_images

Where ``<v0.2.3>`` should be replaced with the v0.2.3 of micapipe that you used to create the ``Singularity image``, ``<user>`` with your ``user name`` on the HPC and ``<hcpserver.edu>`` with the address of the HPC.

Running a Singularity Image
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the data to be preprocessed is also on the HPC, you are ready to run micapipe:

.. code-block:: bash

    $ singularity run --writable-tmpfs --containall \
        -B path/to/your/bids_dataset:/bids \
        -B path/to/your/bids_dataset/derivatives:/out \
        -B path/to/your/working_directory:/tmp \
        -B path/to/your/freesurfer_license_file.txt:opt/licence.txt \
        /path/to/container/micapipe.simg \
        -bids /bids_dataset \
        -out /output_directory \
        -sub HC001 -ses 01 \
        -proc_structural

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

     $ micapipe -help

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

     - **Freesurfer**  7.3.2     (https://surfer.nmr.mgh.harvard.edu/)
     - **FSL**         6.0.2     (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
     - **AFNI**        20.3.03 (https://afni.nimh.nih.gov/download)
     - **MRtrix3**     3.0.1   (https://www.mrtrix.org)
     - **ANTs**        2.3.3   (https://github.com/ANTsX/ANTs)
     - **workbench**   1.3.2   (https://www.humanconnectome.org/software/connectome-workbench)
     - **FIX**         1.06    (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX) *optional*
     - **R**           3.6.3   (https://www.r-project.org)
     - **python**      3.9.16   (https://www.python.org/downloads/)

.. admonition:: Notes on ``FIX`` 🧐

     `FIX <https://www.sciencedirect.com/science/article/abs/pii/S1053811913011956?via%3Dihub>`_ (FMRIB’s ICA-based Xnoiseifier) is used in micapipe for removal of nuisance variable signal in resting-state fMRI data. For bare-metal installations, this portion of the functional processing will only run if FIX is found on the user's system. Note that FIX has several dependencies, specifically FSL, R and one of the following: MATLAB Runtime Component (MCR), full MATLAB or Octave. v0.2.0 1.06 of FIX relies on MATLAB 2017b/MCR v93. Additionally, it requires the following R libraries: 'kernlab','ROCR','class','party','e1071','randomForest'.
