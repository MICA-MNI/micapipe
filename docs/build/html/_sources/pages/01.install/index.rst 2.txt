.. _download_page:

.. title:: Install micapipe

Installation
========================================================

We offer three installation options for micapipe: two containerized versions (Docker and Singularity) as well as a bare-metal installation.

Containerized versions
--------------------------------------------------------
You will need to download `Docker <https://docs.docker.com/engine/install/>`_ or `Singularity <https://singularity.lbl.gov/>`_ in order to run one of the containerized versions of the pipeline.

To pull the micapipe Docker image, simply run the following command from a terminal: ::

    docker pull MICA-MNI/micapipe:latest

Similarly, to pull the micapipe Singularity image, simply run the following command from a terminal: ::

    singularity pull shub://todo/todo

Note that these containers are quite large (13+ GB), given the considerable number of dependencies for micapipe.

"Bare-metal" installation
--------------------------------------------------------
Micapipe can be directly downloaded from Github as follows: ::

    git clone https://github.com/MICA-LAB/micapipe.git

Paths to all dependencies will need to be changed manually to :ref:`set your environment <what_need>`.
