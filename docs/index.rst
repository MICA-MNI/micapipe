.. MICAPIPE documentation master file, created by
   sphinx-quickstart on Wed Jul 15 16:09:38 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**micapipe**
============================================================

.. title:: micapipe

.. raw:: html

   <style type="text/css">
      hr {
      width: 100%;
      height: 1px;
      background-color: #261F4A;
      margin-top: 24px;
      }
   </style>

.. .. image:: https://readthedocs.org/projects/micapipe/badge/?version=latest
  :target: https://micapipe.readthedocs.io/en/latest/?badge=latest
  :alt: Documentation Status

**Welcome to micapipe's documentation!**
============================================================

*Micapipe* is a processing pipeline providing a robust framework to analyze multimodal MRI data. Our pipeline integrates processing streamlines for *T1-weighted*, *microstructure-sensitive*, *diffusion-weighted*, and *resting-state functional imaging* to facilitate the development of multiscale models of neural organization. For this purpose, we leverage several neuroimaging software packages to bring BIDS-formatted raw MRI data to fully-processed surface-based connectivity matrices.

Micapipe is developed by MICA-lab (https://mica-mni.github.io) and collaborators at the McConnell Brain Imaging Center of the Montreal Neurological Institute. 

.. image:: ./figures/micapipe.png
   :scale: 50 %
   :alt: alternate text
   :align: center

.. raw:: html

   <br>

Reproducibility
--------------------------------------------------------
To encourage reproducibility and robustness of investigations, we provide a fully containerized version of Micapipe in the form of a Docker container. Step-by-step tutorials are provided so users can run micapipe through this container. For flexibility, we also provide guidance on bare metal installations. 

.. raw:: html

   <br>

Datasets
--------------------------------------------------------
Micapipe has been tested on several openly available datasets, such as CamCAN, [...]. The pipeline was optimized for the MICs dataset (DOI: ), which includes all data modalities that can be processed in micapipe.

.. raw:: html

   <br>

Development and getting involved
--------------------------------------------------------
GitHub issues, pull requests, mailing list

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting started

   pages/01.install/index
   pages/02.whatyouneed/index
   pages/03.execution/index
   pages/09.whatsnew/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Processing steps

   pages/04.volumetric/index
   pages/05.freesurfer/index
   pages/06.dwi/index
   pages/07.restingstate/index
   pages/08.postmpc/index
   pages/17.geodesic/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Additional tools

   pages/14.micapipe_anonymize/index
   pages/15.micapipe_cleanup/index
   pages/16.mic2bids/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: References & Acknowledgements

   pages/13.writeitdown/index
   pages/10.citingmicapipe/index
   pages/11.references/index
   pages/12.acknowledge/index


--------------------------------------------------------

.. raw:: html

   <br>

Core development team 🧠
--------------------------------------------------------

- **Raúl Rodríguez-Cruces**, *MICA Lab - Montreal Neurological Institute*
- **Jessica Royer**, *MICA Lab - Montreal Neurological Institute*
- **Sara Larivière**, *MICA Lab - Montreal Neurological Institute*
- **Peer Herholz**, McConnell Brain Imaging Centre - Montreal Neurological Institute*
- **Bo-yong Park**, *MICA Lab - Montreal Neurological Institute*
- **Reinder Vos de Wael**, *MICA Lab - Montreal Neurological Institute*
- **Casey Paquola**, *MICA Lab - Montreal Neurological Institute*
- **Oualid Benkarim**, *MICA Lab - Montreal Neurological Institute*
- **Boris Bernhardt**, *MICA Lab - Montreal Neurological Institute*