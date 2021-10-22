.. _micapipe_cleanup:

.. title:: Micapipe Cleanup

``micapipe_cleanup``
================================================

This script was design to erase all outputs, including .log files of a previously run module.
You can use this function in case you encounter an error while processing a module, or if you simply want to rerun a module that you have already started or that is completed. Indeed, if the output files of a module already exist, you will not be able to rerun this module before deleting its outputs. 

.. WARNING:: 
   
   If you want to erase and rerun a core processing module such as ``post_structural`` or ``proc_structural`` you may need to rerun all of its dependent modules, to ensure consistency of interdependent modules!

The basic usage of ``micapipe_cleanup`` requires the same inputs as in mica-pipe.

Terminal:

.. parsed-literal::
   $ micapipe_cleanup **-sub** <subject_id> **-out** <outputDirectory> **-bids** <BIDS-directory> **-ses** <session-name> **-<module-to-erase>**

Docker:

.. parsed-literal::
   $ docker micapipe_cleanup

The module(s) you want to erase should be specified with a flag(s) (e.g. ``-proc_dwi``). Possible optinal arguments are listed below:

.. list-table::
    :widths: 10 1000
    :header-rows: 1
    :class: tight-table

    * - **Optional arguments**
      - **Description**
    * - ``-ses``
      - Optional flag that indicates the session name (if omitted, will be managed as a SINGLE session)
    * - ``-proc_structural``
      - Deletes volumetric processing derivatives
    * - ``-proc_freesurfer``
      - Deletes Freesurfer recon-all processing outputs
    * - ``-post_structural``
      - Deletes post structural volumetric and surface-based processing derivatives
    * - ``-proc_dwi``
      - Deletes diffusion-weighted image processing derivatives
    * - ``-SC``
      - Deletes tractography files and structural connectomes
    * - ``-proc_rsfmri``
      - Deletes resting-state funtional MRI processing derivatives
    * - ``-MPC``
      - Deletes microstructural profile and covariance analysis derivatives
    * - ``-GD``
      - Deletes geodesic distance analysis derivatives