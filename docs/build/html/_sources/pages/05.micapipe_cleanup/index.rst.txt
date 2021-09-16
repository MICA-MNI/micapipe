.. _micapipe_cleanup:

.. title:: Micapipe Cleanup

``micapipe_cleanup``
================================================

This script was design to erase all the outputs including the `log` files of a module that was already processed.
You could use this function in case you have had an error in the processing of a module, or you simply want to rerun a module that you have already started or that is completed,
it is necessary to delete the already generated files.

.. tabs::

    .. tab:: Usage

        The basic usage requires the same inputs as in mica-pipe:

        **Terminal:**

        .. code-block:: bash
           :linenos:
           :emphasize-lines: 2

           micapipe_cleanup -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory> -ses <session-name>
           -<module-to-erase>

        The desired module(s) to erase should be specify with a flag(s).

        **Optional arguments:**

        -ses 	             Optional flag that indicates the session name (if omitted will manage as SINGLE session)
      	-proc_structural   Deletes volumetric processing derivatives
      	-proc_freesurfer   Deletes Freesurfer recon-all processing outputs
      	-proc_dwi          Deletes DWI-Diffusion weighted images processing derivatives
      	-proc_rsfmri       Deletes Resting state Funtional MRI processing derivatives
      	-post_structural   Deletes Post structural volumetric processing derivatives
      	-SC                Deletes Post tractography files and connectomes
      	-MPC               Deletes Microstructural profiles and covariance analysis derivatives

    .. tab:: Warning

        If you want to erase and rerun a core processing module such as ``post_structural`` or ``proc_structural`` you must rerun all its dependent modules, to ensure outputs consistency.
