.. _qc:

.. title:: Quality control report

Quality control
============================================================
Quality of processed outputs can be run at any moment of the processing.
Group and individual integrated quality control reports will allow you to quickly identify missing files, registration performance, as well as outputs requiring further inspection (workflow).


Individual QC
--------------------------------------------------------

This module generates a html report with detailed information for rapid visualization of completed processed steps, core registrations, and data matrices by parcellation scheme and module.

.. code-block:: bash
   :linenos:
   :emphasize-lines: 2

    mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory> \
    -QC_subj

Group level QC
--------------------------------------------------------

The group level quality check, generates a report with all the completed and processed modules by subjects

.. code-block:: bash
   :linenos:
   :emphasize-lines: 2

    mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory> \
    -QC
