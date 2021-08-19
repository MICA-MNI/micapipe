.. _databases:

.. title:: Processing databases

Processing databases
============================================================

In this section you will find examples of how ``mica-pipe`` optional arguments where used to process a variety of databases with different acquisitions parameters.

Microstructure-Informed Connectomis (MICs)
--------------------------------------------------------

`DOI MICs <https://doi.org/10.1101/2021.08.04.454795>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            tree

    .. tab:: mica-pipe arguments

        .. code-block:: bash
           :linenos:

            mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory> -all


Epilepsy and Cognition (EpiC-UNAM)
--------------------------------------------------------

`EpiC-UNAM <https://github.com/rcruces/2020_cognition_connectomics_TLE>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            tree

    .. tab:: mica-pipe arguments

        .. code-block:: bash
           :linenos:

            mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>


Cambridge Centre for Ageing and Neuroscience (Cam-CAN)
--------------------------------------------------------

`Cam-CAN <https://www.cam-can.org/index.php?content=dataset>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            tree

    .. tab:: mica-pipe arguments

        .. code-block:: bash
           :linenos:

            mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>


Midnight Scan Club MSC
--------------------------------------------------------

`MSC  <https://openneuro.org/datasets/ds000224/versions/1.0.3>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            tree

    .. tab:: mica-pipe arguments

        .. code-block:: bash
           :linenos:

            mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>


SUDMEX_CONN
--------------------------------------------------------

`SUDMEX_CONN  <https://openneuro.org/datasets/ds003346/versions/1.1.1>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            tree

    .. tab:: mica-pipe arguments

        .. code-block:: bash
           :linenos:

            mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>


Auditory localization with 7T fMRI (Audiopath)
--------------------------------------------------------

`Audiopath <https://openneuro.org/datasets/ds001942/versions/1.2.0>`_

.. tabs::

    .. tab:: Database structure

        .. parsed-literal::

            tree

    .. tab:: mica-pipe arguments

        .. code-block:: bash
           :linenos:

            mica-pipe -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>
