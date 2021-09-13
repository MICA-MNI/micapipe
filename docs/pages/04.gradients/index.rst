.. _gradient:

.. title:: Computing gradient from output matrices

******************
Building gradients
******************

`BrainSpace tutorials <https://brainspace.readthedocs.io/en/latest/python_doc/auto_examples/index.html>`_

Set the enviroment
============================================================

.. tabs::

   .. code-tab:: py

    # Set the environment
    import os
    import numpy as np
    import matplotlib as plt
    import nibabel as nb
    from nibabel.freesurfer.mghformat import load
    from brainspace.plotting import plot_hemispheres
    from brainspace.mesh.mesh_io import read_surface
    from brainspace.datasets import load_conte69

    # Set the working directory to the 'out' directory
    os.chdir("~/out") # <<<<<<<<<<<< CHANGE THIS PATH

    # This variable will be different for each subject
    subjectID='sub-HC001_ses-01'           # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
    subjectDir='micapipe/sub-HC001/ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's DIRECTORY

    # Here we define the atlas
    atlas='schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS


Structural gradient
============================================================

.. tabs::

   .. code-tab:: py

    # Load the connectome

    # Calculate the gradients

    # Visualize the G1, G2 y G3

    # Visualize the values on surface


Functional gradient
============================================================

.. tabs::

   .. code-tab:: py

    # Load the connectome

    # Calculate the gradients

    # Visualize the G1, G2 y G3

    # Visualize the values on surface


MPC gradient
============================================================

.. tabs::

   .. code-tab:: py

    # Load the connectome

    # Calculate the gradients

    # Visualize the G1, G2 y G3

    # Visualize the values on surface


Geodesic distance gradient
============================================================

.. tabs::

   .. code-tab:: py

    # Load the connectome

    # Calculate the gradients

    # Visualize the G1, G2 y G3

    # Visualize the values on surface



Download code examples: Gradients
--------------------------------------------------------

:download:`Python Jupyter notebook: 'tutorial_gradients.ipynb' <tutorial_gradients.ipynb>`

:download:`Python source code: 'tutorial_gradients.py' <tutorial_gradients.py>`
