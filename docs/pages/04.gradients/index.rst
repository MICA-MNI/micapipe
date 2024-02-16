.. _gradient:

.. title:: Computing gradient from output matrices

******************
Building gradients
******************

.. contents:: Table of Contents

This section describes how to generate macroscale gradient mapping from the output matrices of micapipe. The matrices are the same as in the `Main output matrices <../04.matrices/index.html>`_ tutorial.

For this tutorial we will map each modality of a single subject using `BrainSpace <https://brainspace.readthedocs.io/en/latest/python_doc/auto_examples/index.html>`_, a ``python`` based library.

For further information about how to use `BrainSpace` and macroscale gradient mapping and analysis of neuroimaging and connectome level data visit their documentation here: `BrainSpace documentation <https://brainspace.readthedocs.io/en/latest/index.html>`_

Additionally the libraries ``os``, ``numpy``, ``nibabel`` and ``nibabel`` will be used.

As in previous examples, we will use the subject ``HC001``, session ``01`` from the MICs dataset, and all paths will be relative to the subject directory or ``out/micapipe_v0.2.0/sub-HC001_ses01/`` and the atlas ``schaefer-400``.

In this tutorial we will only plot the first 3 components of the diffusion map embedding (aka gradients).

.. admonition:: diffusion map embedding
    For further references check:
    - Coifman & Lafon. Diffusion maps. Applied and computational harmonic analysis. 2006 Jul 1;21(1):5-30. https://doi.org/10.1016/j.acha.2006.04.006
    - Vos de Wael, Benkarim, et al. BrainSpace: a toolbox for the analysis of macroscale gradients in neuroimaging and connectomics datasets. Commun Biol 3, 103 (2020). https://doi.org/10.1038/s42003-020-0794-7


python environment
------------------------------------------------------------

The first step is to set the python environment and all the variables relative to the subject you're interested to visualize.

.. tabs::

   .. code-tab:: py
      :linenos:

      # Set the environment
      import os
      import glob
      import numpy as np
      import nibabel as nib
      from brainspace.plotting import plot_hemispheres
      from brainspace.mesh.mesh_io import read_surface
      from brainspace.datasets import load_conte69
      from brainspace.gradient import GradientMaps
      from brainspace.utils.parcellation import map_to_labels
      import matplotlib.pyplot as plt

      # Set the working directory to the 'out' directory
      out='/data_/mica3/BIDS_MICs/derivatives'
      os.chdir(out)     # <<<<<<<<<<<< CHANGE THIS PATH

      # This variable will be different for each subject
      sub='HC001' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
      ses='01'    # <<<<<<<<<<<< CHANGE THIS SUBJECT's SESSION
      subjectID=f'sub-{sub}_ses-{ses}'
      subjectDir=f'micapipe_v0.2.0/sub-{sub}/ses-{ses}'

      # Here we define the atlas
      atlas='schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS

      # Path to MICAPIPE from global enviroment
      micapipe=os.popen("echo $MICAPIPE").read()[:-1] # <<<<<<<<<<<< CHANGE THIS PATH

Loading the surfaces
============================================================

.. tabs::

   .. code-tab:: py
     :linenos:

      # Load fsLR-32k
      f32k_lh = read_surface(f'{micapipe}/surfaces/fsLR-32k.L.surf.gii', itype='gii')
      f32k_rh = read_surface(f'{micapipe}/surfaces/fsLR-32k.R.surf.gii', itype='gii')

      # Load fsaverage5
      fs5_lh = read_surface(f'{micapipe}/surfaces/fsaverage5/surf/lh.pial', itype='fs')
      fs5_rh = read_surface(f'{micapipe}/surfaces/fsaverage5/surf/rh.pial', itype='fs')

      # Load LEFT annotation file in fsaverage5
      annot_lh_fs5= nib.freesurfer.read_annot(f'{micapipe}/parcellations/lh.{atlas}_mics.annot')

      # Unique number of labels of a given atlas
      Ndim = max(np.unique(annot_lh_fs5[0]))

      # Load RIGHT annotation file in fsaverage5
      annot_rh_fs5= nib.freesurfer.read_annot(f'{micapipe}/parcellations/rh.{atlas}_mics.annot')[0]+Ndim

      # replace with 0 the medial wall of the right labels
      annot_rh_fs5 = np.where(annot_rh_fs5==Ndim, 0, annot_rh_fs5)

      # fsaverage5 labels
      labels_fs5 = np.concatenate((annot_lh_fs5[0], annot_rh_fs5), axis=0)

      # Mask of the medial wall on fsaverage 5
      mask_fs5 = labels_fs5 != 0

      # Read label for fsLR-32k
      labels_f32k = np.loadtxt(open(f'{micapipe}/parcellations/{atlas}_conte69.csv'), dtype=int)

      # mask of the medial wall
      mask_f32k = labels_f32k != 0

Global variables
============================================================

.. tabs::

   .. code-tab:: py
      :linenos:

      # Number of gradients to calculate
      Ngrad=10

      # Number of gradients to plot
      Nplot=3

      # Labels for plotting based on Nplot
      labels=['G'+str(x) for x in list(range(1,Nplot+1))]

Gradients from atlas based connectomes: single subject
------------------------------------------------------------

Geodesic distance
============================================================

.. tabs::

   .. tab:: Python

        **Load and slice the GD matrix**

        .. code-block:: python
           :linenos:

            # Set the path to the the geodesic distance connectome
            gd_file = f'{subjectDir}/dist/{subjectID}_atlas-{atlas}_GD.shape.gii'

            # Load the cortical connectome
            mtx_gd = nib.load(gd_file).darrays[0].data

            # Remove the Mediall Wall
            mtx_gd = np.delete(np.delete(mtx_gd, 0, axis=0), 0, axis=1)
            GD = np.delete(np.delete(mtx_gd, Ndim, axis=0), Ndim, axis=1)


        **Calculate the GD gradients**

        .. code-block:: python
           :linenos:

            # GD Left hemi
            gm_GD_L = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
            gm_GD_L.fit(GD[0:Ndim, 0:Ndim], sparsity=0.8)

            # GD Right hemi
            gm_GD_R = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
            gm_GD_R.fit(GD[Ndim:Ndim*2, Ndim:Ndim*2], sparsity=0.85, reference=gm_GD_L.gradients_)

        **Plot the GD gradients**

        .. code-block:: python
           :linenos:

            # plot the gradients
            g1=gm_GD_L.gradients_[:, 0]
            g2=gm_GD_L.gradients_[:, 1]
            g3=gm_GD_L.gradients_[:, 2]

            # plot the gradients
            g1R=gm_GD_R.aligned_[:, 0]
            g2R=gm_GD_R.aligned_[:, 1]
            g3R=gm_GD_R.aligned_[:, 2]

            # Creating figure
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111, projection="3d")

            # Creating plot
            ax.scatter3D(g1, g2, g3, color = 'dodgerblue')
            ax.scatter3D(g1R, g2R, g3R, color = 'teal', marker='v')
            plt.title("Structural gradient")
            ax.legend(['Left GD', 'Right GD'])

            ax.set_xlabel('Grad 1')
            ax.set_ylabel('Grad 2')
            ax.set_zlabel('Grad 3')

            # Remove the outer box lines
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            # Show plot
            plt.show()

        .. figure:: gd_scatter.png

        **GD gradients on** ``fsaverage5`` **surface**

        .. code-block:: python
           :linenos:

            # Left and right gradients concatenated
            GD_gradients = np.concatenate((gm_GD_L.gradients_, gm_GD_R.aligned_), axis=0)

            # Map gradients to original parcels
            grad = [None] * Nplot
            for i, g in enumerate(GD_gradients.T[0:Nplot,:]):
                grad[i] = map_to_labels(g, labels_fs5,  fill=np.nan, mask=mask_fs5)

            # Plot Gradients RdYlBu
            plot_hemispheres(fs5_lh, fs5_rh, array_name=grad, size=(1000, 600), cmap='coolwarm',
                             embed_nb=True,  label_text={'left':labels}, color_bar='left',
                             zoom=1.25, nan_color=(1, 1, 1, 1), color_range = 'sym' )

        .. figure:: gd_fs5.png

        **GD gradients to** ``fsLR-32k`` **surface**

        .. code-block:: python
           :linenos:

            # Map gradients to original parcels
            grad = [None] * Nplot
            for i, g in enumerate(GD_gradients.T[0:Nplot,:]):
                grad[i] = map_to_labels(g, labels_f32k, fill=np.nan, mask=mask_f32k)

            # Plot Gradients
            plot_hemispheres(f32k_lh, f32k_rh, array_name=grad, size=(1000, 600), cmap='coolwarm',
                             embed_nb=True,  label_text={'left':labels}, color_bar='left',
                             zoom=1.25, nan_color=(1, 1, 1, 1))


        .. figure:: gd_f32k.png

Structural gradients
============================================================

.. tabs::

   .. tab:: Python

        **Load and slice the structural matrix**

        .. code-block:: python
           :linenos:

            # Set the path to the the structural cortical connectome
            sc_file = f'{subjectDir}/dwi/connectomes/{subjectID}_space-dwi_atlas-{atlas}_desc-iFOD2-40M-SIFT2_full-connectome.shape.gii'

            # Load the cortical connectome
            mtx_sc = nib.load(sc_file).darrays[0].data

            # Fill the lower triangle of the matrix
            mtx_sc = np.log(np.triu(mtx_sc,1)+mtx_sc.T)
            mtx_sc[np.isneginf(mtx_sc)] = 0

            # Slice the connectome to use only cortical nodes
            SC = mtx_sc[49:, 49:]
            SC = np.delete(np.delete(SC, 200, axis=0), 200, axis=1)


        **Calculate the structural gradients**

        .. code-block:: python
           :linenos:

            # SC Left hemi
            gm_SC_L = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
            gm_SC_L.fit(SC[0:Ndim, 0:Ndim], sparsity=0.9)

            # SC Right hemi
            gm_SC_R = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
            gm_SC_R.fit(SC[Ndim:Ndim*2, Ndim:Ndim*2], sparsity=0.9, reference=gm_SC_L.gradients_)

        **Plot the structural gradients**

        .. code-block:: python
           :linenos:

            # plot the left gradients
            g1=gm_SC_L.gradients_[:, 0]
            g2=gm_SC_L.gradients_[:, 1]
            g3=gm_SC_L.gradients_[:, 2]

            # plot the right gradients
            g1R=gm_SC_R.aligned_[:, 0]
            g2R=gm_SC_R.aligned_[:, 1]
            g3R=gm_SC_R.aligned_[:, 2]

            # Creating figure
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111, projection="3d")

            # Creating plot
            ax.scatter3D(g1, g2, g3, color = 'purple')
            ax.scatter3D(g1R, g2R, g3R, color = 'slateblue', marker='v')
            plt.title("Structural gradient")
            ax.legend(['Left SC', 'Right SC'])

            ax.set_xlabel('Grad 1')
            ax.set_ylabel('Grad 2')
            ax.set_zlabel('Grad 3')

            # Remove the outer box lines
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            # Show plot
            plt.show()

        .. figure:: sc_scatter.png

        **Structural gradients on** ``fsLR-32k`` **surface**

        .. code-block:: python
           :linenos:

            # Left and right gradients concatenated
            SC_gradients = np.concatenate((gm_SC_L.gradients_, gm_SC_R.aligned_), axis=0)

            # Map gradients to original parcels
            grad = [None] * Nplot
            for i, g in enumerate(SC_gradients.T[0:Nplot,:]):
            grad[i] = map_to_labels(g, labels_f32k, fill=np.nan, mask=mask_f32k)

            # Plot Gradients
            plot_hemispheres(f32k_lh, f32k_rh, array_name=grad, size=(1000, 600), cmap='coolwarm',
                 embed_nb=True,  label_text={'left':labels}, color_bar='left',
                 zoom=1.25, nan_color=(1, 1, 1, 1), color_range = 'sym' )


        .. figure:: sc_f32k.png

Functional gradients
============================================================

.. tabs::

   .. tab:: Python

        **Load and slice the functional matrix**

        .. code-block:: python
           :linenos:

            # acquisitions
            func_acq='desc-se_task-rest_acq-AP_bold'
            fc_file = f'{subjectDir}/func/{func_acq}/surf/{subjectID}_surf-fsLR-32k_atlas-{atlas}_desc-FC.shape.gii'

            # Load the cortical connectome
            mtx_fs = nib.load(fc_file).darrays[0].data

            # slice the matrix to keep only the cortical ROIs
            FC = mtx_fs[49:, 49:]
            #FC = np.delete(np.delete(FC, Ndim, axis=0), Ndim, axis=1)

            # Fischer transformation
            FCz = np.arctanh(FC)

            # replace inf with 0
            FCz[~np.isfinite(FCz)] = 0

            # Mirror the matrix
            FCz = np.triu(FCz,1)+FCz.T

        **Calculate the functional gradients**

        .. code-block:: python
           :linenos:

            # Calculate the gradients
            gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
            gm.fit(FCz, sparsity=0.85)

        **Plot the functional gradients**

        .. code-block:: python
           :linenos:

            # Plot the gradients
            g1 = gm.gradients_[:, 0]
            g2 = gm.gradients_[:, 1]
            g3 = gm.gradients_[:, 2]

            # Creating figure
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111, projection="3d")

            # Creating plot
            ax.scatter3D(g1, g2, g3, color='red')
            plt.title("Functional gradient")

            ax.set_xlabel('Grad 1')
            ax.set_ylabel('Grad 2')
            ax.set_zlabel('Grad 3')

            # Remove the outer box lines
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            # Show plot
            plt.show()

        .. figure:: fc_scatter.png

        **Functional gradients on** ``fsLR-32k`` **surface**

        .. code-block:: python
           :linenos:

            # Map gradients to original parcels
            grad = [None] * Nplot
            for i, g in enumerate(gm.gradients_.T[0:Nplot,:]):
                grad[i] = map_to_labels(g, labels_f32k, fill=np.nan, mask=mask_f32k)

            # Plot Gradients coolwarm
            plot_hemispheres(f32k_lh, f32k_rh, array_name=grad, size=(1000, 600), cmap='coolwarm',
                             embed_nb=True,  label_text={'left':labels}, color_bar='left',
                             zoom=1.25, nan_color=(1, 1, 1, 1), color_range = 'sym')

        .. figure:: fc_f32k.png

MPC gradients
============================================================

.. tabs::

   .. tab:: Python

        **Function to load MPC**

        .. code-block:: python
           :linenos:

            # Define a function to load and process the MPC matrices
            def load_mpc(File, Ndim):
                """Loads and process a MPC"""

                # Load file
                mpc = nib.load(File).darrays[0].data

                # Mirror the lower triangle
                mpc = np.triu(mpc,1)+mpc.T

                # Replace infinite values with epsilon
                mpc[~np.isfinite(mpc)] = np.finfo(float).eps

                # Replace 0 with epsilon
                mpc[mpc==0] = np.finfo(float).eps

                # Remove the medial wall
                mpc = np.delete(np.delete(mpc, 0, axis=0), 0, axis=1)
                mpc = np.delete(np.delete(mpc, Ndim, axis=0), Ndim, axis=1)

                # retun the MPC
                return(mpc)

        **List and load the MPC matrix**

        .. code-block:: python
           :linenos:

            # Set the path to the the MPC cortical connectome
            mpc_acq='acq-T1map'
            mpc_file = f'{subjectDir}/mpc/{mpc_acq}/{subjectID}_atlas-{atlas}_desc-MPC.shape.gii'

            # Load the cortical connectome
            mpc = load_mpc(mpc_file, Ndim)


        **Calculate the MPC gradients**

        .. code-block:: python
           :linenos:

            # Calculate the gradients
            gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
            gm.fit(mpc, sparsity=0)


        **Plot the MPC gradients**

        .. code-block:: python
           :linenos:

            # Plot the gradients
            g1 = gm.gradients_[:, 0]
            g2 = gm.gradients_[:, 1]
            g3 = gm.gradients_[:, 2]

            # Creating figure
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111, projection="3d")

            # Creating plot
            ax.scatter3D(g1, g2, g3, color = 'green')
            plt.title("MPC gradient")

            ax.set_xlabel('Grad 1')
            ax.set_ylabel('Grad 2')
            ax.set_zlabel('Grad 3')

            # Remove the outer box lines
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            # Show plot
            plt.show()

        .. figure:: mpc_scatter.png

        **MPC gradients on** ``fsLR-32k`` **surface**

        .. code-block:: python
           :linenos:

            # Map gradients to original parcels
            grad = [None] * Nplot
            for i, g in enumerate(gm.gradients_.T[0:Nplot,:]):
                grad[i] = map_to_labels(g, labels_f32k, fill=np.nan, mask=mask_f32k)

            # Plot Gradients
            plot_hemispheres(f32k_lh, f32k_rh, array_name=grad, size=(1000, 600), cmap='coolwarm',
                             embed_nb=True,  label_text={'left':labels}, color_bar='left',
                             zoom=1.25, nan_color=(1, 1, 1, 1), color_range = 'sym' )

        .. figure:: mpc_f32k.png

Load all matrices from a dataset processed
------------------------------------------------------------

1. Start by generating a list of files using regular expressions for matrices with a consistent structure. Specifically, we'll focus on loading the ``T1map MPC`` connectome data for ``schaefer-400`` from the MPC directory.

2. Create an empty three-dimensional array with dimensions ``{ROI * ROI * subjects}``.

3. Load each matrix iteratively and populate the array with the data.

4. Once the array is populated, perform computations on it. In this case, we'll calculate the group mean connectome.

5. Use the group mean connectome to compute the group mean diffusion map for the ``T1map MPC``.

6. Finally, visualize the results by plotting the first three gradients (eigen vectors) of the group mean diffusion map on a surface ``fsLR-32k``.

MPC gradients: ALL subjects mean
============================================================

.. tabs::

   .. tab:: Python

        **Load all the MPC matrices**

        .. code-block:: python
           :linenos:

            # MPC T1map acquisition
            mpc_acq='T1map'

            # 1. List all the matrices from all subjects
            mpc_files = sorted(glob.glob(f'micapipe_v0.2.0/sub-PX*/ses-*/mpc/acq-{mpc_acq}/*_atlas-{atlas}_desc-MPC.shape.gii'))
            N = len(mpc_files)
            print(f"Number of subjects's MPC: {N}")

            # 2. Empty 3D array to load the data
            mpc_all=np.empty([Ndim*2, Ndim*2, len(mpc_files)], dtype=float)

            # 3. Load all the  MPC matrices
            for i, f in enumerate(mpc_files):
                mpc_all[:,:,i] = load_mpc(f, Ndim)

            # Print the shape of the 3D-array: {roi * roi * subjects}
            mpc_all.shape

        **Calculate the mean group MPC gradients**

        .. code-block:: python
           :linenos:

            # 4. Mean group MPC across all subjects (z-axis)
            mpc_all_mean = np.mean(mpc_all, axis=2)

            # Calculate the gradients
            gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
            gm.fit(mpc_all_mean, sparsity=0)

        **Plot the mean group MPC gradients**

        .. code-block:: python
           :linenos:

            # Plot the gradients
            g1 = gm.gradients_[:, 0]
            g2 = gm.gradients_[:, 1]
            g3 = gm.gradients_[:, 2]

            # Creating figure
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111, projection="3d")

            # Creating plot
            ax.scatter3D(g1, g2, g3, color = 'green')
            plt.title("MPC gradient")

            ax.set_xlabel('Grad 1')
            ax.set_ylabel('Grad 2')
            ax.set_zlabel('Grad 3')

            # Remove the outer box lines
            ax.xaxis.pane.fill = False
            ax.yaxis.pane.fill = False
            ax.zaxis.pane.fill = False

            # Show plot
            plt.show()

        .. figure:: mpc-all_scatter.png

        **Mean group MPC gradients on** ``fsLR-32k`` **surface**

        .. code-block:: python
           :linenos:

            # Map gradients to original parcels
            grad = [None] * Nplot
            for i, g in enumerate(gm.gradients_.T[0:Nplot,:]):
                grad[i] = map_to_labels(g, labels_f32k, fill=np.nan, mask=mask_f32k)

            # Plot Gradients
            plot_hemispheres(f32k_lh, f32k_rh, array_name=grad, size=(1000, 600), cmap='coolwarm',
                             embed_nb=True,  label_text={'left':labels}, color_bar='left',
                             zoom=1.25, nan_color=(1, 1, 1, 1), color_range = 'sym' )

        .. figure:: mpc-all_f32k.png

Download the code!: atlas based gradients
------------------------------------------------------------

:download:`Python Jupyter notebook: 'tutorial_gradients.ipynb' <tutorial_gradients.ipynb>`

:download:`Python source code: 'tutorial_gradients.py' <tutorial_gradients.py>`

*********************
``fsLR-5k`` gradients
*********************

Set the environment
------------------------------------------------------------

.. tabs::

   .. tab:: Python

        .. code-block:: python
           :linenos:

            # Set the environment
            import os
            import glob
            import numpy as np
            import nibabel as nib
            from brainspace.plotting import plot_hemispheres
            from brainspace.mesh.mesh_io import read_surface
            from brainspace.datasets import load_conte69
            from brainspace.gradient import GradientMaps
            from brainspace.utils.parcellation import map_to_labels
            import matplotlib.pyplot as plt

            # Set the working directory to the 'out' directory
            out='/data_/mica3/BIDS_MICs/derivatives'  # <<<<<<<<<<<< CHANGE THIS PATH
            os.chdir(f'{out}/micapipe_v0.2.0')

            # This variable will be different for each subject
            sub='HC001' # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
            ses='01'    # <<<<<<<<<<<< CHANGE THIS SUBJECT's SESSION
            subjectID=f'sub-{sub}_ses-{ses}'
            subjectDir=f'micapipe_v0.2.0/sub-{sub}/ses-{ses}'

            # Path to MICAPIPE from global enviroment
            micapipe=os.popen("echo $MICAPIPE").read()[:-1] # <<<<<<<<<<<< CHANGE THIS PATH

        **Load the surfaces**

        .. code-block:: python
           :linenos:

            # Load fsLR-5k inflated surface
            micapipe='/data_/mica1/01_programs/micapipe-v0.2.0'
            f5k_lhi = read_surface(micapipe + '/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
            f5k_rhi = read_surface(micapipe + '/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')

            # fsLR-5k mask
            mask_lh = nib.load(micapipe + '/surfaces/fsLR-5k.L.mask.shape.gii').darrays[0].data
            mask_rh = nib.load(micapipe + '/surfaces/fsLR-5k.R.mask.shape.gii').darrays[0].data
            mask_5k = np.concatenate((mask_lh, mask_rh), axis=0)

        **Functions to load** ``fsLR-5k`` **connectomes**

        .. code-block:: python
           :linenos:

            # Define functions to load GD, SC, FC and MPC fsLR-32k
            def load_mpc(File):
                """Loads and process a MPC"""

                # Load file
                mpc = nib.load(File).darrays[0].data

                # Mirror the lower triangle
                mpc = np.triu(mpc,1)+mpc.T

                # Replace infinite values with epsilon
                mpc[~np.isfinite(mpc)] = np.finfo(float).eps

                # Replace 0 with epsilon
                mpc[mpc==0] = np.finfo(float).eps

                # retun the MPC
                return(mpc)

            def load_gd(File):
                """Loads and process a GD"""

                # load the matrix
                mtx_gd = nib.load(File).darrays[0].data

                return mtx_gd

            def load_fc(File):
                """Loads and process a functional connectome"""

                # load the matrix
                FC = nib.load(File).darrays[0].data

                # Fisher transform
                FCz = np.arctanh(FC)

                # replace inf with 0
                FCz[~np.isfinite(FCz)] = 0

                # Mirror the matrix
                FCz = np.triu(FCz,1)+FCz.T
                return FCz

            def load_sc(File):
                """Loads and process a structura connectome"""

                # load the matrix
                mtx_sc = nib.load(File).darrays[0].data

                # Mirror the matrix
                mtx_sc = np.triu(mtx_sc,1)+mtx_sc.T

                return mtx_sc

        **Functions to calculate** ``fsLR-5k`` **diffusion maps**

        .. code-block:: python
           :linenos:

            # Gradients aka eigen vector of the diffusion map embedding
            def fslr5k_dm_lr(mtx, mask_5k, Ngrad=3, log=True, S=0):
                """
                Create the gradients from the SC or GD matrices.
                Use log=False for GD gradients
                """
                if log != True:
                    mtx_log = mtx
                else:
                    # log transform the connectome
                    mtx_log = np.log(mtx)

                # Replace NaN with 0
                mtx_log[np.isnan(mtx_log)] = 0

                # Replace negative infinite with 0
                mtx_log[np.isneginf(mtx_log)] = 0

                # Replace infinite with 0
                mtx_log[~np.isfinite(mtx_log)] = 0

                # replace 0 values with almost 0
                mtx_log[mtx_log==0] = np.finfo(float).eps

                # Left and right mask
                indx_L = np.where(mask_5k[0:4842]==1)[0]
                indx_R = np.where(mask_5k[4842:9684]==1)[0]

                # Left and right SC
                mtx_L = mtx_log[0:4842, 0:4842]
                mtx_R = mtx_log[4842:9684, 4842:9684]

                # Slice the matrix
                mtx_L_masked = mtx_L[indx_L, :]
                mtx_L_masked = mtx_L_masked[:, indx_L]
                mtx_R_masked = mtx_R[indx_R, :]
                mtx_R_masked = mtx_R_masked[:, indx_R]

                # mtx Left hemi
                mtx_L = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
                mtx_L.fit(mtx_L_masked, sparsity=S)

                # mtx Right hemi
                mtx_R = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
                mtx_R.fit(mtx_R_masked, sparsity=S, reference=mtx_L.gradients_)

                # Left and right gradients concatenated
                mtx_gradients = np.concatenate((mtx_L.gradients_, mtx_R.aligned_), axis=0)

                # Boolean mask
                mask_surf = mask_5k != 0

                # Get the index of the non medial wall regions
                indx = np.where(mask_5k==1)[0]

                # Map gradients to surface
                grad = [None] * Ngrad
                for i, g in enumerate(mtx_gradients.T[0:Ngrad,:]):
                    # create a new array filled with NaN values
                    g_nan = np.full(mask_surf.shape, np.nan)
                    g_nan[indx] = g
                    grad[i] = g_nan

                return(mtx_gradients, grad)

            def fslr5k_dm(mtx, mask, Ngrad=3, S=0.9):
                """Create the gradients from the MPC matrix
                    S=sparcity, by default is 0.9
                """
                # Cleanup before diffusion embeding
                mtx[~np.isfinite(mtx)] = 0
                mtx[np.isnan(mtx)] = 0
                mtx[mtx==0] = np.finfo(float).eps

                # Get the index of the non medial wall regions
                indx = np.where(mask==1)[0]

                # Slice the matrix
                mtx_masked = mtx[indx, :]
                mtx_masked = mtx_masked[:, indx]

                # Calculate the gradients
                gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
                gm.fit(mtx_masked, sparsity=S)

                # Map gradients to surface
                grad = [None] * Ngrad

                # Boolean mask
                mask_surf = mask != 0

                for i, g in enumerate(gm.gradients_.T[0:Ngrad,:]):

                    # create a new array filled with NaN values
                    g_nan = np.full(mask_surf.shape, np.nan)
                    g_nan[indx] = g
                    grad[i] = g_nan

                return(gm, grad)

        **Global variables**

        .. code-block:: python
           :linenos:

            # Number of vertices of the fsLR-5k matrices (per hemisphere)
            N5k = 9684

            # Number of gradients to calculate
            Ngrad=10

            # Number of gradients to plot
            Nplot=3

            # Labels for plotting based on Nplot
            labels=['G'+str(x) for x in list(range(1,Nplot+1))]

Geodesic distance: single subject ``fsLR-5k``
------------------------------------------------------------

.. tabs::

   .. tab:: Python

        .. code-block:: python
           :linenos:

            # List the file
            gd_file = glob.glob(f"sub-{sub}/ses-{ses}/dist/*_surf-fsLR-5k_GD.shape.gii")

            # Loads the GD fsLR-5k matrix
            gd_5k = load_gd(gd_file[0])

            # Calculate the gradients
            gd_dm, grad = fslr5k_dm_lr(gd_5k, mask_5k, Ngrad=Ngrad, log=False, S=0.85)

            # plot the gradients
            plot_hemispheres(f5k_lhi, f5k_rhi, array_name=grad[0:Nplot], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
              zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym',
              color_bar='right', label_text={'left': labels})

        .. figure:: gd_f5k.png

Structual connectome: single subject ``fsLR-5k``
------------------------------------------------------------

.. tabs::

   .. tab:: Python

        .. code-block:: python
           :linenos:

            # List the file
            sc_file = sorted(glob.glob(f"sub-{sub}/ses-{ses}/dwi/connectomes/*_surf-fsLR-5k_desc-iFOD2-40M-SIFT2_full-connectome.shape.gii"))

            # Loads the SC fsLR-5k matrix
            sc_5k = load_sc(sc_file[0])

            # Calculate the gradients
            sc_dm, grad = fslr5k_dm_lr(sc_5k, mask_5k, Ngrad=Ngrad, log=False, S=0.9)

            # PLot the gradients (G2-G4)
            plot_hemispheres(f5k_lhi, f5k_rhi, array_name=grad[1:Nplot+1], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
              zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym',
              color_bar='right', label_text={'left': labels})

        .. figure:: sc_f5k.png

Functional connectome: single subject ``fsLR-5k``
------------------------------------------------------------

.. tabs::

   .. tab:: Python

        .. code-block:: python
           :linenos:

            # List the file
            func_acq='desc-se_task-rest_acq-AP_bold'
            fc_file = sorted(glob.glob(f"sub-{sub}/ses-{ses}/func/{func_acq}/surf/*_surf-fsLR-5k_desc-FC.shape.gii"))

            # Loads the FC fsLR-5k matrix
            fc_5k = load_fc(fc_file[0])

            # Calculate the gradients
            fc_dm, grad = fslr5k_dm(fc_5k, mask_5k, Ngrad=Ngrad, S=0.9)

            # plot the gradients
            plot_hemispheres(f5k_lhi, f5k_rhi, array_name=grad[0:Nplot], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
              zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym',
              color_bar='right', label_text={'left': labels})

        .. figure:: fc_f5k.png

MPC T1map: single subject ``fsLR-5k``
------------------------------------------------------------

.. tabs::

   .. tab:: Python

        .. code-block:: python
           :linenos:

            # MPC T1map acquisition and file
            mpc_acq='T1map'
            mpc_file = sorted(glob.glob(f"sub-{sub}/ses-{ses}/mpc/acq-{mpc_acq}/*surf-fsLR-5k_desc-MPC.shape.gii"))

            # Loads the MPC fsLR-5k matrix
            mpc_5k = load_mpc(mpc_file[0])

            # Calculate the gradients (diffusion map)
            mpc_dm, grad = fslr5k_dm(mpc_5k, mask_5k, Ngrad=Ngrad, Smooth=True, S=0.9)

            # Plot the gradients
            plot_hemispheres(f5k_lhi, f5k_rhi, array_name=grad[0:Nplot], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
              zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym',
              color_bar='right', label_text={'left': labels})

        .. figure:: mpc_f5k.png

MPC T1map: ALL subjects ``fsLR-5k``
------------------------------------------------------------

**Load all matrices from a dataset processed**

1. Start by generating a list of files using regular expressions for matrices with a consistent structure. Specifically, we'll focus on loading the ``T1map MPC`` connectome data for ``fsLR-5k`` from the MPC directory.

2. Create an empty three-dimensional array with dimensions ``{ROI * ROI * vertices}``.

3. Load each matrix iteratively and populate the array with the data.

4. Once the array is populated, perform computations on it. In this case, we'll calculate the group mean connectome.

5. Use the group mean connectome to compute the group mean diffusion map for the ``T1map MPC``.

6. Finally, visualize the results by plotting the first three gradients (eigen vectors) of the group mean diffusion map on a surface ``fsLR-5k``.

.. tabs::

   .. tab:: Python

        .. code-block:: python
           :linenos:

            # MPC T1map acquisition
            mpc_acq='T1map'

            # List all the matrices from all subjects
            mpc_file = sorted(glob.glob(f"sub-PX*/ses-01/mpc/acq-{mpc_acq}/*surf-fsLR-5k_desc-MPC.shape.gii"))
            N = len(mpc_file)
            print(f"Number of subjects's MPC: {N}")

            # Loads all the MPC fsLR-5k matrices
            mpc_5k_all=np.empty([N5k, N5k, len(mpc_file)], dtype=float)
            for i, f in enumerate(mpc_file):
                mpc_5k_all[:,:,i] = load_mpc(f)

            # Print the shape of the array: {vertices * vertices * subjects}
            mpc_5k_all.shape

            # Mean group MPC across all subjects (z-axis)
            mpc_5k_mean = np.mean(mpc_5k_all, axis=2)

            # Calculate the gradients (diffusion map)
            mpc_dm, grad = fslr5k_dm(mpc_5k_mean, mask_5k, Ngrad=Ngrad, S=0)

            # Plot the gradients
            plot_hemispheres(f5k_lhi, f5k_rhi, array_name=grad[0:Nplot], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
              zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym',
              color_bar='right', label_text={'left': labels})


        .. figure:: mpc-all_f5k.png

Download the code!: ``fsLR-5k`` gradients
------------------------------------------------------------

:download:`Python Jupyter notebook: 'tutorial_fsLR-5k.ipynb' <tutorial_fsLR-5k.ipynb>`

:download:`Python source code: 'tutorial_fsLR-5k.py' <tutorial_fsLR-5k.py>`
