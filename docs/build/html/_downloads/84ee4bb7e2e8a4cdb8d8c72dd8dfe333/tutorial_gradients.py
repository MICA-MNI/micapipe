#!/usr/bin/env python
# coding: utf-8

# # Building gradients
# ## Set the enviroment

# In[5]:


# python notebook
#
# Tutorial 3 - Building gradients
# micapipe v0.1.1
#
# Created by RRC on September 2021 (the second year of the pademic)

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


# ## Structural gradient

# In[ ]:


# Set the path to the the structural cortical connectome
cnt_sc_cor = subjectDir + '/dwi/connectomes/' + subjectID + '_space-dwi_atlas-' + atlas + '_desc-iFOD2-40M-SIFT2_full-connectome.txt'

# Load the cortical connectome
mtx_sc = np.loadtxt(cnt_sc_cor, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_scSym = np.triu(mtx_sc,1)+mtx_sc.T


# ## Functional gradient

# In[ ]:


# Set the path to the the functional cortical connectome
cnt_fs = subjectDir + '/func/surfaces/' + subjectID + '_rsfmri_space-fsnative_atlas-' + atlas + '_desc-FC.txt'

# Load the cortical connectome
mtx_fs = np.loadtxt(cnt_fs, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_fcSym = np.triu(mtx_fs,1)+mtx_fs.T


# ## MPC gradient

# In[ ]:


# Set the path to the the MPC cortical connectome
cnt_mpc = subjectDir + '/anat/surfaces/micro_profiles/' + subjectID + '_space-fsnative_atlas-' + atlas + '_desc-MPC.txt'

# Load the cortical connectome
mtx_mpc = np.loadtxt(cnt_mpc, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_mpcSym = np.triu(mtx_mpc,1)+mtx_mpc.T


# ## Geodesic distance

# In[ ]:


# Set the path to the the geodesic distance connectome
cnt_gd = subjectDir + '/anat/surfaces/geo_dist/' + subjectID + '_space-fsnative_atlas-' + atlas + '_GD.txt'

# Load the cortical connectome
mtx_gd = np.loadtxt(cnt_gd, dtype=np.float, delimiter=' ')

