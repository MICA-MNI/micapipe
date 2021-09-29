#!/usr/bin/env python
# coding: utf-8

# # Organization of the outputs

# In[1]:


# python notebook
#
# Tutorial 1 - Main output matrices 
# micapipe v0.1.1
#
# Created by RRC on September 2021 (the second year of the pademic)

# Load required packages
import os
import numpy as np
from nilearn import plotting
import matplotlib as plt

# Set the working directory to the 'out' directory
os.chdir("~/out") # <<<<<<<<<<<< CHANGE THIS PATH

# This variable will be different for each subject
subjectID='sub-HC001_ses-01'           # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
subjectDir='micapipe/sub-HC001/ses-01' # <<<<<<<<<<<< CHANGE THIS SUBJECT's DIRECTORY

# Here we define the atlas 
atlas='schaefer-400' # <<<<<<<<<<<< CHANGE THIS ATLAS


# ## Structural connectomes
# 
# ### Full structural connectome

# In[4]:


# Set the path to the the structural cortical connectome
cnt_sc_cor = subjectDir + '/dwi/connectomes/' + subjectID + '_space-dwi_atlas-' + atlas + '_desc-iFOD2-40M-SIFT2_full-connectome.txt'

# Load the cortical connectome
mtx_sc = np.loadtxt(cnt_sc_cor, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_scSym = np.triu(mtx_sc,1)+mtx_sc.T

# Plot the log matrix
corr_plot = plotting.plot_matrix(np.log(mtx_scSym), figure=(10, 10), labels=None, cmap='Purples', vmin=0, vmax=10)


# ### Full structural connectome edge lengths

# In[5]:


# Set the path to the the structural cortical connectome
cnt_sc_EL = cnt_sc_cor= subjectDir + '/dwi/connectomes/' + subjectID + '_space-dwi_atlas-' + atlas + '_desc-iFOD2-40M-SIFT2_full-edgeLengths.txt'

# Load the cortical connectome
mtx_scEL = np.loadtxt(cnt_sc_EL, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_scELSym = np.triu(mtx_scEL,1)+mtx_scEL.T

# Plot the log matrix
corr_plot = plotting.plot_matrix(mtx_scELSym, figure=(10, 10), labels=None, cmap='Purples', vmin=0, vmax=200)


# ## Resting state functional connectome

# In[6]:


# Set the path to the the functional cortical connectome
cnt_fs = subjectDir + '/func/surfaces/' + subjectID + '_rsfmri_space-fsnative_atlas-' + atlas + '_desc-FC.txt'

# Load the cortical connectome
mtx_fs = np.loadtxt(cnt_fs, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_fcSym = np.triu(mtx_fs,1)+mtx_fs.T

# Plot the matrix
corr_plot = plotting.plot_matrix(mtx_fcSym, figure=(10, 10), labels=None, cmap='Reds', vmin=0, vmax=1)


# ### Time series (ROI x time)

# In[7]:


# Set the path to the the time series file
cnt_time = subjectDir + '/func/surfaces/' + subjectID + '_rsfmri_space-fsnative_atlas-' + atlas + '_desc-timeseries.txt'

# Load the time series
mtx_time = np.loadtxt(cnt_time, dtype=np.float, delimiter=' ')

# Plot as a matrix
corr_plot = plotting.plot_matrix(mtx_time.T, figure=(12, 5), labels=None, cmap='plasma', vmin=-100, vmax=100)


# ## MPC connectomes

# In[8]:


# Set the path to the the MPC cortical connectome
cnt_mpc = subjectDir + '/anat/surfaces/micro_profiles/' + subjectID + '_space-fsnative_atlas-' + atlas + '_desc-MPC.txt'

# Load the cortical connectome
mtx_mpc = np.loadtxt(cnt_mpc, dtype=np.float, delimiter=' ')

# Fill the lower triangle of the matrix
mtx_mpcSym = np.triu(mtx_mpc,1)+mtx_mpc.T

# Plot the matrix
corr_plot = plotting.plot_matrix(mtx_mpcSym, figure=(10, 10), labels=None, cmap='Greens')


# ### Intensity profiles (Profile x ROI)

# In[9]:


# Set the path to the the time series file
cnt_int = subjectDir + '/anat/surfaces/micro_profiles/' + subjectID + '_space-fsnative_atlas-' + atlas + '_desc-intensity_profiles.txt'

# Load the time series
mtx_int = np.loadtxt(cnt_int, dtype=np.float, delimiter=' ')

# Plot as a matrix
corr_plot = plotting.plot_matrix(mtx_int, figure=(20,10), labels=None, cmap='Greens', colorbar=False)


# ## Geodesic distance connectomes

# In[10]:


# Set the path to the the geodesic distance connectome
cnt_gd = subjectDir + '/anat/surfaces/geo_dist/' + subjectID + '_space-fsnative_atlas-' + atlas + '_GD.txt'

# Load the cortical connectome
mtx_gd = np.loadtxt(cnt_gd, dtype=np.float, delimiter=' ')

# Plot the matrix
corr_plot = plotting.plot_matrix(mtx_gd, figure=(10, 10), labels=None, cmap='Blues')


# In[ ]:




