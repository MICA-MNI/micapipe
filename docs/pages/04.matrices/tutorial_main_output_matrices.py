#!/usr/bin/env python
# coding: utf-8

# # Organization of the outputs

# In[1]:


# python notebook
#
# Tutorial 1 - Main output matrices 
# micapipe v0.2.3

# Load required packages
import os
import numpy as np
import nibabel as nib
from nilearn import plotting
import matplotlib as plt

# Set the working directory to the 'out' directory
os.chdir("/data_/mica3/BIDS_MICs/derivatives") # <<<<<<<<<<<< CHANGE THIS PATH TO YOUR OUT DIRECTORY

# This variable will be different for each subject
sub='HC001'           # <<<<<<<<<<<< CHANGE THIS SUBJECT's ID
ses='01'              # <<<<<<<<<<<< CHANGE THIS SESSION
subjectID=f'sub-{sub}_ses-{ses}'           
subjectDir=f'micapipe_v0.2.0/sub-{sub}/ses-{ses}' 

# Here we define the atlas 
atlas='schaefer-400'


# ## Structural connectomes
# 
# ### Full structural connectome

# In[2]:


# Set the path to the the structural cortical connectome
cnt_sc_cor = f'{subjectDir}/dwi/connectomes/{subjectID}_space-dwi_atlas-{atlas}_desc-iFOD2-40M-SIFT2_full-connectome.shape.gii'

# Load the cortical connectome
mtx_sc = nib.load(cnt_sc_cor).darrays[0].data

# Fill the lower triangle of the matrix
mtx_scSym = np.triu(mtx_sc,1)+mtx_sc.T

# Plot the log matrix
corr_plot = plotting.plot_matrix(np.log(mtx_scSym), figure=(10, 10), labels=None, cmap='Purples', vmin=0, vmax=10)


# ### Full structural connectome edge lengths

# In[3]:


# Set the path to the the structural cortical connectome
cnt_sc_EL = cnt_sc_cor= f'{subjectDir}/dwi/connectomes/{subjectID}_space-dwi_atlas-{atlas}_desc-iFOD2-40M-SIFT2_full-edgeLengths.shape.gii'

# Load the cortical connectome
mtx_scEL = nib.load(cnt_sc_EL).darrays[0].data

# Fill the lower triangle of the matrix
mtx_scELSym = np.triu(mtx_scEL,1)+mtx_scEL.T

# Plot the log matrix
corr_plot = plotting.plot_matrix(mtx_scELSym, figure=(10, 10), labels=None, cmap='Purples', vmin=0, vmax=200)


# ## Resting state functional connectome

# In[4]:


# Set the path to the the functional connectome
# acquisitions
func_acq='desc-se_task-rest_acq-AP_bold'
cnt_fs = subjectDir + f'/func/{func_acq}/surf/{subjectID}_surf-fsLR-32k_atlas-{atlas}_desc-FC.shape.gii'

# Load the cortical connectome
mtx_fs = nib.load(cnt_fs).darrays[0].data

# Fill the lower triangle of the matrix
mtx_fcSym = np.triu(mtx_fs,1)+mtx_fs.T

# Plot the matrix
corr_plot = plotting.plot_matrix(mtx_fcSym, figure=(10, 10), labels=None, cmap='Reds', vmin=0, vmax=1)


# ### Time series (ROI x time)

# In[5]:


# Set the path to the the time series file
cnt_time = subjectDir + f'/func/{func_acq}/surf/{subjectID}_surf-fsLR-32k_desc-timeseries_clean.shape.gii'

# Load the time series
mtx_time = nib.load(cnt_time).darrays[0].data

# Plot as a matrix
corr_plot = plotting.plot_matrix(mtx_time, figure=(300, 10), labels=None, cmap='plasma', vmin=-100, vmax=100)


# ## MPC connectomes

# In[6]:


# Set the path to the the MPC cortical connectome
mpc_acq='acq-T1map'
cnt_mpc = subjectDir + f'/mpc/{mpc_acq}/{subjectID}_atlas-{atlas}_desc-MPC.shape.gii'

# Load the cortical connectome
mtx_mpc = nib.load(cnt_mpc).darrays[0].data

# Fill the lower triangle of the matrix
mtx_mpcSym = np.triu(mtx_mpc,1)+mtx_mpc.T

# Plot the matrix
corr_plot = plotting.plot_matrix(mtx_mpcSym, figure=(10, 10), labels=None, cmap='Greens')


# ### Intensity profiles (Profile x ROI)

# In[7]:


# Set the path to the the Intensity profiles file
cnt_int = subjectDir + f'/mpc/{mpc_acq}/{subjectID}_atlas-{atlas}_desc-intensity_profiles.shape.gii'

# Load the Intensity profiles
mtx_int = nib.load(cnt_int).darrays[0].data

# Plot as a matrix
corr_plot = plotting.plot_matrix(mtx_int, figure=(20,10), labels=None, cmap='Greens', colorbar=False)


# ## Geodesic distance connectomes

# In[8]:


# Set the path to the the geodesic distance connectome
cnt_gd = f'{subjectDir}/dist/{subjectID}_atlas-{atlas}_GD.shape.gii'

# Load the cortical connectome
mtx_gd = nib.load(cnt_gd).darrays[0].data

# Plot the matrix
corr_plot = plotting.plot_matrix(mtx_gd, figure=(10, 10), labels=None, cmap='Blues')


# In[ ]:




