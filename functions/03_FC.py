#!/usr/bin/env python3.7
# -*- coding: utf-8 -*-

"""
micapipe Functional Connectome and cofounding regression script.

Generates a surface-based clean data [func~spike+wm+csf+gsr] (optional)
Generates multiple functional connectomes based on micapipe's provided parcellations.

    Positional arguments
    ----------

    subject     :  str
                    BIDS id (including ses if necessary, eg. sub-01_ses-01).

    funcDir     :  str
                    Path to subject func directory.

    labelDir    :  str
                    Path to subject surface labels.

    parcDir     :  str
                    Path to micapipe parcellations directory.

    volmDir     :  str
                    Path to subject anat/parc directory (parcellations).

    performNSR  :  int
                    YES: 1, NO: 0. Performs Nuisance Signal Regression

    performGSR  :  int
                    YES: 1, NO: 0. Perform global signal regression to the clean data

    func_lab    :  str
                    Identifier of type of sequence multi/single echo.

    noFC      :  str
                    [TRUE, FALSE]. If True skipps the functional connectomes generation.

Created and modified from 2019 to 2022
@author: A collaborative effort of the MICA lab  :D
"""

import sys
import os
import glob
import numpy as np
import nibabel as nib
from sklearn.linear_model import LinearRegression
import warnings

warnings.simplefilter('ignore')

subject = sys.argv[1]
funcDir = sys.argv[2]
labelDir = sys.argv[3]
parcDir = sys.argv[4]
volmDir = sys.argv[5]
performNSR = sys.argv[6]
performGSR = sys.argv[7]
func_lab = sys.argv[8]
noFC = sys.argv[9]
GSRtag = sys.argv[10]

# Rename outputs with GSR
if GSRtag == 'TRUE':
    gsr='_gsr'
else:
    gsr=''

# check if surface directory exist; exit if false
if os.listdir(funcDir+'/surf/'):
    print('')
    print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-')
    print('surf directory found; lets get the party started!')
    print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-')
    print('')
else:
    print('')
    print(':( sad face :( sad face :( sad face :(')
    print('No surface directory. Exiting. Bye-bye')
    print('Run the module -proc_freesurfer first!')
    print(':( sad face :( sad face :( sad face :(')
    print('')
    exit()


# Load subcortical and cerebellar timeseries
# Subcortex
sctx = np.loadtxt(funcDir+'/volumetric/' + subject + func_lab + '_timeseries_subcortical.txt')
if sctx.shape:
    n_sctx = sctx.shape[1]
else:
    print('Uh oh, your subcortical timeseries file is empty; exiting. Bye-bye')
    exit()

# Cerebellum
# A little hacky because the co-registration to fMRI space make some cerebellar ROIs disappear
# Missing label indices are replaced by zeros in timeseries and FC matrix
cereb_tmp = np.loadtxt(funcDir+'/volumetric/' + subject + func_lab + '_timeseries_cerebellum.txt')
f = open(funcDir+'/volumetric/' + subject + func_lab + '_cerebellum_roi_stats.txt', "rt")
cerebLabels = f.read()
s1 = cerebLabels.find("nii.gz")
startROIs = s1 + len("nii.gz") + 6
values = cerebLabels[startROIs:].split("\t")
roiLabels = values[0::2]

def missing_elements(roiList):
    start, end = 0, 33
    return sorted(set(range(start, end + 1)).difference(roiList))

if len(roiLabels) == 34: # no labels disappeared in the co-registration
    print('All cerebellar labels found in parcellation!')
    cereb = cereb_tmp
    exclude_labels = [np.nan]
else:
    print('Some cerebellar ROIs were lost in co-registration to fMRI space')
    cereb = np.zeros((cereb_tmp.shape[0], 34), dtype=np.int8)
    roiLabelsInt = np.zeros((1,len(roiLabels)), dtype=np.int8)
    for (ii, _) in enumerate(roiLabels):
        roiLabelsInt[0,ii] = int(float(roiLabels[ii]))
    roiLabelsInt = roiLabelsInt - 1
    for (ii, _) in enumerate(roiLabels):
        cereb[:,roiLabelsInt[0,ii]] = cereb_tmp[:,ii]
    exclude_labels = missing_elements(roiLabelsInt[0])
    print('Matrix entries for following ROIs will be zero: ', exclude_labels)

# Calculate number of non cortical rows/colunms in matrix and concatenate
n = sctx.shape[1] + cereb.shape[1] # so we know data.shape[1] - n = num of ctx vertices only

# Load confound files
os.chdir(funcDir+'/volumetric/')
x_spike = " ".join(glob.glob(funcDir+'/volumetric/'+'*spikeRegressors_FD.1D'))
x_dof = " ".join(glob.glob(funcDir+'/volumetric/*'+func_lab+'.1D'))
# x_refmse = " ".join(glob.glob(funcDir+'/volumetric/'+'*metric_REFMSE.1D'))
x_fd = " ".join(glob.glob(funcDir+'/volumetric/'+'*metric_FD*'))
x_csf = " ".join(glob.glob(funcDir+'/volumetric/'+'*CSF*'))
x_wm = " ".join(glob.glob(funcDir+'/volumetric/'+'*WM*'))
x_gs = " ".join(glob.glob(funcDir+'/volumetric/'+'*global*'))


def expand_dim(Data):
    if Data.ndim == 1:
        Data = np.expand_dims(Data, axis=1)
    return Data

# Function that calculates the regressed data
def get_regressed_data(x_spike, Data, performNSR, performGSR, Data_name):
    """ Nuisance signal regression: spikes and (optional) WM/CSF:
        This function loads and generates the regression matrix according to the input parameters
        and returns the data corrected by the selected regresors (wm, csf, gs).
        By default will regress motion parameters and spikeRegressors

    Parameters
    ----------
    x_spike : str (spikeRegressors_FD file)
    Data : array
    performNSR : str (0,1)
    performGSR : str (0,1)

    Return
    ------
    Data_corr : array of data corrected
    """
    dof = np.loadtxt(x_dof)
    csf = expand_dim(np.loadtxt(x_csf))
    wm = expand_dim(np.loadtxt(x_wm))
    gs = expand_dim(np.loadtxt(x_gs))
    spike = []
    print('')
    def check_arrays():
        if np.array_equal(mdl, ones) != True:
            print('apply regression')
            slm = LinearRegression().fit(mdl, Data)
            Data_res = Data - np.dot(mdl, slm.coef_.T)
        else:
            Data_res = Data
        return Data_res
    if x_spike:
        spike = expand_dim(np.loadtxt(x_spike))
        ones = np.ones((spike.shape[0], 1))
        mdl = []

        # regress out spikes from individual timeseries
        if performNSR == "1":
            print(Data_name + ' model : func ~ spikes + dof + wm + csf')
            mdl = np.append(np.append(np.append(np.append(ones, spike, axis=1), dof, axis=1), wm, axis=1), csf, axis=1)
        elif performGSR == "1":
            print(Data_name + ' model : func ~ spikes + dof + wm + csf + gs')
            mdl = np.append(np.append(np.append(np.append(np.append(ones, spike, axis=1), dof, axis=1), wm, axis=1), csf, axis=1), gs, axis=1)
        else:
            print(Data_name + 'Default model : func ~ spikes')
            mdl = np.append(ones, spike, axis=1)
        # apply regression
        Data_corr = check_arrays()
    else:
        ones = np.ones((wm.shape[0], 1))
        print('NO spikeRegressors_FD file, will skip loading: ' + Data_name)
        if performNSR == 1:
            print(Data_name + ', model : func ~ dof + wm + csf')
            mdl = np.append(np.append(np.append(ones, dof, axis=1), wm, axis=1), csf, axis=1)
        elif performGSR == 1:
            print(Data_name + ', model : func ~ dof + wm + csf + gs')
            mdl = np.append(np.append(np.append(np.append(ones, dof, axis=1), wm, axis=1), csf, axis=1), gs, axis = 1)
        else:
            print(Data_name + ', model : none')
            mdl = ones
        # apply regression
        Data_corr = check_arrays()
    return Data_corr

# Process subcortex and cerebellum
sctx_cereb = np.append(sctx, cereb, axis=1)
del sctx
del cereb
sctx_cereb_corr = get_regressed_data(x_spike, sctx_cereb, performNSR, performGSR, 'sctx_cereb')

# ------------------------------------------
#     C O R T E X  processing
# ------------------------------------------

def funcgii_load(gii):
    out = np.zeros([len(gii.darrays),gii.darrays[0].data.shape])
    for n in range(len(gii.darrays)):
        out[n,:] = gii.darrays[n].data.shape
    return out

# Find and load surface-registered cortical timeseries
os.chdir(funcDir+'/surf/')
x_lh = " ".join(glob.glob(funcDir+'/surf/'+'*_hemi-L_surf-fsnative.func.gii'))
x_rh = " ".join(glob.glob(funcDir+'/surf/'+'*_hemi-R_surf-fsnative.func.gii'))
lh_data = funcgii_load(nib.load(x_lh))
rh_data = funcgii_load(nib.load(x_rh))

# exit if more than one scan exists
if len(x_lh.split(" ")) == 1:
    print('')
    print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-')
    print('only one scan found; all good in the hood')
    print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-')
    print('')
else:
    print('more than one scan found; exiting. Bye-bye')
    exit()

# Reformat data
data = []
data = np.transpose(np.append(lh_data, rh_data, axis=0))
n_vertex_ctx = data.shape[1]
del lh_data
del rh_data

# correlation matrices
data_corr = get_regressed_data(x_spike, data, performNSR, performGSR, 'fsLR')

# save spike regressed and concatenanted timeseries (subcortex, cerebellum, cortex)
np.savetxt(funcDir+'/surf/' + subject + '_surf-fsnative_desc-timeseries_clean' + gsr + '.txt', data_corr, fmt='%.6f')

# Read the processed parcellations
parcellationList = os.listdir(volmDir)

# Slice the file names and remove nii*
parcellationList=[sub.split('atlas-')[1].split('.nii')[0] for sub in parcellationList]

# Remove cerebellum and subcortical strings
parcellationList.remove('subcortical')
parcellationList.remove('cerebellum')

if noFC!="TRUE":
    for parcellation in parcellationList:
        # Load left and right annot files
        fname_lh = 'lh.' + parcellation + '_mics.annot'
        ipth_lh = os.path.join(labelDir, fname_lh)
        [labels_lh, ctab_lh, names_lh] = nib.freesurfer.io.read_annot(ipth_lh, orig_ids=True)
        fname_rh = 'rh.' + parcellation + '_mics.annot'
        ipth_rh = os.path.join(labelDir, fname_rh)
        [labels_rh, ctab_rh, names_rh] = nib.freesurfer.io.read_annot(ipth_rh, orig_ids=True)
        # Join hemispheres
        nativeLength = len(labels_lh)+len(labels_rh)
        if dataNative_corr.shape[1] != nativeLength:
            print('ERROR..' + parcellation + '_mics.annot' + ' has a mismatch between the native surface and the annot file!!')
        else:
            native_parc = np.zeros((nativeLength))
            for (x, _) in enumerate(labels_lh):
                native_parc[x] = np.where(ctab_lh[:,4] == labels_lh[x])[0][0]
            for (x, _) in enumerate(labels_rh):
                native_parc[x + len(labels_lh)] = np.where(ctab_rh[:,4] == labels_rh[x])[0][0] + len(ctab_lh)

        # Generate connectome on native space parcellation
        # Parcellate cortical timeseries
        uparcel = np.unique(native_parc)
        ts_native_ctx = np.zeros([data_corr.shape[0], len(uparcel)])
        for (lab, _) in enumerate(uparcel):
            tmpData = data_corr[:, native_parc == int(uparcel[lab])]
            ts_native_ctx[:,lab] = np.mean(tmpData, axis = 1)

        ts = np.append(sctx_cereb_corr, ts_native_ctx, axis=1)
        np.savetxt(funcDir + '/surf/' + subject + '_atlas-' + parcellation + '_desc-timeseries' + gsr + '.txt', ts, fmt='%.12f')

        ts_r = np.corrcoef(np.transpose(ts))

        if np.isnan(exclude_labels[0]) == False:
            for i in exclude_labels:
                ts_r[:, i + n_sctx] = 0
                ts_r[i + n_sctx, :] = 0
            ts_r = np.triu(ts_r)
        else:
            ts_r = np.triu(ts_r)

        np.savetxt(funcDir + '/surf/' + subject + '_atlas-' + parcellation + '_desc-FC' + gsr + '.txt', ts_r, fmt='%.6f')
        del ts_native_ctx
        del native_parc
        del ts_r
        del ts
else:
    print('')
    print('...... no FC was selected, will skipp the functional connectome generation')

# Clean up
del data_corr
del data

# fsLR-5k FC



# ------------------------------------------
# Additional QC
# ------------------------------------------
print('')
print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
print('Calculating tSNR and framewise displacement')
print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
# mean framewise displacement + save plot
fd = np.loadtxt(x_fd)
title = 'mean FD: ' + str(np.mean(fd))
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(16, 6))
plt.plot(fd, color="#2171b5")
plt.title(title, fontsize=16)
ax.set(xlabel='')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(funcDir+'/volumetric/' + subject + func_lab + '_framewiseDisplacement.png', dpi=300)

del fd

# tSNR
lh_nat_noHP = " ".join(glob.glob(funcDir+'/surf/'+'*hemi-L_surf-fsnative_NoHP.func.gii'))
lh_nat_noHP_data = funcgii_load(nib.load(lh_nat_noHP))
rh_nat_noHP = " ".join(glob.glob(funcDir+'/surf/'+'*hemi-R_surf-fsnative_NoHP.func.gii'))
rh_nat_noHP_data = funcgii_load(nib.load(rh_nat_noHP))

# delete NoHP (no longer needed)
os.remove(lh_nat_noHP)
os.remove(rh_nat_noHP)

lhM = np.mean(lh_nat_noHP_data, axis = 1)
lhSD = np.std(lh_nat_noHP_data, axis = 1)
lh_tSNR = np.divide(lhM, lhSD)

rhM = np.mean(rh_nat_noHP_data, axis = 1)
rhSD = np.std(rh_nat_noHP_data, axis = 1)
rh_tSNR = np.divide(rhM, rhSD)

tSNR = np.append(lh_tSNR, rh_tSNR)
tSNR = np.expand_dims(tSNR, axis=1)

np.savetxt(funcDir+'/volumetric/' + subject + func_lab + '_tSNR' + gsr + '.txt', tSNR, fmt='%.12f')
print('')
print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
print('func regression and FC ran successfully')
print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
