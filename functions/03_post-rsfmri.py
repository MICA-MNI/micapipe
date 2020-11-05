import sys
import os
import glob
import numpy as np
import nibabel as nib
from sklearn.linear_model import LinearRegression
import warnings

warnings.simplefilter('ignore')
#import enigmatoolbox.datasets
#from enigmatoolbox.plotting import plot_cortical
#from enigmatoolbox.utils.parcellation import parcel_to_surface, surface_to_parcel

subject = sys.argv[1]
funcDir = sys.argv[2]
labelDir = sys.argv[3]
parcDir = sys.argv[4]

# check if surface directory exist; exit if false
if os.listdir(funcDir+'/surfaces/'):
    print('')
    print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-')
    print('surfaces directory found; lets get the party started!')
    print('-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-')
    print('')
else:
    print('')
    print(':( sad face :( sad face :( sad face :(')
    print('No surface directory. Exiting. Bye-bye')
    print(':( sad face :( sad face :( sad face :(')
    print('')
    exit()

# find necessary files
os.chdir(funcDir+'/surfaces/')
x_lh = " ".join(glob.glob(funcDir+'/surfaces/'+'*lh*c69-32k_10mm*'))
x_rh = " ".join(glob.glob(funcDir+'/surfaces/'+'*rh*c69-32k_10mm*'))
os.chdir(funcDir+'/volumetric/')
x_spike = " ".join(glob.glob(funcDir+'/volumetric/'+'*spikeRegressors_FD*'))
x_dof = " ".join(glob.glob(funcDir+'/volumetric/'+'*singleecho.1D'))
x_refrms = " ".join(glob.glob(funcDir+'/volumetric/'+'*metric_REFRMS.1D'))
x_fd = " ".join(glob.glob(funcDir+'/volumetric/'+'*metric_FD*'))

# Grab subcortical and cerebellar timeseries and merge them to conte69 and native timeseries
sctx = np.loadtxt(funcDir+'/volumetric/'+subject+'_singleecho_timeseries_subcortical.txt')
cereb = np.loadtxt(funcDir+'/volumetric/'+subject+'_singleecho_timeseries_cerebellum.txt')
n = sctx.shape[1] + cereb.shape[1]                 # so we know data.shape[1] - n = num of ctx vertices only
n_sctx = sctx.shape[1]

# We want to clean up native timeseries too!
x_lh_nat = " ".join(glob.glob(funcDir+'/surfaces/'+'*fmri2fs_lh_10mm.mgh'))
x_rh_nat = " ".join(glob.glob(funcDir+'/surfaces/'+'*fmri2fs_rh_10mm.mgh'))

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

data = []
dataNative = []
spike = []

# load and concatenate surface-registered ctx + sctx + cerebellum timeseries
# conte
lh_data = nib.load(x_lh)
lh_data = np.squeeze(lh_data.get_fdata())
rh_data = nib.load(x_rh)
rh_data = np.squeeze(rh_data.get_fdata())
data = np.transpose(np.append(lh_data, rh_data, axis=0))
data = np.append(np.append(data, sctx, axis=1), cereb, axis=1)

# native
lh_data_nat = nib.load(x_lh_nat)
lh_data_nat = np.squeeze(lh_data_nat.get_fdata())
rh_data_nat = nib.load(x_rh_nat)
rh_data_nat = np.squeeze(rh_data_nat.get_fdata())
dataNative = np.transpose(np.append(lh_data_nat, rh_data_nat, axis=0))
dataNative = np.append(np.append(dataNative, sctx, axis=1), cereb, axis=1)

# load confound files
if x_spike:
    spike = np.loadtxt(x_spike)
    
    if spike.ndim == 1:
        spike = np.expand_dims(spike, axis=1)
    else:
        print("spike file loaded as 2D")
    
    # regress out spikes from individual timeseries
    ones = np.ones((spike.shape[0], 1))
    mdl = []
    mdl = np.append(ones, spike, axis=1)

    # conte
    slm = LinearRegression().fit(data, mdl)
    data_corr = data-np.dot(mdl, slm.coef_)

    # native
    slm = LinearRegression().fit(dataNative, mdl)
    dataNative_corr = dataNative-np.dot(mdl, slm.coef_)
    
else:
    del spike
    print('no spikey, no spikey, will skippy')
    
    # conte
    data_corr = data

    # native
    dataNative_corr = dataNative


# save timeseries in conte69 format
np.savetxt(funcDir+'/surfaces/' + subject + '_rsfMRI-timeseries_conte69_clean.txt', data_corr, fmt='%.6f')


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
plt.savefig(funcDir+'/surfaces/' + subject + '_rsfMRI-framewiseDisplacement.png', dpi=300)


# tSNR
lh_nat_noHP = " ".join(glob.glob(funcDir+'/surfaces/'+'*fmri2fs_lh_NoHP.mgh'))
lh_nat_noHP_data = nib.load(lh_nat_noHP)
lh_nat_noHP_data = np.squeeze(lh_nat_noHP_data.get_fdata())
rh_nat_noHP = " ".join(glob.glob(funcDir+'/surfaces/'+'*fmri2fs_rh_NoHP.mgh'))
rh_nat_noHP_data = nib.load(rh_nat_noHP)
rh_nat_noHP_data = np.squeeze(rh_nat_noHP_data.get_fdata())

lhM = np.mean(lh_nat_noHP_data, axis = 1)
lhSD = np.std(lh_nat_noHP_data, axis = 1)
lh_tSNR = np.divide(lhM, lhSD)

rhM = np.mean(rh_nat_noHP_data, axis = 1)
rhSD = np.std(rh_nat_noHP_data, axis = 1)
rh_tSNR = np.divide(rhM, rhSD)

tSNR = np.append(lh_tSNR, rh_tSNR)
tSNR = np.expand_dims(tSNR, axis=1)

np.savetxt(funcDir+'/surfaces/' + subject + '_tSNR.txt', tSNR, fmt='%.12f')


# Parcellate the data to like so many different parcellations, !¡!¡!¡ ôôô-my-god ¡!¡!¡!
# Start with conte parcellations
parcellationList = ['glasser-360_conte69',
                    'vosdewael-100_conte69', 'vosdewael-200_conte69',
                    'vosdewael-300_conte69', 'vosdewael-400_conte69',
                    'schaefer-100_conte69', 'schaefer-200_conte69', 'schaefer-300_conte69',
                    'schaefer-400_conte69', 'schaefer-500_conte69', 'schaefer-600_conte69',
                    'schaefer-700_conte69', 'schaefer-800_conte69', 'schaefer-900_conte69',
                    'schaefer-1000_conte69', 'aparc_conte69']
for parcellation in parcellationList:
    parcPath = os.path.join(parcDir, parcellation) + '.csv'
    thisparc = np.loadtxt(parcPath)

    # QC file generated on schaefer400
    #if parcellation is "schaefer-400_conte69":
        
        # # map the ctx timeseries to the parcels
        # data_corr_ctx = data_corr[:, :-n]
        # ts = surface_to_parcel(data_corr_ctx, parcellation)
        # ts = np.append(ts, data_corr[:, -n:], axis=1)
        # ts_r = np.corrcoef(np.transpose(ts))
        # ts_r[0, :] = 0
        # ts_r[:, 0] = 0
        # np.savetxt(funcDir + '/surfaces/' + subject + '_rsfMRI-connectome_' + parcOutputName + '_conte69_clean.txt',
                   # ts_r)

        # # seed-based correlation and save figure | C O R T I C A L
        # seeds = (19, 57, 82, 98, 120, 137, 152)
        # seed_conn = [None] * len(seeds)
        # for ii in range(len(seeds)):
            # seed_conn[ii] = parcel_to_surface(ts_r[seeds[ii], :-n], parcellation, fill=np.nan)
        # # save cortical figure
        # plot_cortical(array_name=seed_conn, surface_name="conte69", size=(1200, 1200),
                         # cmap=['Purples', 'Blues', 'Greens', 'PuRd', 'YlOrBr', 'Oranges', 'Reds'],
                         # color_bar=True, label_text=['VN', 'SMN', 'DAN', 'SN', 'LMBC', 'FPN', 'DMN'],
                         # transparent_bg=False, screenshot=True, color_range=(0.35, 1), filename=funcDir+'/surfaces/' +
                         # subject + 'rsfMRI-QC_' + parcOutputName + '.png')

        # # seed-based correlation and save figure | S U B C O R T I C A L
        # seeds = range(7)
        # seed_conn = [None] * len(seeds)
        # for ii in seeds:
            # seed_conn[ii] = parcel_to_surface(ts_r[-n + ii, :-n], parcellation, fill=np.nan)
        # # save cortical figure
        # plot_cortical(array_name=seed_conn, surface_name="conte69", size=(1200, 1200),
                      # cmap=['Reds', 'Reds', 'Reds', 'Reds', 'Reds', 'Reds', 'Reds'],
                      # color_bar=True, label_text=['L-', 'L-Amy', 'L-Caud', 'L-Hip', 'L-Pal', 'L-Put', 'L-Thal'],
                      # transparent_bg=False, screenshot=True, color_range=(0, .5), filename=funcDir + '/surfaces/' +
                      # subject + 'rsfMRI-QC_' + parcOutputName + '_left_sctx.png')

        # # seed-based correlation and save figure | C E R E B E L L A R
        # seeds = range(0, 14, 2)
        # seed_conn = [None] * len(seeds)
        # for ii in range(len(seeds)):
            # seed_conn[ii] = parcel_to_surface(ts_r[-n + n_sctx + ii, :-n], parcellation, fill=np.nan)
        # # save cortical figure
        # plot_cortical(array_name=seed_conn, surface_name="conte69", size=(1200, 1200),
                      # cmap=['Reds', 'Reds', 'Reds', 'Reds', 'Reds', 'Reds', 'Reds'],
                      # color_bar=True, label_text=['', '', '', '', '', '', ''],
                      # transparent_bg=False, screenshot=True, color_range=(0, .5), filename=funcDir + '/surfaces/' +
                      # subject + 'rsfMRI-QC_' + parcOutputName + '_cereb.png')

    #else:
    
    # Parcellate cortical timeseries
    data_corr_ctx = data_corr[:, :-n]
    uparcel = np.unique(thisparc)
    ts_ctx = np.zeros([data_corr_ctx.shape[0], len(uparcel)])
    for lab in range(len(uparcel)):
        tmpData = data_corr_ctx[:, thisparc == lab]
        ts_ctx[:,lab] = np.mean(tmpData, axis = 1)
    
    ts = np.append(ts_ctx, data_corr[:, -n:], axis=1)
    ts_r = np.corrcoef(np.transpose(ts))
    ts_r[0, :] = 0
    ts_r[:, 0] = 0
    ts_r = np.triu(ts_r)
    np.savetxt(funcDir + '/surfaces/' + subject + '_rsfMRI-connectome_' + parcellation + '_clean.txt',
               ts_r, fmt='%.6f')


# Now generate native surface connectomes
# These files are saved directly to the freesurfer directory through micapipe postStruct
parcellationList = ['glasser-360',
                    'vosdewael-100', 'vosdewael-200',
                    'vosdewael-300', 'vosdewael-400',
                    'schaefer-100', 'schaefer-200', 'schaefer-300',
                    'schaefer-400', 'schaefer-500', 'schaefer-600',
                    'schaefer-700', 'schaefer-800', 'schaefer-900',
                    'schaefer-1000',
                    'aparc', 'aparc-a2009s',
                    'economo']
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
    native_parc = np.zeros((nativeLength))
    for x in range(len(labels_lh)):
        native_parc[x] = np.where(ctab_lh[:,4] == labels_lh[x])[0][0]
    for x in range(len(labels_rh)):
        native_parc[x + len(labels_lh)] = np.where(ctab_rh[:,4] == labels_rh[x])[0][0] + len(ctab_lh)

    # Generate connectome on native space parcellation
    # Parcellate cortical timeseries
    dataNative_corr_ctx = dataNative_corr[:, :-n]
    uparcel = np.unique(native_parc)
    ts_native_ctx = np.zeros([dataNative_corr_ctx.shape[0], len(uparcel)])
    for lab in range(len(uparcel)):
        tmpData = dataNative_corr_ctx[:, native_parc == lab]
        ts_native_ctx[:,lab] = np.mean(tmpData, axis = 1)

    ts = np.append(ts_native_ctx, dataNative_corr[:, -n:], axis=1)
    np.savetxt(funcDir + '/surfaces/' + subject + '_rsfMRI-timeseries_' + parcellation + '_clean.txt', ts, fmt='%.12f')
    ts_r = np.corrcoef(np.transpose(ts))
    ts_r = np.triu(ts_r)
    np.savetxt(funcDir + '/surfaces/' + subject + '_rsfMRI-connectome_' + parcellation + '_clean.txt', ts_r, fmt='%.6f')
