#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1. Saves surfaces func.gii with the intensities mapped as a matrix
  for each surface template {fsnative,fsaverage5, fsLR-5k, fsLR-32k}
2. Creates a vertex-wise MPC from fsLR-5k

    Parameters
    ----------
    dataDir : str, Path to the micapipe output directory
    sub     : str, Subject ID
    ses     : str, Session, is SINGLE if no session exist
    acq     : str, Acquisition name of the qMRI
    mpc_dir : str, name of the output directory ["mpc", "mpc-swm"]
    num_surf: int, number of surface to processed

    Usage
    -----
    build_mpc-vertex.py "$out" "$id" "$SES" "${mpc_p}"

@author: rcruces
"""

# Import packages
import sys
import os
import glob
import numpy as np
import nibabel as nb
import scipy

# Get the MICAPIPE path from the environment
micapipe = os.environ["MICAPIPE"]

# Define input arguments
dataDir = sys.argv[1]
sub = sys.argv[2]
ses = sys.argv[3]
acq = sys.argv[4]
mpc_dir = sys.argv[5]
num_surf = sys.argv[6]

def save_gii(data_array, file_name):
    # Initialize gifti: NIFTI_INTENT_SHAPE - 2005, FLOAT32 - 16
    gifti_data = nb.gifti.GiftiDataArray(data=data_array, intent=2005, datatype=16)

    # this is the GiftiImage class
    gifti_img = nb.gifti.GiftiImage(meta=None, darrays=[gifti_data])

    # Save the new GIFTI file
    nb.save(img=gifti_img, filename=file_name)

# Manage single session
if ses=="SINGLE":
    subject_dir=f"{dataDir}/sub-{sub}"
    bids_id=f"sub-{sub}"
else:
    subject_dir=f"{dataDir}/sub-{sub}/{ses}"
    bids_id=f"sub-{sub}_{ses}"

# setting output directory
if acq=="DEFAULT":
    OPATH = f"{subject_dir}/{mpc_dir}/"
else:
    OPATH = f"{subject_dir}/{mpc_dir}/{acq}/"

# Load all the surfaces and create a numpy.array
def get_hemisphere(surface_number, hemi, surf):
    thisname_gii = "{output}{bids_id}_hemi-{hemi}_surf-{surf}_label-MPC-{surface_number:d}.func.gii".format(output=OPATH, bids_id=bids_id, hemi=hemi, surface_number=surface_number+1, surf=surf)
    data = nb.load(thisname_gii).darrays[0].data
    return data

def get_feature_array(surf, Save=True):
    BBl = np.vstack( [get_hemisphere(ii, 'L', surf) for ii in range(int(num_surf))] )
    BBr = np.vstack( [get_hemisphere(ii, 'R', surf) for ii in range(int(num_surf))] )

    # Concatenate hemispheres and flip so pial surface is at the top
    BB = np.flipud(np.concatenate((BBl, BBr), axis = 1))

    if Save==True:
        fileName="{output}{bids_id}_surf-{surf}_desc-intensity_profiles.shape.gii".format(output=OPATH, bids_id=bids_id, surf=surf)
        print('[INFO]... saving '+surf+' intensities as gifti array')
        save_gii(BB, fileName)
    else:
        return(BB)

# Function to build the MPC from an intencity profile
def build_mpc(data, mask):
    # Calculate mean across columns, excluding mask and any excluded labels input
    I_M = np.nanmean(np.float32(np.where(mask, data, np.nan)), axis=1)

    # Get residuals of all columns (controlling for mean)
    I_resid = np.zeros(data.shape)
    for c in range(data.shape[1]):
        y = data[:,c]
        x = I_M
        slope, intercept, _, _, _ = scipy.stats.linregress(x,y)
        y_pred = intercept + slope*x
        I_resid[:,c] = y - y_pred

    # Calculate correlation coefficient of the intesities with residuals
    R = np.corrcoef(I_resid, rowvar=False)

    # Log transform
    MPC = 0.5 * np.log( np.divide(1 + R, 1 - R) )
    MPC[np.isnan(MPC)] = 0
    MPC[np.isinf(MPC)] = 0

    # CLEANUP: correct diagonal and round values to reduce file size
    # Replace all values in diagonal by zeros to account for floating point error
    for i in range(0,MPC.shape[0]):
            MPC[i,i] = 0

    # Output MPC, microstructural profiles, and problem nodes
    return (MPC)

# Load the fsLR-5k mask
mask_lh = nb.load(micapipe + '/surfaces/fsLR-5k.L.mask.shape.gii').darrays[0].data
mask_rh = nb.load(micapipe + '/surfaces/fsLR-5k.R.mask.shape.gii').darrays[0].data
fsLR5k_mask = np.concatenate((mask_lh, mask_rh), axis=0)

# Save the surfaces.array as a func.gii file
surfaces=['fsnative', 'fsaverage5', 'fsLR-5k', 'fsLR-32k']
for x in surfaces: get_feature_array(x)

# Create a vertex-wise MPC from fsLR-5k
surf_array_fsLR5k = get_feature_array('fsLR-5k', Save=False)
MPC_fsLR5k = build_mpc(surf_array_fsLR5k, fsLR5k_mask)
fileName="{output}{bids_id}_surf-fsLR-5k_desc-MPC.shape.gii".format(output=OPATH, bids_id=bids_id)

# Save it as shape GIFTI
print('[INFO]... saving '+fileName)
save_gii(MPC_fsLR5k, fileName)

# cleanup - remove all feature-surf
tmp_files=sorted(glob.glob("{output}/*label-MPC-*.func.gii".format(output=OPATH)))
for x in tmp_files: os.remove(x)
