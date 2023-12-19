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
from build_mpc import build_mpc

# Define input arguments
dataDir = sys.argv[1]
sub = sys.argv[2]
ses_num = sys.argv[3]
acq = sys.argv[4]

def save_gii(data_array, file_name):
    # Initialize gifti: NIFTI_INTENT_SHAPE - 2005, FLOAT32 - 16
    gifti_data = nb.gifti.GiftiDataArray(data=data_array, intent=2005, datatype=16)

    # this is the GiftiImage class
    gifti_img = nb.gifti.GiftiImage(meta=None, darrays=[gifti_data])

    # Save the new GIFTI file
    nb.save(img=gifti_img, filename=file_name)

# Number of surfaces (Harcoded since the begining)
num_surf = 14

# Manage single session
if ses_num=="SINGLE":
    ses_str="{dataDir}/sub-{sub}".format(dataDir=dataDir, sub=sub)
    bids_id="sub-{sub}".format(sub=sub)
else:
    ses_str="{dataDir}/sub-{sub}/{ses}".format(dataDir=dataDir, sub=sub, ses=ses_num)
    bids_id="sub-{sub}_{ses}".format(sub=sub, ses=ses_num)

# setting output directory
if acq=="DEFAULT":
    OPATH = "{subject_dir}/mpc/".format(subject_dir=ses_str)
else:
    OPATH = "{subject_dir}/mpc/{acq}/".format(subject_dir=ses_str, acq=acq)

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

# Save the surfaces.array as a func.gii file
surfaces=['fsnative', 'fsaverage5', 'fsLR-5k', 'fsLR-32k']
for x in surfaces: get_feature_array(x)

# Create a vertex-wise MPC from fsLR-5k
surf_array_fsLR5k = get_feature_array('fsLR-5k', Save=False)
(MPC_fsLR5k, I, problemNodes) = build_mpc(surf_array_fsLR5k)
fileName="{output}{bids_id}_surf-fsLR-5k_desc-MPC.shape.gii".format(output=OPATH, bids_id=bids_id)

# Save it as shape GIFTI
print('[INFO]... saving '+fileName)
save_gii(MPC_fsLR5k, fileName)

# cleanup - remove all feature-surf
tmp_files=sorted(glob.glob("{output}/*label-MPC-*.func.gii".format(output=OPATH)))
for x in tmp_files: os.remove(x)
