# code from https://github.com/khanlab/hippunfold/blob/master/hippunfold/workflow/scripts/laplace_coords.py

# solves Laplace equation over the domain of white matter, using grey matter as the source and ventricles as the sink. inputs are expected to be Free/FastSurfer aparc+aseg.mgz or .nii.gz files


import nibabel as nib
import numpy as np
import skfmm
from scipy.ndimage import binary_dilation, convolve
import sys
import os
import glob

print('starting laplace solver')
in_seg = sys.argv[1]
out_laplace = sys.argv[2]

# parameters
convergence_threshold = 1e-4
max_iters = 10000
#fg_labels = [2, 4, 11, 12, 26, 17, 31, 10, 5, 28, 13, 30, 41, 43, 50, 51, 58, 53, 63, 49, 44, 60, 52, 62, 77, 255, 254, 253, 252, 251, 72, 80]
fg_labels = [41, 2]
#src_labels = np.hstack(([54, 18], np.arange(1000,2999)))
src_labels = np.arange(1000,2999)
# sink_labels = [4, 43, 31, 63]#, 5, 44]

lbl_nib = nib.load(in_seg)
lbl = lbl_nib.get_fdata()
print('loaded data and parameters')

# initialize foreground , source, and sink
fg = np.isin(lbl,fg_labels)
fg = binary_dilation(fg) # dilate to make sure we always "catch" neighbouring surfaces in our gradient
source = np.isin(lbl,src_labels)
source[fg] = 0
#sink = np.isin(lbl,sink_labels)
sink = 1-fg-source

# initialize solution with fast marching
# fast march forward
phi = np.ones_like(lbl)
phi[source == 1] = 0
mask = np.ones_like(lbl)
mask[fg == 1] = 0
mask[source == 1] = 0
phi = np.ma.MaskedArray(phi, mask)
forward = skfmm.travel_time(phi, np.ones_like(lbl))
forward = forward.data
# fast match backward
phi = np.ones_like(lbl)
phi[sink == 1] = 0
mask = np.ones_like(lbl)
mask[fg == 1] = 0
mask[sink == 1] = 0
phi = np.ma.MaskedArray(phi, mask)
backward = skfmm.travel_time(phi, np.ones_like(lbl))
backward = backward.data
# combine
forward = forward / np.max(forward)
backward = backward / np.max(backward)
backward = -backward + 1
init_coords = (forward + backward) / 2
init_coords = init_coords-np.min(init_coords)
init_coords = init_coords / np.max(init_coords)
init_coords[fg == 0] = 0

# set up filter (27NN)
hl = np.ones([3, 3, 3])
hl = hl / np.sum(hl)

# initialize coords
coords = init_coords
coords[source == 1] = 0
coords[sink == 1] = 1

print('initialized solution')

upd_coords = coords.copy()

# iterate until the solution doesn't change anymore (or reach max iters)
for i in range(max_iters):

    upd_coords = nan_convolve(coords, hl, fill_value=np.nan, preserve_nan=True)

    upd_coords[source == 1] = 0
    upd_coords[sink == 1] = 1

    # check difference between last
    diff_coords = coords - upd_coords
    diff_coords[np.isnan(diff_coords)] = 0
    ssd = (diff_coords * diff_coords).sum(axis=None)
    print(f'itaration {i}, convergence: {ssd}')
    if ssd < convergence_threshold:
        break
    coords = upd_coords


# save file
print('saving')
coords_nib = nib.Nifti1Image(coords, lbl_nib.affine, lbl_nib.header)
nib.save(coords_nib, out_laplace)
