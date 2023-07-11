# code from https://github.com/khanlab/hippunfold/blob/master/hippunfold/workflow/scripts/create_warps.py

# shifts a wm surface inward along a Laplace field


import nibabel as nib
import numpy as np
from scipy.interpolate import griddata

in_surf = sys.argv[1]
in_laplace = sys.argv[2]
out_surf = sys.argv[3]
depth = sys.argv[4]

in_surf = '/data/mica3/BIDS_PNI/derivatives/micapipe_v0.2.0/sub-PNC003/ses-01/surf/sub-PNC003_ses-01_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii'
surf = nib.load(in_surf)
V = surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
in_laplace='tmp.nii.gz'
laplace = nib.load(in_laplace)

# laplace to gradient
dx,dy,dz = np.gradient(lp)
mask = np.logical_or(np.logical_or(abs(dx>0), abs(dy>0)), abs(dz>0))
i,j,k = np.where(mask)
# matrix space to world coordinates
realworldcoords = np.stack((i,j,k,np.ones((len(i)))))
realworldcoords = laplace.affine @ realworldcoords
realworldcoords = realworldcoords[:3,:]
# build interpolant to sample gradient at any given world point
interp_x = LinearNDInterpolator(realworldcoords.T,dx[i,j,k])
interp_y = LinearNDInterpolator(realworldcoords.T,dy[i,j,k])
interp_z = LinearNDInterpolator(realworldcoords.T,dz[i,j,k])

# shift the surface vertices
for i in range(depth):
    for v in range(len(V)):
        V[v,0] += interp_x[V[v,:]]
        V[v,1] += interp_y[V[v,:]]
        V[v,2] += interp_z[V[v,:]]

nib.save(surf, out_surf)

