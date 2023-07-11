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

dx,dy,dz = np.gradient(laplace)

for i in range(depth):
    V = V + gradient[V[:,0],V[:,1],V[:,2]]

nib.save(surf, out_surf)

