# code from https://github.com/khanlab/hippunfold/blob/master/hippunfold/workflow/scripts/create_warps.py

# shifts a wm surface inward along a Laplace field

import copy
import nibabel as nib
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import sys

print('starting surface shift')

in_surf = sys.argv[1]
in_laplace = sys.argv[2]
out_surf = sys.argv[3]
depth = float(sys.argv[4])

convergence_threshold = 1e-3
max_iters = int(100000*depth)

surf = nib.load(in_surf)
V = surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
F = surf.get_arrays_from_intent('NIFTI_INTENT_TRIANGLE')[0].data
laplace = nib.load(in_laplace)
lp = laplace.get_fdata()
print('loaded data and parameters')

# apply inverse affine to surface to get to matrix space
V[:,:] = V - laplace.affine[:3,3].T
for d in range(3):
    V[:,d] = V[:,d]*(1/laplace.affine[d,d])
# laplace to gradient
dx,dy,dz = np.gradient(lp)
# make interpolator of gradients
points = (range(lp.shape[0]), range(lp.shape[1]), range(lp.shape[2]))
interp_x = RegularGridInterpolator(points, dx)
interp_y = RegularGridInterpolator(points, dy)
interp_z = RegularGridInterpolator(points, dz)
print('gradient interpolator ready')


# shift the surface vertices
for i in range(max_iters):
    Vnew = copy.deepcopy(V)
    pts = lp[V[:,0].astype(int),V[:,1].astype(int),V[:,2].astype(int)]<depth
    Vnew[pts,0] += interp_x(V[pts,:])
    Vnew[pts,1] += interp_y(V[pts,:])
    Vnew[pts,2] += interp_z(V[pts,:])
    ssd = np.sum((V-Vnew)**2,axis=None)
    print(f'itaration {i}, convergence: {ssd}, still moving: {np.sum(pts)}')
    if ssd < convergence_threshold:
        break
    V[:,:] = Vnew[:,:]

# return to world coords
for d in range(3):
    V[:,d] = V[:,d]*(laplace.affine[d,d])
V[:,:] = V + laplace.affine[:3,3].T

nib.save(surf, out_surf)

