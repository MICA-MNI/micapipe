# code from https://github.com/khanlab/hippunfold/blob/master/hippunfold/workflow/scripts/create_warps.py

# shifts a wm surface inward along a Laplace field

import copy
import nibabel as nib
import numpy as np
from scipy.interpolate import NearestNDInterpolator
import sys

print('starting surface shift')

in_surf = sys.argv[1]
in_laplace = sys.argv[2]
out_surf_prefix = sys.argv[3]
def arg2float_list(arg):
    return list(map(float, arg.split(',')))
if len(sys.argv)>4:
    depth = arg2float_list(sys.argv[4])
else:
    depth = [1,2,3] # default depths

convergence_threshold = 1e-4
step_size = 10.0
max_iters = int(1e4)


# load data
surf = nib.load(in_surf)
V = surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
F = surf.get_arrays_from_intent('NIFTI_INTENT_TRIANGLE')[0].data
laplace = nib.load(in_laplace)
lp = laplace.get_fdata()
print('loaded data and parameters')

# laplace to gradient
dx,dy,dz = np.gradient(lp*step_size)

# apply affine to move from voxels to world
mask = np.logical_or(dx!=0, dy!=0, dz!=0)
x,y,z = np.where(mask)
points = np.vstack((x,y,z,np.ones((len(x)))))
points = laplace.affine @ points
points = points[:3,:].T

# make interpolator of gradients
# Note Linear is better, but very slow!
interp_x = NearestNDInterpolator(points, dx[mask])
interp_y = NearestNDInterpolator(points, dy[mask])
interp_z = NearestNDInterpolator(points, dz[mask])
print('gradient interpolator ready')

distance_travelled = np.zeros((len(V)))
for d in depth:
    # shift the surface vertices
    for i in range(max_iters):
        Vnew = copy.deepcopy(V)
        pts = distance_travelled < d
        stepx = interp_x(V[pts,:])
        stepy = interp_y(V[pts,:])
        stepz = interp_z(V[pts,:])
        Vnew[pts,0] += stepx
        Vnew[pts,1] += stepy
        Vnew[pts,2] += stepz
        distance_travelled[pts] += np.sqrt(stepx**2 + stepy**2 + stepz**2)
        ssd = np.sum((V-Vnew)**2,axis=None)
        if i%100 == 0:
            print(f'itaration {i}, convergence: {ssd}, still moving: {np.sum(pts)}')
        if ssd < convergence_threshold:
            break
        V[:,:] = Vnew[:,:]
    nib.save(surf, out_surf_prefix + str(d) + 'mm.surf.gii')
