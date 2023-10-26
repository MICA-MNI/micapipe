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
out_surf_prefix = sys.argv[3]
def arg2float_list(arg):
    """Convert command-line argument into a list of numbers.
    >>> arg2float_list("1,2,3")
    [1, 2, 3]
    """
    return list(map(float, arg.split(',')))
if len(sys.argv)>4:
    depth = arg2float_list(sys.argv[4])
else:
    depth = [0.01,0.02,0.03] # default depths


convergence_threshold = 1e-3
max_iters = int(1e4)

surf = nib.load(in_surf)
V = surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
F = surf.get_arrays_from_intent('NIFTI_INTENT_TRIANGLE')[0].data
laplace = nib.load(in_laplace)
lp = laplace.get_fdata()
print('loaded data and parameters')

for d in depth:

    # apply inverse affine to surface to get to matrix space
    print(laplace.affine)
    V[:,:] = V - laplace.affine[:3,3].T
    for xyz in range(3):
        V[:,xyz] = V[:,xyz]*(1/laplace.affine[xyz,xyz])

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
        pts = lp[V[:,0].astype(int),V[:,1].astype(int),V[:,2].astype(int)]<d
        Vnew[pts,0] += interp_x(V[pts,:])
        Vnew[pts,1] += interp_y(V[pts,:])
        Vnew[pts,2] += interp_z(V[pts,:])
        ssd = np.sum((V-Vnew)**2,axis=None)
        if i%100 == 0:
            print(f'itaration {i}, convergence: {ssd}, still moving: {np.sum(pts)}')
        if ssd < convergence_threshold:
            break
        V[:,:] = Vnew[:,:]

    # return to world coords
    for xyz in range(3):
        V[:,xyz] = V[:,xyz]*(laplace.affine[xyz,xyz])
    V[:,:] = V + laplace.affine[:3,3].T

    nib.save(surf, out_surf_prefix + str(d).replace("0.","") + '.surf.gii')
