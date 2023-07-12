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

convergence_threshold = 1e-5
max_iters = int(10000*depth)

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
#mask = np.logical_or(np.logical_or(abs(dx>0), abs(dy>0)), abs(dz>0))
#points = np.where(mask)
# make interpolator of gradients
points = (range(lp.shape[0]), range(lp.shape[1]), range(lp.shape[2]))
interp_x = RegularGridInterpolator(points, dx)
interp_y = RegularGridInterpolator(points, dy)
interp_z = RegularGridInterpolator(points, dz)
print('gradient interpolator ready')


from joblib import Parallel, delayed
import warnings

def avg_neighbours(invar):
    '''Averages vertex-wise data at vertex n with its neighbouring vertices. F, cdat, n should be passed as a tuple (for easier parallel).'''
    F,cdat,n = invar
    frows = np.where(F==n)[0]
    v = np.unique(F[frows,:])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = np.nanmean(cdat[v])
    return out
def surfdat_smooth(F,cdata,iters=1,cores=8):
    '''Smoothes surface data across neighbouring vertices. This assumes that vertices are evenly spaced and evenly connected.
    TODO: convert to mm vis calibration curves in https://github.com/MELDProject/meld_classifier/blob/9d3d364de86dc207d3a1e5ec11dcab3ef012ebcb/meld_classifier/mesh_tools.py#L17'''
    cdat = copy.deepcopy(cdata)
    for i in range(iters):
        cdat = Parallel(n_jobs=cores)(delayed(avg_neighbours)((F,cdat,n)) for n in range(len(cdat)))
        cdata_smooth = np.array(cdat)
        cdat = copy.deepcopy(cdata_smooth)
    return cdata_smooth


# shift the surface vertices
for i in range(max_iters):
    Vnew = copy.deepcopy(V)
    pts = lp[V[:,0].astype(int),V[:,1].astype(int),V[:,2].astype(int)]<depth
    Vnew[pts,0] += interp_x(V[pts,:])
    Vnew[pts,1] += interp_y(V[pts,:])
    Vnew[pts,2] += interp_z(V[pts,:])
    # smooth in between iterations
    Vnew[:,0] = surfdat_smooth(F,Vnew[:,0])
    Vnew[:,1] = surfdat_smooth(F,Vnew[:,1])
    Vnew[:,2] = surfdat_smooth(F,Vnew[:,2])
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

