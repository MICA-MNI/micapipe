####################################################################################################

# Translated from matlab: 
# Original script by Boris Bernhardt and modified by Casey Paquola
# Translated to python by Jessica Royer
#
# geoDistMapper   Compute geodesic distance between parcel centroids and all other vertices of
#                 native surface.
#  
# INPUT
# lh_surf         left surface to map along, surfstat or gifti compatible
# rh_surf         left surface to map along, surfstat or gifti compatible
# outPath         output path
# lh_annot        path to left annotation files. Centroid vertices will be selected for
#                 each parcel and distance will only be propogated from these
#                 points.
# rh_annot        path to left annotation files.
#
# OUTPUT
# GD              Geodesic distance matrix, nParcel x nParcel 

####################################################################################################


# Import packages
import sys
import subprocess
import numpy as np
import nibabel as nib
import scipy as sp
from scipy import spatial


# Set up arguments
lh_surf = sys.argv[1]
rh_surf = sys.argv[2]
outPath = sys.argv[3]
lh_annot = sys.argv[4] 
rh_annot = sys.argv[5] 
wbPath = sys.argv[6]

# Load surface
lh = nib.load(lh_surf)
vertices_lh = lh.agg_data('NIFTI_INTENT_POINTSET')
faces_lh = lh.agg_data('NIFTI_INTENT_TRIANGLE')

rh = nib.load(rh_surf)
vertices_rh = rh.agg_data('NIFTI_INTENT_POINTSET')
faces_rh = rh.agg_data('NIFTI_INTENT_TRIANGLE') + len(vertices_lh)

vertices = np.append(vertices_lh, vertices_rh, axis = 0)
faces = np.append(faces_lh, faces_rh, axis = 0)


# Read annotation & join hemispheres
[labels_lh, ctab_lh, names_lh] = nib.freesurfer.io.read_annot(lh_annot, orig_ids=True)
[labels_rh, ctab_rh, names_rh] = nib.freesurfer.io.read_annot(rh_annot, orig_ids=True)
nativeLength = len(labels_lh)+len(labels_rh)
parc = np.zeros((nativeLength))
for x in range(len(labels_lh)):
    parc[x] = np.where(ctab_lh[:,4] == labels_lh[x])[0][0]
for x in range(len(labels_rh)):
    parc[x + len(labels_lh)] = np.where(ctab_rh[:,4] == labels_rh[x])[0][0] + len(ctab_lh)


# Find centre vertex
uparcel = np.unique(parc)
voi = np.zeros([1, len(uparcel)])

print("Finings centre vertex for each parcel")
for n in range(len(uparcel)):
    this_parc = np.where(parc == uparcel[n])[0]
    distances = sp.spatial.distance.pdist(np.squeeze(vertices[this_parc,:]), 'euclidean') # Returns condensed matrix of distances
    distancesSq = sp.spatial.distance.squareform(distances) # convert to square form
    sumDist = np.sum(distancesSq, axis = 1) # sum distance across columns
    index = np.where(sumDist == np.min(sumDist)) # minimum sum distance index
    voi[0, n] = this_parc[index[0][0]]


# Initialize distance matrix
GD = np.empty((uparcel.shape[0], uparcel.shape[0]))


# Left hemisphere
parcL = parc[0:len(labels_lh)]
print("Running geodesic distance in the left hemisphere")
for ii in range(len(np.unique(labels_lh))):
    vertex = int(voi[0,ii])
    voiStr = str(vertex)
    
    cmdStr = "{wbPath} -surface-geodesic-distance {lh_surf} {voiStr} {outPath}_this_voi.func.gii".format(wbPath=wbPath, lh_surf=lh_surf, voiStr=voiStr, outPath=outPath)
    subprocess.run(cmdStr.split())
    
    tmpname = outPath + '_this_voi.func.gii'
    tmp = nib.load(tmpname).agg_data()
    parcGD = np.empty((1, len(np.unique(labels_lh))))
    for n in range(len(np.unique(labels_lh))):
        tmpData = tmp[parcL == uparcel[n]]
        tmpMean = np.mean(tmpData)
        parcGD[0, n] = tmpMean
    
    GD[ii,:] = np.append(parcGD, np.zeros((parcGD.shape[0], parcGD.shape[1])), axis = 1)


# Right hemisphere
parcR = parc[-len(labels_rh):]
print("Running geodesic distance in the right hemisphere")
for ii in range(len(np.unique(labels_rh))):
    ii_rh = int(ii + len(uparcel)/2)
    vertex = int(voi[0,ii_rh] - len(vertices_lh))
    voiStr = str(vertex)
        
    cmdStr = "{wbPath} -surface-geodesic-distance {rh_surf} {voiStr} {outPath}_this_voi.func.gii".format(wbPath=wbPath, rh_surf=rh_surf, voiStr=voiStr, outPath=outPath)
    subprocess.run(cmdStr.split())
    
    tmpname = outPath + '_this_voi.func.gii'
    tmp = nib.load(tmpname).agg_data()
    parcGD = np.empty((1, len(np.unique(labels_rh))))
    for n in range(len(np.unique(labels_rh))):
        n_rh = int(n + len(uparcel)/2)
        tmpData = tmp[parcR == uparcel[n_rh]]
        tmpMean = np.mean(tmpData)
        parcGD[0, n] = tmpMean
    
    GD[ii_rh,:] = np.append(np.zeros((parcGD.shape[0], parcGD.shape[1])), parcGD, axis = 1)

np.savetxt(outPath + '_GD.txt', GD, fmt='%.12f')


