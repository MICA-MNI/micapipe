#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute geodesic distance matrix vertex-wise or parcel-wise (centroids) of native surface.
For parcel-wise GD, centroid vertices will be selected for each parcel and distance will only be propagated from these points.

    Parameters
    ----------

    lh_surf : str
        Path to left surface (GIFTI).

    rh_surf : str
        Path to right surface (GIFTI).

    outPath : str
        Output GD matrix path and file name.

    lh_annot : str, mandatory with -parcel_wise
        Path to left annotation files in the same surface as lh_surf.

    rh_annot: str, mandatory with -parcel_wise
        Path to right annotation files in the same surface as rh_surf.

    parcel_wise  : optional
        Calculate distance parcel-wise. Default is vertex-wise.

    Returns
    -------
    GD  : numpy.ndarray file
        Geodesic distance matrix (.shape.gii) { (Vertex x Vertex) | (nParcel x nParcel) }

    Usage
    -----
    geoDistMapper   Compute geodesic distance between parcel centroids and all other vertices of

Original matlab script by Boris Bernhardt and modified by Casey Paquola
Translated to python by Jessica Royer
Updated by rcruces march 2023
"""

# Import packages
import os
import sys
import argparse
import subprocess
import numpy as np
import nibabel as nb
from scipy import spatial
import pygeodesic.geodesic as geodesic

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-lh_surf',
                    dest='lh_surf',
                    type=str,
                    help='Path to left surface (GIFTI).',
                    nargs='+'
                    )
parser.add_argument('-rh_surf',
                    dest='rh_surf',
                    type=str,
                    help='Path to right surface (GIFTI).',
                    nargs='+'
                    )
parser.add_argument('-outPath',
                    dest='outPath',
                    type=str,
                    help='Output GD matrix path and file name'
                    )
parser.add_argument('-parcel_wise',
                    default=False,
                    action='store_true',
                    help='Calculate distance parcel-wise'
                    )
parser.add_argument('-lh_annot',
                    dest='lh_annot',
                    type=str,
                    help='Path to left annotation files in the same surface as lh_surf.'
                    )
parser.add_argument('-rh_annot',
                    dest='rh_annot',
                    type=str,
                    help='Path to right annotation files in the same surface as rh_surf.'
                    )
args = parser.parse_args()

# Function save as gifti
def save_gii(data_array, file_name):
    # Initialize gifti: NIFTI_INTENT_SHAPE - 2005, FLOAT32 - 16
    gifti_data = nb.gifti.GiftiDataArray(data=data_array, intent=2005, datatype=16)

    # this is the GiftiImage class
    gifti_img = nb.gifti.GiftiImage(meta=None, darrays=[gifti_data])

    # Save the new GIFTI file
    nb.save(img=gifti_img, filename=file_name)

# get the real paths
lh_surf = args.lh_surf[0]
rh_surf = args.rh_surf[0]
outPath = args.outPath
# Print the arguments
print("\nCheck inputs")
print(f'  -lh_surf      : "{lh_surf}"')
print(f'  -rh_surf      : "{rh_surf}"')
print(f'  -outPath      : "{outPath}"')
print(f'  -parcel_wise  : "{args.parcel_wise}"')
if args.parcel_wise==True:
    lh_annot = args.lh_annot
    rh_annot = args.rh_annot
    print(f'  -lh_annot   : "{lh_annot}"')
    print(f'  -rh_annot :  "{rh_annot}"')
# ----------------------------------------------------

# Load LEFT surfaces
lh = nb.load(lh_surf)
vertices_lh = lh.agg_data('NIFTI_INTENT_POINTSET')
faces_lh = lh.agg_data('NIFTI_INTENT_TRIANGLE')
# Load RIGHT surfaces
rh = nb.load(rh_surf)
vertices_rh = rh.agg_data('NIFTI_INTENT_POINTSET')
faces_rh = rh.agg_data('NIFTI_INTENT_TRIANGLE')

# Calculate GD parcel-wise
if args.parcel_wise == True:
    # WORKBENCH COMMAND PATH
    wbPath = os.popen("which wb_command").read()
    wbPath = wbPath.replace('wb_command\n', 'wb_command')

    vertices = np.append(vertices_lh, vertices_rh, axis = 0)
    faces = np.append(faces_lh, faces_rh+len(vertices_lh), axis = 0)

    # Read annotation & join hemispheres
    [labels_lh, ctab_lh, names_lh] = nb.freesurfer.io.read_annot(lh_annot, orig_ids=True)
    [labels_rh, ctab_rh, names_rh] = nb.freesurfer.io.read_annot(rh_annot, orig_ids=True)
    nativeLength = len(labels_lh)+len(labels_rh)
    parc = np.zeros((nativeLength))
    for (x, _) in enumerate(labels_lh):
        parc[x] = np.where(ctab_lh[:,4] == labels_lh[x])[0][0]
    for (x, _) in enumerate(labels_rh):
        parc[x + len(labels_lh)] = np.where(ctab_rh[:,4] == labels_rh[x])[0][0] + len(ctab_lh)

    # Find centre vertex
    uparcel = np.unique(parc)
    voi = np.zeros([1, len(uparcel)])

    print("[ INFO ]..... Finding central vertex for each parcel")
    for (n, _) in enumerate(uparcel):
        this_parc = np.where(parc == uparcel[n])[0]
        distances = spatial.distance.pdist(np.squeeze(vertices[this_parc,:]), 'euclidean') # Returns condensed matrix of distances
        distancesSq = spatial.distance.squareform(distances) # convert to square form
        sumDist = np.sum(distancesSq, axis = 1) # sum distance across columns
        index = np.where(sumDist == np.min(sumDist)) # minimum sum distance index
        voi[0, n] = this_parc[index[0][0]]

    # Initialize distance matrix
    GD = np.empty((uparcel.shape[0], uparcel.shape[0]))

    # ---------------------------------------------------------------------
    # Calculate distance from VERTEX to all other central VERTICES
    cpus = os.environ['OMP_NUM_THREADS']
    # Left hemisphere
    parcL = parc[0:len(labels_lh)]
    print("[ INFO ]..... Running geodesic distance using {cpus} threads".format(cpus=cpus))
    N = len(np.unique(labels_lh))
    # Iterate over each central vertex
    for ii in range(N):
        vertex = int(voi[0,ii])
        voiStr = str(vertex)
        # Calculate the GD from tha vertex to the rest of the vertices
        cmdStr = "{wbPath} -surface-geodesic-distance {lh_surf} {voiStr} {outPath}_this_voi.func.gii".format(wbPath=wbPath, lh_surf=lh_surf, voiStr=voiStr, outPath=outPath)
        subprocess.run(cmdStr.split())
        # Load file with GD column format
        tmpname = outPath + '_this_voi.func.gii'
        tmp = nb.load(tmpname).agg_data()
        parcGD = np.empty((1, N))
        for n in range(N):
            tmpData = tmp[parcL == uparcel[n]]
            tmpMean = np.mean(tmpData)
            parcGD[0, n] = tmpMean
        # Save column on GD matrix
        GD[ii,:] = np.append(parcGD, np.zeros((parcGD.shape[0], parcGD.shape[1])), axis = 1)

    # Right hemisphere
    N = len(np.unique(labels_rh))
    parcR = parc[-len(labels_rh):]
    # Iterate over each central vertex
    for ii in range(N):
        ii_rh = int(ii + len(uparcel)/2)
        vertex = int(voi[0,ii_rh] - len(vertices_lh))
        voiStr = str(vertex)
        # Calculate the GD from tha vertex to the rest of the vertices
        cmdStr = "{wbPath} -surface-geodesic-distance {rh_surf} {voiStr} {outPath}_this_voi.func.gii".format(wbPath=wbPath, rh_surf=rh_surf, voiStr=voiStr, outPath=outPath)
        subprocess.run(cmdStr.split())
        # Load file with GD column format
        tmpname = outPath + '_this_voi.func.gii'
        tmp = nb.load(tmpname).agg_data()
        parcGD = np.empty((1, N))
        for n in range(N):
            n_rh = int(n + len(uparcel)/2)
            tmpData = tmp[parcR == uparcel[n_rh]]
            tmpMean = np.mean(tmpData)
            parcGD[0, n] = tmpMean
        # Save column on GD matrix
        GD[ii_rh,:] = np.append(np.zeros((parcGD.shape[0], parcGD.shape[1])), parcGD, axis = 1)
else:
    print("[ INFO ]..... Running geodesic distance vertex-wise")
    # Initialize class to calculate geodesic distance
    geoalg_L = geodesic.PyGeodesicAlgorithmExact(vertices_lh,faces_lh)
    geoalg_R = geodesic.PyGeodesicAlgorithmExact(vertices_rh,faces_rh)

    N = vertices_lh.shape[0]
    GD = np.empty((N*2,N*2), dtype=float)
    for x in range(0,N,1):
        # Left vertices
        dist_L = geoalg_L.geodesicDistances(np.array([x]))[0]
        GD[x,:] = np.concatenate((dist_L, np.zeros(N)), axis = 0)
        # Right vertices
        dist_R = geoalg_R.geodesicDistances(np.array([x]))[0]
        GD[x+N,:] =  np.concatenate((np.zeros(N), dist_R), axis = 0)

save_gii(GD, outPath+'.shape.gii')
print("[ INFO ]..... Geodesic distance completed")
