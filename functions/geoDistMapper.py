#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute geodesic distance matrix vertex-wise or parcel-wise (centroids) of native surface.
For parcel-wise GD, centroid vertices will be selected for each parcel and distance will only be propogated from these points.

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
        Geodesic distance matrix (.txt) { (Vertex x Vertex) | (nParcel x nParcel) }

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

# get the real paths
lh_surf = args.lh_surf[0]
rh_surf = args.rh_surf[0]

# Print the arguments
print("\nCheck inputs")
print(f'  -lh_surf      : "{lh_surf}"')
print(f'  -rh_surf      : "{rh_surf}"')
print(f'  -outPath      : "{args.outPath}"')
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

# Initialize class to calculate geodesic distance
geoalg_L = geodesic.PyGeodesicAlgorithmExact(vertices_lh,faces_lh)
geoalg_R = geodesic.PyGeodesicAlgorithmExact(vertices_rh,faces_rh)

# Calculate GD parcel-wise
if args.parcel_wise == True:
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
    # Left hemisphere
    parcL = parc[0:len(labels_lh)]
    print("[ INFO ]..... Running geodesic distance on the left hemisphere")
    N = len(np.unique(labels_lh))
    # Define array with left center vertices
    target_vertices = voi[:,0:N][0]
    # Iterate over each central vertex
    for ii in range(N):
        vertex = np.array([ int(voi[0,ii]) ])
        parcGD = geoalg_L.geodesicDistances(vertex, target_vertices)[0]
        GD[ii,:] = np.concatenate((parcGD, np.zeros(N)), axis = 0)

    # Right hemisphere
    print("[ INFO ]..... Running geodesic distance on the right hemisphere")
    N = len(np.unique(labels_rh))
    # Define array with right center vertices
    target_vertices = voi[:,N:][0] - len(vertices_lh)
    # Iterate over each central vertex
    for ii in range(len(np.unique(labels_rh))):
        ii_rh = int(ii + len(uparcel)/2)
        vertex = np.array([ int(voi[0,ii_rh] - len(vertices_lh)) ])
        parcGD = geoalg_R.geodesicDistances(vertex, target_vertices)[0]
        GD[ii_rh,:] = np.concatenate((np.zeros(N), parcGD), axis = 0)
else:
    print("[ INFO ]..... Running geodesic distance vertex-wise")
    N = vertices_lh.shape[0]
    GD = np.empty((N*2,N*2), dtype=float)
    for x in range(0,N,1):
        # Left vertices
        dist_L = geoalg_L.geodesicDistances(np.array([x]))[0]
        GD[x,:] = np.concatenate((dist_L, np.zeros(N)), axis = 0)
        # Right vertices
        dist_R = geoalg_R.geodesicDistances(np.array([x]))[0]
        GD[x+N,:] =  np.concatenate((np.zeros(N), dist_R), axis = 0)

np.savetxt(args.outPath + '.txt', GD, fmt='%.12f')
print("[ INFO ]..... Geodesic distance completed")
