"""
micapipe distance mappper.

Generates a surface-based/atlas-base based ED or GD matrices.
ED = Euclidena Distance Matrix
GD = Gerodesic distance Matrix

    Positional arguments
    ----------

    inFile   :   mtx = txt input matrix of the mass center data to calculate the EDM
                     surf = surf.gii file used to calculate the GD

    outName      :  str
                    Path to the output file.

    feature      :  str
                    Algrith to calculate the distance {EDM, GM}

    str_feat     :  str
                    Identifier of the file to be process

@author: MICA lab
"""

import sys
import os
import numpy as np
import warnings
import pygeodesic
import pygeodesic.geodesic as geodesic
from brainspace.mesh.mesh_io import read_surface
import nibabel as nb
import vtk

warnings.simplefilter('ignore')

inFile = sys.argv[1]
outName = sys.argv[2]
dist_feat = sys.argv[3]
str_feat = sys.argv[4]

def getPointsAndCellsFromPolydata(polydata):
    """
    Extract points and cells from polydata object
    Args:
        polydata (vtk.vtkPolyData): polydata object representing mesh

    Returns:
        ndarray: Array of points
        ndarrdy: Array of cells
    """

    # Use numpy wrapper to easily extract points and cells
    polydata = dsa.WrapDataObject(polydata) if isinstance(polydata, vtk.vtkPolyData) else polydata

    # Get points
    points = np.array(polydata.GetPoints(), dtype=np.float64)

    # Get cells
    polygons = np.array(polydata.GetPolygons(), dtype=np.int32)
    n = polygons[0]+1
    polygons = np.resize(polygons, (polygons.size//n, n))
    polygons = polygons[:, 1:n]

    return points, polygons

print('[INFO].... Distance calculation')
# ------------------------------------------
# Calculate GD
# ------------------------------------------
if dist_feat == 'GD':
    print('[INFO].... Calculating GD of ' + str_feat)

    print('Load surface')
    GD_surf = read_surface(inFile, itype='gii')

    print('Calculate the GD: Initialise the PyGeodesicAlgorithmExact class instance')
    points, faces = getPointsAndCellsFromPolydata(GD_surf)
    geoalg = geodesic.PyGeodesicAlgorithmExact(points,faces)

    N = GD_surf.points.shape[0]
    GD = np.empty((N,N), dtype=float)
    for x in range(0,N,1):
        distances, best_source = geoalg.geodesicDistances(np.array([x]))
        GD[x,:] = distances

    print('Save the GD file')
    np.savetxt(outName, GD)

elif dist_feat == 'EDM':
    print('[INFO].... Calculating EDM of' + str_feat)

    # Read distance table
    coord = np.loadtxt(inFile,delimiter=",")

    # Calculate distance matri
    ED = np.array([ np.linalg.norm(coord - p, axis=1) for p in coord])

    # Save file
    np.save(outName, ED)
