#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
micapipe surface QC.

Generates PNG surface images for Quality Check in micapipe.

    Parameters
    ----------

    subBIDS  :  str
                BIDS id (including ses if necessary, eg. sub-01_ses-01).

    subDir   :  str
                Path to subject derivatives directory.

    Examples
    --------
    >>> qc_surf.py -subBIDS sub-00_ses-01 -subDir ~/derivatives/micapipe/sub-00/ses-01

Created on Tue July 9 2021 (the second year of the pademic).
@author: rcruces.
"""

# packages
import glob
import os
import nibabel as nb
import argparse
from nibabel.freesurfer.mghformat import load
import numpy as np
import matplotlib as plt
import matplotlib.pyplot as pltpy
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.mesh.mesh_operations import combine_surfaces
from brainspace.vtk_interface import wrap_vtk, serial_connect
from brainspace.datasets import load_conte69
from vtk import vtkPolyDataNormals

def load_surface(lh, rh, with_normals=True, join=False):
    """
    Loads surfaces.

    Parameters
    ----------
    with_normals : bool, optional
        Whether to compute surface normals. Default is True.
    join : bool, optional.
        If False, return one surface for left and right hemispheres. Otherwise,
        return a single surface as a combination of both left and right.
        surfaces. Default is False.

    Returns
    -------
    surf : tuple of BSPolyData or BSPolyData.
        Surfaces for left and right hemispheres. If ``join == True``, one
        surface with both hemispheres.
    """

    surfs = [None] * 2
    for i, side in enumerate([lh, rh]):
        surfs[i] = read_surface(side)
        if with_normals:
            nf = wrap_vtk(vtkPolyDataNormals, splitting=False,
                          featureAngle=0.1)
            surfs[i] = serial_connect(surfs[i], nf)

    if join:
        return combine_surfaces(*surfs)
    return surfs[0], surfs[1]

def cmap_gradient(N, base_cmaps=['inferno', 'Dark2', 'Set1', 'Set2']):
    """
    Creates a gradient color map of a defined lenght.

    Parameters
    ----------
    N     : int
         Number of colors to extract from each of the base_cmaps below.
    cmaps : list
        List of color map(s) to merge.

    Returns
    -------
    cmap : numpy ndarray
        Surfaces for left and right hemispheres. If ``join == True``, one
        surface with both hemispheres.
    """
    # number of colors
    N = 75

    # we go from 0.2 to 0.8 below to avoid having several whites and blacks in the resulting cmaps
    colors = np.concatenate([pltpy.get_cmap(name)(np.linspace(0,1,N)) for name in base_cmaps])
    cmap = plt.colors.ListedColormap(colors)
    return cmap

# Arguments
parser = argparse.ArgumentParser()

parser.add_argument('-subBIDS',
                    dest='subBIDS',
                    type=str,
                    help='BIDS id (including ses if necessary, eg. sub-01_ses-01)'
                    )

parser.add_argument('-subDir',
                    dest='subDir',
                    type=str,
                    help='Path to subject derivatives directory'
                    )

args = parser.parse_args()

# Arguments
subBIDS = args.subBIDS
subDir = os.path.realpath(args.subDir)
derivatives = subDir.split('/micapipe')[0]

# Print the arguments
print("\nCheck inputs:")
print(f'  -subBIDS  : "{subBIDS}"')
print(f'  -subDir   : "{subDir}"')

# Set paths and variables
dir_fS = derivatives + '/freesurfer/'
dir_fs_label = dir_fS + subBIDS + '/label/'
dir_conte = subDir + '/anat/surfaces/conte69/'
dir_morph = subDir + '/anat/surfaces/morphology/'
dir_mpc = subDir + '/anat/surfaces/micro_profiles/'
dir_QC_png = subDir + '/QC/png/'

# Colormap
ColCurv= plt.colors.ListedColormap(['#A2CD5A', '#A0CA5B', '#9FC85C', '#9EC55D', '#9DC35E', '#9CC05F', '#9BBE61', '#9ABB62', '#99B963', '#98B664', '#96B465', '#95B166', '#94AF68', '#93AC69', '#92AA6A', '#91A76B', '#90A56C', '#8FA26D', '#8EA06F', '#8C9D70', '#8B9B71', '#8A9972', '#899673', '#889475', '#879176', '#868F77', '#858C78', '#848A79', '#82877A', '#81857C', '#80827D', '#7F807E', '#807D7D', '#827A7A', '#857777', '#877575', '#8A7272', '#8C6F6F', '#8F6C6C', '#916969', '#946666', '#966464', '#996161', '#9B5E5E', '#9D5B5B', '#A05858', '#A25656', '#A55353', '#A75050', '#AA4D4D', '#AC4A4A', '#AF4747', '#B14545', '#B44242', '#B63F3F', '#B93C3C', '#BB3939', '#BE3636', '#C03434', '#C33131', '#C52E2E', '#C82B2B', '#CA2828', '#CD2626'])

# Load fsaverage5 inflated
fs5I_lh = read_surface(dir_fS+'fsaverage5/surf/lh.inflated', itype='fs')
fs5I_rh = read_surface(dir_fS+'fsaverage5/surf/rh.inflated', itype='fs')

# Load conte69
c69_lh, c69_rh = load_conte69()

# Load native mid surface
mid_lh, mid_rh = load_surface(dir_fS+subBIDS+'/surf/lh.midthickness.surf.gii',
                                dir_fS+subBIDS+'/surf/rh.midthickness.surf.gii', with_normals=True, join=False)

# Load native surface
surf_lh = read_surface(dir_fS+subBIDS+'/surf/lh.pial', itype='fs')
surf_rh = read_surface(dir_fS+subBIDS+'/surf/rh.pial', itype='fs')

# Load native white matter surface
wm_lh = read_surface(dir_fS+subBIDS+'/surf/lh.white', itype='fs')
wm_rh = read_surface(dir_fS+subBIDS+'/surf/rh.white', itype='fs')

# Load native inflated surface
inf_lh = read_surface(dir_fS+subBIDS+'/surf/lh.inflated', itype='fs')
inf_rh = read_surface(dir_fS+subBIDS+'/surf/rh.inflated', itype='fs')

# ------------------------------------------------------------------------------ #
# proc_freesurfer
try:
    # Freesurfer native thickness
    th = np.concatenate((nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/lh.thickness'), nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/rh.thickness')), axis=0)
    plot_hemispheres(surf_lh, surf_rh, array_name=th, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno",transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_thickness.png')

    # Freesurfer native curvature
    cv = np.concatenate((nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/lh.curv'), nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/rh.curv')), axis=0)
    plot_hemispheres(wm_lh, wm_rh, array_name=cv, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv,transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_curv.png')

    # Freesurfer native sulcal depth
    sd = np.concatenate((nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/lh.sulc'), nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/rh.sulc')), axis=0)
    plot_hemispheres(wm_lh, wm_rh, array_name=sd, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(-5, 5), cmap='cividis',transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_sulc.png')

    # FS Destrieux atlas (aparc.a2009s)
    parc = np.concatenate((nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/lh.aparc.a2009s.annot')[0], nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/rh.aparc.a2009s.annot')[0]), axis=0)
    plot_hemispheres(surf_lh, surf_rh, array_name=parc, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), cmap=cmap_gradient(len(np.unique(parc)), ['inferno', 'hsv', 'hsv', 'tab20b']),transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_a2009s.png')

    # Desikan-Killiany Atlas (aparc)
    parcDK = np.concatenate((nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/lh.aparc.annot')[0], nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/rh.aparc.annot')[0]), axis=0)
    plot_hemispheres(surf_lh, surf_rh, array_name=parcDK, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), cmap=cmap_gradient(len(np.unique(parcDK)), ['inferno', 'hsv', 'hsv', 'tab20b']), transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_aparc.png')
except ValueError:
    print("[WARNING].... some freesurfer processing is missing")
else:
    print("[INFO].... Creating PNGs of freesurfer data")

# ------------------------------------------------------------------------------ #
# post_structural
# Extract items containing numbers in name
atlas = glob.glob(dir_fs_label + 'lh.*_mics.annot', recursive=True)
atlas = [f.replace(dir_fs_label, '').replace('.annot','').replace('lh.','') for f in atlas]

# Loop to print the filenames
for annot in atlas:
    nom = subBIDS + "_atlas-" + annot + "_desc-surf.png"
    fileL= dir_fS + subBIDS + '/label/lh.' + annot + '.annot'
    fileR= dir_fS + subBIDS + '/label/rh.' + annot + '.annot'
    nomPng = dir_QC_png + nom
    label = np.concatenate((nb.freesurfer.read_annot(fileL)[0], nb.freesurfer.read_annot(fileR)[0]), axis=0)
    if os.path.exists(fileL) and os.path.exists(fileR):
        print("[INFO].... Creating PNG of " + annot + " on native surface")
        plot_hemispheres(surf_lh, surf_rh, array_name=label, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), cmap=cmap_gradient(len(np.unique(label)), ['inferno', 'hsv', 'hsv', 'tab20b']), transparent_bg=False,
                         screenshot = True, filename = nomPng)

# ------------------------------------------------------------------------------ #
# MPC
mask =np.hstack( np.where(th < 0.5, 0, 1) )
for k in range(1,15):
    k=str(k)
    k0=k.rjust(2,'0')
    mpc_lh = dir_mpc + subBIDS + '_space-fsnative_desc-lh_MPC-' + k + '.mgh'
    mpc_rh = dir_mpc + subBIDS + '_space-fsnative_desc-rh_MPC-' + k + '.mgh'
    if os.path.exists(mpc_lh) and os.path.exists(mpc_rh):
        nom = subBIDS + "_space-fsnative_desc-surf_MPC-" + k0 + ".png"
        nomPng = dir_QC_png + nom
        print("[INFO].... Creating PNG of MPC-" + k0 )
        mpc = np.hstack(np.concatenate((np.array(load(mpc_lh).get_fdata()), np.array(load(mpc_rh).get_fdata())), axis=0))*mask
        Qt = (round(np.quantile(mpc[np.nonzero(mpc)],0.05),0), round(np.quantile(mpc[np.nonzero(mpc)],0.95),0))
        plot_hemispheres(surf_lh, surf_rh, array_name=mpc, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
            nan_color=(0, 0, 0, 1), color_range=Qt, cmap="viridis",transparent_bg=False,
            screenshot = True, filename = nomPng)

# ------------------------------------------------------------------------------ #
# Morphology
try:
    # Curvature native inflated surface
    cv_lh = dir_morph + subBIDS + '_space-fsnative_desc-lh_curvature.mgh'
    cv_rh = dir_morph + subBIDS + '_space-fsnative_desc-rh_curvature.mgh'
    if os.path.exists(cv_lh) and os.path.exists(cv_rh):
        CV = np.hstack(np.concatenate((np.array(load(cv_lh).get_fdata()), np.array(load(cv_rh).get_fdata())), axis=0))
        plot_hemispheres(inf_lh, inf_rh, array_name=CV, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv, transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_curv_morph.png')

    # Curvature inflated fsaverage 5
    cv_lh_fs5 = dir_morph + subBIDS + '_space-fsaverage5_desc-lh_curvature.mgh'
    cv_rh_fs5 = dir_morph + subBIDS + '_space-fsaverage5_desc-rh_curvature.mgh'
    if os.path.exists(cv_lh_fs5) and os.path.exists(cv_rh_fs5):
        CV_fs5 = np.hstack(np.concatenate((np.array(load(cv_lh_fs5).get_fdata()), np.array(load(cv_rh_fs5).get_fdata())), axis=0))
        plot_hemispheres(fs5I_lh, fs5I_rh, array_name=CV_fs5, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv, transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsaverage5_desc-surf_curv_morph.png')

    # Curvature fsaverage 5 fwhm=10mm
    cv_lh_fs5S = dir_morph + subBIDS + '_space-fsaverage5_desc-lh_curvature_10mm.mgh'
    cv_rh_fs5S = dir_morph + subBIDS + '_space-fsaverage5_desc-rh_curvature_10mm.mgh'
    if os.path.exists(cv_lh_fs5S) and os.path.exists(cv_rh_fs5S):
        CV_fs5S = np.hstack(np.concatenate((np.array(load(cv_lh_fs5S).get_fdata()), np.array(load(cv_rh_fs5S).get_fdata())), axis=0))
        plot_hemispheres(fs5I_lh, fs5I_rh, array_name=CV_fs5S, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv, transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsaverage5_desc-surf_curv_10mm_morph.png')

    # Curvature conte69
    cv_lh_c69 = dir_morph + subBIDS + '_space-conte69-32k_desc-lh_curvature.mgh'
    cv_rh_c69 = dir_morph + subBIDS + '_space-conte69-32k_desc-rh_curvature.mgh'
    if os.path.exists(cv_lh_c69) and os.path.exists(cv_rh_c69):
        CV_c69 = np.hstack(np.concatenate((np.array(load(cv_lh_c69).get_fdata()), np.array(load(cv_rh_c69).get_fdata())), axis=0))
        plot_hemispheres(c69_lh, c69_rh, array_name=CV_c69, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv, transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_curv_morph.png')

    # Curvature conte69 fwhm=10mm
    cv_lh_c69S = dir_morph + subBIDS + '_space-conte69-32k_desc-lh_curvature_10mm.mgh'
    cv_rh_c69S = dir_morph + subBIDS + '_space-conte69-32k_desc-rh_curvature_10mm.mgh'
    if os.path.exists(cv_lh_c69S) and os.path.exists(cv_rh_c69S):
        CV_c69S = np.hstack(np.concatenate((np.array(load(cv_lh_c69S).get_fdata()), np.array(load(cv_rh_c69S).get_fdata())), axis=0))
        plot_hemispheres(c69_lh, c69_rh, array_name=CV_c69S, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv, transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_curv_10mm_morph.png')

    # Thickness native inflated
    th_lh = dir_morph + subBIDS + '_space-fsnative_desc-lh_thickness.mgh'
    th_rh = dir_morph + subBIDS + '_space-fsnative_desc-rh_thickness.mgh'
    if os.path.exists(th_lh) and os.path.exists(th_rh):
        TH = np.hstack(np.concatenate((np.array(load(th_lh).get_fdata()), np.array(load(th_rh).get_fdata())), axis=0))*mask
        plot_hemispheres(inf_lh, inf_rh, array_name=th, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno", transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_th_morph.png')

    # Thickness fsaverage 5 inflated
    th_lh_fs5 = dir_morph + subBIDS + '_space-fsaverage5_desc-lh_thickness.mgh'
    th_rh_fs5 = dir_morph + subBIDS + '_space-fsaverage5_desc-rh_thickness.mgh'
    if os.path.exists(th_lh_fs5) and os.path.exists(th_rh_fs5):
        th_fs5 = np.hstack(np.concatenate((np.array(load(th_lh_fs5).get_fdata()), np.array(load(th_rh_fs5).get_fdata())), axis=0))
        plot_hemispheres(fs5I_lh, fs5I_rh, array_name=th_fs5, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno", transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsaverage5_desc-surf_th_morph.png')

    # Thickness fsaverage 5 fwhm=10mm
    th_lh_fs10 = dir_morph + subBIDS + '_space-fsaverage5_desc-lh_thickness_10mm.mgh'
    th_rh_fs10 = dir_morph + subBIDS + '_space-fsaverage5_desc-rh_thickness_10mm.mgh'
    if os.path.exists(th_lh_fs10) and os.path.exists(th_rh_fs10):
        th_fs10 = np.hstack(np.concatenate((np.array(load(th_lh_fs10).get_fdata()), np.array(load(th_rh_fs10).get_fdata())), axis=0))
        plot_hemispheres(fs5I_lh, fs5I_rh, array_name=th_fs10, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno", transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsaverage5_desc-surf_th_10mm_morph.png')

    # thickness conte69
    th_lh_c69 = dir_morph + subBIDS + '_space-conte69-32k_desc-lh_thickness.mgh'
    th_rh_c69 = dir_morph + subBIDS + '_space-conte69-32k_desc-rh_thickness.mgh'
    if os.path.exists(th_lh_c69) and os.path.exists(th_rh_c69):
        th_c69 = np.hstack(np.concatenate((np.array(load(th_lh_c69).get_fdata()), np.array(load(th_rh_c69).get_fdata())), axis=0))
        plot_hemispheres(c69_lh, c69_rh, array_name=th_c69, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno", transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_th_morph.png')

    # thickness conte69 fwhm=10mm
    th_lh_c69S = dir_morph + subBIDS + '_space-conte69-32k_desc-lh_thickness_10mm.mgh'
    th_rh_c69S = dir_morph + subBIDS + '_space-conte69-32k_desc-rh_thickness_10mm.mgh'
    if os.path.exists(th_lh_c69S) and os.path.exists(th_rh_c69S):
        th_c69S = np.hstack(np.concatenate((np.array(load(th_lh_c69S).get_fdata()), np.array(load(th_rh_c69S).get_fdata())), axis=0))
        plot_hemispheres(c69_lh, c69_rh, array_name=th_c69S, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno", transparent_bg=False,
                         screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_th_10mm_morph.png')
except ValueError:
    print("[WARNING].... some Morphology surfaces are missing")
else:
    print("[INFO].... Creating PNGs of -Morphology")

# ------------------------------------------------------------------------------ #
# conte69 surfaces
try:
    c69_lhS, c69_rhS = load_surface(dir_conte+subBIDS+'_lh_sphereReg.surf.gii',
                                  dir_conte+subBIDS+'_rh_sphereReg.surf.gii', with_normals=True, join=False)
    c69_lhM, c69_rhM = load_surface(dir_conte+subBIDS+'_space-conte69-32k_desc-lh_midthickness.surf.gii',
                                  dir_conte+subBIDS+'_space-conte69-32k_desc-rh_midthickness.surf.gii', with_normals=True, join=False)
    c69_lhP, c69_rhP = load_surface(dir_conte+subBIDS+'_space-conte69-32k_desc-lh_pial.surf.gii',
                                  dir_conte+subBIDS+'_space-conte69-32k_desc-rh_pial.surf.gii', with_normals=True, join=False)
    c69_lhW, c69_rhW = load_surface(dir_conte+subBIDS+'_space-conte69-32k_desc-lh_white.surf.gii',
                                  dir_conte+subBIDS+'_space-conte69-32k_desc-rh_white.surf.gii', with_normals=True, join=False)
    Val = np.repeat(0, c69_lhM.n_points + c69_rhM.n_points, axis=0)
    grey = plt.colors.ListedColormap(np.full((256, 4), [0.65, 0.65, 0.65, 1]))

    # Sphere native
    plot_hemispheres(c69_lhS, c69_rhS, array_name=th, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="gray", transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-fsnative_desc-surf_sphere.png')

    # conte69 midthickness
    plot_hemispheres(c69_lhM, c69_rhM, array_name=Val, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(-1,1), cmap=grey, transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_mid.png')

    # conte69 pial
    plot_hemispheres(c69_lhP, c69_rhP, array_name=Val, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap=grey, transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_pial.png')

    # conte69 white
    plot_hemispheres(c69_lhW, c69_rhW, array_name=Val, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap=grey, transparent_bg=False,
                     screenshot = True, filename = dir_QC_png + subBIDS + '_space-conte69_desc-surf_white.png')
except ValueError:
    print("[WARNING].... some post_structural conte69 files are missing")
else:
    print("[INFO].... Creating PNGs of conte69 surfaces")
