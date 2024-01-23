#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MICA pipe Quality Check script:

Generates a pdf file for QC of the processing

    Parameters
    ----------

    sub      : str, subject identification

    out       : Output directory for the processed files <derivatives>

    bids      : Path to BIDS Directory

    ses       : OPTIONAL flag that indicates the session name (if omitted will manage ase SINGLE session)

    tmpDir    : OPTIONAL specifcation of temporary directory location <path> (Default is /tmp)

    quiet     : OPTIONAL do not print comments

    nocleanup : OPTIONAL do not delete temporary directory at script completion

    version   : OPTIONAL print software version

    USAGE
    --------
    >>> QC.py -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>

  McGill University, MNI, MICA-lab, Created 13 September 2022
  https://github.com/MICA-MNI/micapipe
  http://mica-mni.github.io/

  author: alexander-ngo
"""

from xhtml2pdf import pisa
import os
import argparse
import json
import glob
import nibabel as nb
import numpy as np
import matplotlib as plt
import matplotlib.pyplot as pltpy
import seaborn
from brainspace.datasets import load_mask
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.mesh.mesh_operations import combine_surfaces
from brainspace.utils.parcellation import map_to_labels
from brainspace.vtk_interface import wrap_vtk, serial_connect
from vtk import vtkPolyDataNormals
from pyvirtualdisplay import Display

# Arguments
parser = argparse.ArgumentParser()

# Required
parser.add_argument('-sub',
                    dest='sub',
                    type=str,
                    help='Subject identification',
                    required=True
                    )

parser.add_argument('-out',
                    dest='out',
                    type=str,
                    help='Output directory for the processed files <derivatives>',
                    required=True
                    )

parser.add_argument('-bids',
                    dest='bids',
                    type=str,
                    help='Path to BIDS Directory',
                    required=True
                    )

parser.add_argument('-ses',
                    dest='ses',
                    type=str,
                    default='',
                    help='Flag that indicates the session name'
                    )

parser.add_argument('-tmpDir',
                    dest='tmpDir',
                    type=str,
                    default='/tmp',
                    help='Path to temporary directory location'
                    )

parser.add_argument('-quiet',
                    action='store_false',
                    help='Print comments? {True, False}'
                    )

parser.add_argument('-nocleanup',
                    action='store_false',
                    help='Delete temporary directory at script completion? {True, False}'
                    )

parser.add_argument('-micapipe',
                    dest='MICAPIPE',
                    type=str,
                    help='Path to MICAPIPE Directory',
                    required=True
                    )


parser.add_argument('-version',
                    action='store_false',
                    help='Print software version? {True, False}'
                    )

args = parser.parse_args()

# Arguments
sub = args.sub
out = os.path.realpath(args.out)
bids = os.path.realpath(args.bids)
ses = args.ses
tmpDir = args.tmpDir
quiet = args.nocleanup
version = args.version
MICAPIPE = args.MICAPIPE

# Optional inputs:
# Session
if ses == "" or ses == "SINGLE":
    ses_number = "SINGLE"
    subj_dir = "%s/%s"%(out,sub)
    sbids = sub
else:
    ses_number = ses
    ses = "ses-" + ses_number
    subj_dir = "%s/%s/%s"%(out,sub,ses)
    sbids = sub + "_" + ses

derivatives = out.split('/micapipe_v0.2.0')[0]
derivatives = derivatives+'/micapipe_v0.2.0'

## ------------------------------------------------------------------------- ##
##                                                                           ##
##                      Helper functions to generate PDF                     ##
##                                                                           ##
## ------------------------------------------------------------------------- ##
def file_exists(file_path):
    return os.path.isfile(file_path)

def check_json_exist(jsonPath=''):
    json_exist = os.path.isfile(jsonPath)

    return json_exist

def check_json_complete(jsonPath=''):
    module_description = os.path.realpath(jsonPath)
    with open( module_description ) as f:
        module_description = json.load(f)
    json_complete = module_description["Status"] == "COMPLETED"

    return json_complete


def report_header_template(sub='', ses_number='', dataset_name='', MICAPIPE=''):
    # Header
    report_header = (
        # Micapipe banner
        '<img id=\"top\" src=\"{MICAPIPE}/img/micapipe_long.png\" alt=\"micapipe\">'

        # Dataset name
        '<h1 style="color:#343434;font-family:Helvetica, sans-serif !important;text-align:center;margn-bottom:0">'
        '{dataset_name} <h1>'

        # Subject's ID | Session
        '<h3 style="color:#343434;font-family:Helvetica, sans-serif;text-align:center;margin-bottom:25px">'
        '<b>Subject</b>: {sub} &nbsp | &nbsp <b>Session</b>: {ses_number} </h3>'
    )
    return report_header.format(sub=sub, ses_number=ses_number, dataset_name=dataset_name, MICAPIPE=MICAPIPE)

# Header template
def report_module_header_template(module=''):
    # Module header:
    report_module_header = (
        '<p style="border:2px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica, '
        'sans-serif;font-size:14px">'
        '<b>Module: {module}</b> <p>'
    )
    return report_module_header.format(module=module)

# QC Summary
def report_qc_summary_template(jsonPath=''):
    module_description = os.path.realpath(jsonPath)
    with open( module_description ) as f:
        module_description = json.load(f)
    module = module_description["Module"]
    status = module_description["Status"]
    progress = module_description["Progress"]
    time = module_description["Processing.time"]
    threads = module_description["Threads"]
    micapipe_version = module_description["micapipeVersion"]
    date = module_description["Date"]

    # QC summary table
    report_qc_summary = (
        '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
        '<b>QC summary</b> </p>'

        '<table style="border:1px solid #666;width:100%">'
            # Status
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Status</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{status}: {progress} steps completed</td></tr>'
            # Processing time
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Processing time</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{time} minutes</td></tr>'
            # Thread
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Number of threads</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{threads}</td></tr>'
            # Micapipe version
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Micapipe version</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{micapipe_version}</td></tr>'
            # Date
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Date</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{date}</td></tr>'
        '</table>'
    )
    return report_qc_summary.format(status=status, progress=progress, time=time, threads=threads, micapipe_version=micapipe_version, date=date)

# INPUT TEMPLATE
def report_module_input_template(inputs=[]):
    # Module inputs
    report_module_input = (
        # Module inputs:
        '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0">'
        '<b>Inputs</b> </p>'
    )

    _input = ''
    for i in inputs:
        _input_template = ('<ul><li>{i}</li></ul>')
        _input += _input_template.format(i=i)
    report_module_input += _input

    return report_module_input

# OUTPUT TEMPLATE
def report_module_output_template(outName='', outPath='', figPath=''):
    # Module Outputs
    report_module_output = (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b> {outName} </b> </p>'

        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin:0px">'
        'Filepath: <wbr>{outPath}</wbr> </p>'

        '<center> <img style="width:500px%;margin-top:0px" src="{figPath}"> </center>'
    )
    return report_module_output.format(outName=outName, outPath=outPath, figPath=figPath)


# NIFTI_CHECK (generate qc images from volumes)
def nifti_check(outName='', outPath='', refPath='', roi=False, figPath=''):

    if os.path.exists(outPath):
        if os.path.exists(refPath):
            ROI = '-roi' if roi else ''
            os.system("%s/functions/nifti_capture.py -img %s %s %s -out %s"%(MICAPIPE, refPath, outPath, ROI, figPath))
        else:
            os.system("%s/functions/nifti_capture.py -img %s -out %s"%(MICAPIPE, outPath, figPath))

        _static_block = report_module_output_template(outName=outName, outPath=outPath, figPath=figPath)
    else:
        _static_block = ('<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> {outName} </b> </p>'
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin:0px">'
            'Filepath: does not exist </p>').format(outName=outName)

    return _static_block

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
    cmap = plt.colors.ListedColormap(colors.tolist())
    return cmap

derivatives = out.split('/micapipe_v0.2.0')[0]

ColCurv= plt.colors.ListedColormap(['#A2CD5A', '#A0CA5B', '#9FC85C', '#9EC55D', '#9DC35E', '#9CC05F', '#9BBE61', '#9ABB62', '#99B963', '#98B664', '#96B465', '#95B166', '#94AF68', '#93AC69', '#92AA6A', '#91A76B', '#90A56C', '#8FA26D', '#8EA06F', '#8C9D70', '#8B9B71', '#8A9972', '#899673', '#889475', '#879176', '#868F77', '#858C78', '#848A79', '#82877A', '#81857C', '#80827D', '#7F807E', '#807D7D', '#827A7A', '#857777', '#877575', '#8A7272', '#8C6F6F', '#8F6C6C', '#916969', '#946666', '#966464', '#996161', '#9B5E5E', '#9D5B5B', '#A05858', '#A25656', '#A55353', '#A75050', '#AA4D4D', '#AC4A4A', '#AF4747', '#B14545', '#B44242', '#B63F3F', '#B93C3C', '#BB3939', '#BE3636', '#C03434', '#C33131', '#C52E2E', '#C82B2B', '#CA2828', '#CD2626'])

fs5I_lh = read_surface(MICAPIPE+'/surfaces/fsaverage5/surf/lh.inflated', itype='fs')
fs5I_rh = read_surface(MICAPIPE+'/surfaces/fsaverage5/surf/rh.inflated', itype='fs')
c69_32k_I_lh = read_surface(MICAPIPE+'/surfaces/fsLR-32k.L.inflated.surf.gii', itype='gii')
c69_32k_I_rh = read_surface(MICAPIPE+'/surfaces/fsLR-32k.R.inflated.surf.gii', itype='gii')
c69_5k_I_lh = read_surface(MICAPIPE+'/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
c69_5k_I_rh = read_surface(MICAPIPE+'/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')
mask_32k = load_mask(join=True)

## ------------------------------------------------------------------------- ##
##                                                                           ##
##                              Build QC report                              ##
##                                                                           ##
## ------------------------------------------------------------------------- ##

## ---------------------------- MICAPIPE header ---------------------------- ##
def qc_header():
    dataset_description = os.path.realpath("%s/dataset_description.json"%(bids))
    with open( dataset_description ) as f:
        dataset_description = json.load(f)
    dataset_name = dataset_description["Name"]
    _static_block = report_header_template(sub=sub, ses_number=ses_number, dataset_name=dataset_name, MICAPIPE=MICAPIPE)

    return _static_block


## ------------------------ PROC-STRUCTURAL MODULE ------------------------- ##
def qc_proc_structural(proc_structural_json=''):

    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='proc_structural')

    # QC summary
    _static_block += report_qc_summary_template(proc_structural_json)

    if not check_json_complete(proc_structural_json):
        return _static_block

    # Inputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Inputs</b> </p>'
    )

    nativepro_json = os.path.realpath("%s/anat/%s_space-nativepro_T1w.json"%(subj_dir,sbids))
    with open( nativepro_json ) as f:
        nativepro_json = json.load(f)

    inputs = nativepro_json["inputsRawdata"].split(' ')
    for i, s in enumerate(inputs):
        outName = "Main scan" if len(inputs) == 1 else "Main scan (run %s)"%(i+1)
        figPath = "%s/t1w_mainScan%s.png"%(tmpDir,i+1)
        _static_block += nifti_check(outName=outName, outPath=s, figPath=figPath)

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs </b> </p>'
    )

    T1w_nativepro = "%s/anat/%s_space-nativepro_T1w.nii.gz"%(subj_dir,sbids)

    figPath = "%s/nativepro_T1w_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="T1w nativepro", outPath=T1w_nativepro, figPath=figPath)

    outPath = "%s/anat/%s_space-nativepro_T1w_brain_mask.nii.gz"%(subj_dir,sbids)
    figPath = "%s/nativepro_T1w_brain_mask_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="T1w nativepro brain mask", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

    outPath = "%s/nativepro_T1w_brain_5tt.nii.gz"%(tmpDir)
    figPath = "%s/nativepro_T1w_brain_5tt_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="T1w nativepro 5 tissue segmentation (5TT)", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

    outPath = "%s/%s_space-MNI152_0.8_T1w_brain.nii.gz"%(tmpDir,sbids)
    MNI152_0_8mm = MICAPIPE + "/MNI152Volumes/MNI152_T1_0.8mm_brain_mask.nii.gz"
    figPath = "%s/nativepro_T1w_brain_mni152_08_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: T1w nativepro in MNI152 0.8mm", outPath=MNI152_0_8mm, refPath=outPath, figPath=figPath)

    outPath = "%s/%s_space-MNI152_2_T1w_brain.nii.gz"%(tmpDir,sbids)
    MNI152_2mm = MICAPIPE + "/MNI152Volumes/MNI152_T1_2mm_brain_mask.nii.gz"
    figPath = "%s/nativepro_T1w_brain_mni152_2_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: T1w nativepro in MNI152 2mm", outPath=MNI152_2mm, refPath=outPath, figPath=figPath)

    outPath =  "%s/anat/%s_space-nativepro_T1w_brain_pve_2.nii.gz"%(subj_dir,sbids)
    figPath = "%s/nativepro_T1w_brain_pve_2_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Partial volume: white matter", outPath=outPath, figPath=figPath)

    return _static_block


## --------------------------- PROC-SURF MODULE ---------------------------- ##
def qc_proc_surf(proc_surf_json=''):

    # QC header
    _static_block = qc_header()

    processing = proc_surf_json.split('proc_surf-')[1].split('.json')[0]
    _static_block +=  report_module_header_template(module="proc_surf (%s)"%(processing))

    # QC summary
    _static_block += report_qc_summary_template(proc_surf_json)

    if not check_json_complete(proc_surf_json):
        return _static_block

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b>Native surfaces</b> </p>'
    )

    global recon, surfaceDir, surf_lh, surf_rh, wm_lh, wm_rh, inf_lh, inf_rh

    # Load native surface
    surfaceDir = "%s/%s"%(out.split('/micapipe_v0.2.0')[0],processing)
    surf_lh = read_surface(surfaceDir+'/'+sbids+'/surf/lh.pial', itype='fs')
    surf_rh = read_surface(surfaceDir+'/'+sbids+'/surf/rh.pial', itype='fs')
    wm_lh = read_surface(surfaceDir+'/'+sbids+'/surf/lh.white', itype='fs')
    wm_rh = read_surface(surfaceDir+'/'+sbids+'/surf/rh.white', itype='fs')
    inf_lh = read_surface(surfaceDir+'/'+sbids+'/surf/lh.inflated', itype='fs')
    inf_rh = read_surface(surfaceDir+'/'+sbids+'/surf/rh.inflated', itype='fs')

    # Native thickness
    dsize = (900, 250)
    display = Display(visible=0, size=dsize)
    display.start()
    th = np.concatenate((nb.freesurfer.read_morph_data(surfaceDir + '/' + sbids + '/surf/lh.thickness'), nb.freesurfer.read_morph_data(surfaceDir + '/' + sbids + '/surf/rh.thickness')), axis=0)
    plot_hemispheres(surf_lh, surf_rh, array_name=th, size=dsize, color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno",transparent_bg=False,
                     screenshot = True, offscreen=True, filename = tmpDir + '/' + sbids + '_space-fsnative_desc-surf_thickness.png')
    # Native curvature
    cv = np.concatenate((nb.freesurfer.read_morph_data(surfaceDir + '/' + sbids + '/surf/lh.curv'), nb.freesurfer.read_morph_data(surfaceDir + '/' + sbids + '/surf/rh.curv')), axis=0)
    plot_hemispheres(inf_lh, inf_rh, array_name=cv, size=dsize, color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap='cividis',transparent_bg=False,
                     screenshot = True, offscreen=True, filename = tmpDir + '/' + sbids + '_space-fsnative_desc-surf_curv.png')
    # Native sulcal depth
    sd = np.concatenate((nb.freesurfer.read_morph_data(surfaceDir + '/' + sbids + '/surf/lh.sulc'), nb.freesurfer.read_morph_data(surfaceDir + '/' + sbids + '/surf/rh.sulc')), axis=0)
    plot_hemispheres(wm_lh, wm_rh, array_name=sd, size=dsize, color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(-5, 5), cmap=ColCurv,transparent_bg=False,
                     screenshot = True, offscreen=True, filename = tmpDir + '/' + sbids + '_space-fsnative_desc-surf_sulc.png')
    display.stop()
    native_surface_table = (
        '<table style="border:1px solid #666;width:100%">'
            # Thickness
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><b>Thickness</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{thPath}"></td></tr>'
            # Curvature
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><b>Curvature</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{cvPath}"></td></tr>'
            # Sulcal depth
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><b>Sulcal depth</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{sdPath}"></td></tr>'
        '</table>'
    )

    _static_block += native_surface_table.format(thPath=tmpDir+'/'+sbids+'_space-fsnative_desc-surf_thickness.png',
        cvPath=tmpDir+'/'+sbids+'_space-fsnative_desc-surf_curv.png',
        sdPath=tmpDir+'/'+sbids+'_space-fsnative_desc-surf_sulc.png',
        parcPath=tmpDir+'/'+sbids+'_space-fsnative_desc-surf_a2009s.png',
        parcDKPath=tmpDir+'/'+sbids+'_space-fsnative_desc-surf_aparc.png'
    )


    surf_json = os.path.realpath("%s/surf/%s_proc_surf-%s.json"%(subj_dir,sbids,processing))
    with open( surf_json ) as f:
        surf_description = json.load(f)
    recon = surf_description["SurfRecon"]
    surfaceDir = surf_description["SurfaceDir"]

    return _static_block


## ------------------------- POST-STRUCTURAL MODULE ------------------------ ##
def qc_post_structural(post_structural_json=''):

    post_struct_json = os.path.realpath("%s/anat/%s_post_structural.json"%(subj_dir,sbids))
    with open( post_struct_json ) as f:
        post_struct_description = json.load(f)
    recon = post_struct_description["SurfRecon"]
    surfaceDir = post_struct_description["SurfaceDir"]

    # Load native surface
    surf_lh = read_surface(surfaceDir+'/'+sbids+'/surf/lh.pial', itype='fs')
    surf_rh = read_surface(surfaceDir+'/'+sbids+'/surf/rh.pial', itype='fs')

    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='post_structural (%s)'%(recon))

    # QC summary
    _static_block += report_qc_summary_template(post_structural_json)

    if not check_json_complete(post_structural_json):
        return _static_block

    # Main outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs </b> </p>'
    )

    # Regitration/atlases
    outPath = post_struct_description["NativeSurfSpace"]["fileName"]
    refPath = "%s/anat/%s_space-fsnative_T1w.nii.gz"%(subj_dir,sbids)
    figPath = "%s/nativepro_T1w_fsnative_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: T1w nativepro in %s native space"%(recon), outPath=outPath, refPath=refPath, figPath=figPath)

    outPath = "%s/parc/%s_space-nativepro_T1w_atlas-cerebellum.nii.gz"%(subj_dir,sbids)
    refPath = "%s/anat/%s_space-nativepro_T1w.nii.gz"%(subj_dir,sbids)
    figPath = "%s/nativepro_T1w_cerebellum_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="T1w nativepro cerebellum atlas", outPath=outPath, refPath=refPath, figPath=figPath, roi=True)

    outPath = "%s/parc/%s_space-nativepro_T1w_atlas-subcortical.nii.gz"%(subj_dir,sbids)
    refPath = "%s/anat/%s_space-nativepro_T1w.nii.gz"%(subj_dir,sbids)
    figPath = "%s/nativepro_T1w_subcortical_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="T1w nativepro subcortical atlas", outPath=outPath, refPath=refPath, figPath=figPath, roi=True)

    # Parcellations
    _static_block += (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b> Parcellations </b> </p>'
    )

    parcellation_table = (
        '<table style="border:1px solid #666;width:100%">'
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:center><b>Parcellation</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Surface labels</b></td></tr>'
    )

    label_dir = os.path.realpath("%s/parc/"%(subj_dir))
    annot_dir = os.path.realpath(surfaceDir+'/'+sbids+'/label/')
    files = os.listdir(label_dir)
    filtered_files = [file for file in files if "cerebellum" not in file and "subcortical" not in file and "fsLR-5k" not in file]
    atlas = sorted([file.split("atlas-")[1].split(".nii.gz")[0] for file in filtered_files])

    for annot in atlas:
        fig = sbids + "_atlas-" + annot + "_desc-surf.png"
        fileL= "%s/lh.%s_mics.annot"%(annot_dir,annot)
        fileR= "%s/rh.%s_mics.annot"%(annot_dir,annot)
        figPath = tmpDir + '/' + fig
        label = np.concatenate((nb.freesurfer.read_annot(fileL)[0], nb.freesurfer.read_annot(fileR)[0]), axis=0)
        if os.path.exists(fileL) and os.path.exists(fileR):
            print("[INFO].... Creating PNG of " + annot + " on native surface")
            dsize = (900, 750)
            display = Display(visible=0, size=dsize)
            display.start()
            plot_hemispheres(surf_lh, surf_rh, array_name=label, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                             nan_color=(0, 0, 0, 1), cmap=cmap_gradient(len(np.unique(label)), ['inferno', 'hsv', 'hsv', 'tab20b']), transparent_bg=False,
                             screenshot = True, offscreen=True, filename = figPath)
            display.stop()
        parcellation_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{annot}</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{figPath}"></td></tr>'
        ).format(annot=annot,figPath=figPath)

    parcellation_table += "</table>"

    _static_block += parcellation_table

    # Surfaces
    _static_block += (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b> Surfaces </b> </p>'
    )

    surf_dir = "%s/surf/"%(subj_dir)
    for i, surf in enumerate(['fsnative', 'fsaverage5', 'fsLR-5k', 'fsLR-32k']):
        lhM, rhM = load_surface(surf_dir+sbids+'_hemi-L_space-nativepro_surf-'+surf+'_label-midthickness.surf.gii',
                                surf_dir+sbids+'_hemi-R_space-nativepro_surf-'+surf+'_label-midthickness.surf.gii', with_normals=True, join=False)
        lhP, rhP = load_surface(surf_dir+sbids+'_hemi-L_space-nativepro_surf-'+surf+'_label-pial.surf.gii',
                                surf_dir+sbids+'_hemi-R_space-nativepro_surf-'+surf+'_label-pial.surf.gii', with_normals=True, join=False)
        lhW, rhW = load_surface(surf_dir+sbids+'_hemi-L_space-nativepro_surf-'+surf+'_label-white.surf.gii',
                                surf_dir+sbids+'_hemi-R_space-nativepro_surf-'+surf+'_label-white.surf.gii', with_normals=True, join=False)
        Val = np.repeat(0, lhM.n_points + rhM.n_points, axis=0)
        grey = plt.colors.ListedColormap(np.full((256, 4), [0.65, 0.65, 0.65, 1]))

        # Sphere native
        dsize = (900, 250)
        display = Display(visible=0, size=dsize)
        display.start()
        plot_hemispheres(lhM, rhM, array_name=Val, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-1,1), cmap=grey, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = tmpDir + '/' + sbids + '_space-nativepro_surf-' + surf + '_label-midthickness.png')
        plot_hemispheres(lhP, rhP, array_name=Val, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap=grey, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = tmpDir + '/' + sbids + '_space-nativepro_surf-' + surf + '_label-pial.png')
        plot_hemispheres(lhW, rhW, array_name=Val, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap=grey, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = tmpDir + '/' + sbids + '_space-nativepro_surf-' + surf + '_label-white.png')
        display.stop()
        surface_table = '<br />' if i != 0 else ''
        surface_table += (
            '<table style="border:1px solid #666;width:100%">'
                 '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>{surf}</b></td></tr>'
                 # Pial
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>Pial surface</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{pPath}"></td></tr>'
                 # Middle
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>Middle surface</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{mPath}"></td></tr>'
                 # White
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>White surface</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{wPath}"></td></tr>'
            '</table>'
        )


        _static_block += surface_table.format(surf=surf,
            pPath=tmpDir+'/'+sbids+'_space-nativepro_surf-'+surf+'_label-pial.png',
            mPath=tmpDir+'/'+sbids+'_space-nativepro_surf-'+surf+'_label-midthickness.png',
            wPath=tmpDir+'/'+sbids+'_space-nativepro_surf-'+surf+'_label-white.png'
        )

    # Morphological outputs:
    _static_block += (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b>Morphological features</b></p>'
    )
    display = Display(visible=0, size=(900, 250))
    display.start()
    for feature in ['curv', 'thickness']:

        feature_title = 'Curvature' if feature=='curv' else 'Thickness'
        feature_cmap = 'cividis' if feature=='curv' else 'rocket'
        feature_crange = (-0.2, 0.2) if feature=='curv' else (1.5,4)

        feature_fsn_lh = "%s/maps/%s_hemi-L_surf-fsnative_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_fsn_rh = "%s/maps/%s_hemi-R_surf-fsnative_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_fsn_png = "%s/%s_surf-fsnative_label-%s.png"%(tmpDir,sbids,feature)
        f = np.concatenate((nb.load(feature_fsn_lh).darrays[0].data, nb.load(feature_fsn_rh).darrays[0].data), axis=0)
        plot_hemispheres(inf_lh, inf_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=feature_crange, cmap=feature_cmap, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = feature_fsn_png)

        feature_fs5_lh = "%s/maps/%s_hemi-L_surf-fsaverage5_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_fs5_rh = "%s/maps/%s_hemi-R_surf-fsaverage5_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_fs5_png = "%s/%s_surf-fsaverage5_label-%s.png"%(tmpDir,sbids,feature)
        f = np.concatenate((nb.load(feature_fs5_lh).darrays[0].data, nb.load(feature_fs5_rh).darrays[0].data), axis=0)
        plot_hemispheres(fs5I_lh, fs5I_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=feature_crange, cmap=feature_cmap, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = feature_fs5_png)

        feature_c69_5k_lh = "%s/maps/%s_hemi-L_surf-fsLR-5k_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_c69_5k_rh = "%s/maps/%s_hemi-R_surf-fsLR-5k_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_c69_5k_png = "%s/%s_surf-fsLR-5k_label-%s.png"%(tmpDir,sbids,feature)
        f = np.concatenate((nb.load(feature_c69_5k_lh).darrays[0].data, nb.load(feature_c69_5k_rh).darrays[0].data), axis=0)
        plot_hemispheres(c69_5k_I_lh, c69_5k_I_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=feature_crange, cmap=feature_cmap, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = feature_c69_5k_png)

        feature_c69_32k_lh = "%s/maps/%s_hemi-L_surf-fsLR-32k_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_c69_32k_rh = "%s/maps/%s_hemi-R_surf-fsLR-32k_label-%s.func.gii"%(subj_dir,sbids,feature)
        feature_c69_32k_png = "%s/%s_surf-fsLR-32k_label-%s.png"%(tmpDir,sbids,feature)
        f = np.concatenate((nb.load(feature_c69_32k_lh).darrays[0].data, nb.load(feature_c69_32k_rh).darrays[0].data), axis=0)
        plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=feature_crange, cmap=feature_cmap, transparent_bg=False,
                         screenshot = True, offscreen=True, filename = feature_c69_32k_png)

        morph_table = '<br />' if feature == 'thickness' else ''
        morph_table += (
            '<table style="border:1px solid #666;width:100%">'
                 '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>{feature_title}</b></td></tr>'
                 # Fsnative
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>Fsnative</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{feature_fsn_png}"></td></tr>'
                 # Fsaverage5
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>Fsaverage5</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{feature_fs5_png}"></td></tr>'
                 # fsLR-5k
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>FsLR (5k)</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{feature_c69_5k_png}"></td></tr>'
                 # fsLR-32k
                 '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>FsLR (32k)</b></td>'
                 '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{feature_c69_32k_png}"></td></tr>'
            '</table>'
        )

        _static_block += morph_table.format(feature_title=feature_title,
            feature_fsn_png=feature_fsn_png,
            feature_fs5_png=feature_fs5_png,
            feature_c69_5k_png=feature_c69_5k_png,
            feature_c69_32k_png=feature_c69_32k_png
        )
    display.stop()
    return _static_block

## --------------------------- PROC_FLAIR MODULE --------------------------- ##
def qc_proc_flair(proc_flair_json=''):

    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='proc_flair')

    # QC summary
    _static_block += report_qc_summary_template(proc_flair_json)

    preproc_flair_json = os.path.realpath("%s/maps/%s_space-nativepro_map-flair.json"%(subj_dir,sbids))
    with open( preproc_flair_json ) as f:
        preproc_flair_json = json.load(f)

    if not check_json_complete(proc_flair_json):
        return _static_block

    # Inputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Inputs</b> </p>'
    )

    outPath = preproc_flair_json["inputNIFTI"]["Name"]
    figPath = "%s/flair_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="FLAIR", outPath=outPath, figPath=figPath)

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
    )

    outPath = preproc_flair_json["fileName"]
    refPath = "%s/anat/%s_space-nativepro_T1w.nii.gz"%(subj_dir,sbids)
    figPath = "%s/flair_nativepro_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: FLAIR in T1w nativepro space", outPath=outPath, refPath=refPath, figPath=figPath)

    _static_block += (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b> FLAIR surface mapping</b> </p>'
    )

    flair_table = (
        '<table style="border:1px solid #666;width:100%">'
    )
    for i, surf in enumerate(['fsnative', 'fsaverage5', 'fsLR-5k', 'fsLR-32k']):
        flair_lh = nb.load('%s/maps/%s_hemi-L_surf-%s_label-midthickness_flair.func.gii'%(subj_dir,sbids,surf)).darrays[0].data
        flair_rh = nb.load('%s/maps/%s_hemi-R_surf-%s_label-midthickness_flair.func.gii'%(subj_dir,sbids,surf)).darrays[0].data
        flair = np.concatenate((flair_lh, flair_rh), axis=0)
        flair_fig = tmpDir + '/' + sbids + '_space-nativepro_surf-' + surf + '_label-midthickness_flair.png'
        if surf == 'fsnative':
            surf_lh = inf_lh
            surf_rh = inf_rh
        elif surf == 'fsaverage5':
            surf_lh = fs5I_lh
            surf_rh = fs5I_rh
        elif surf == 'fsLR-5k':
            surf_lh = c69_5k_I_lh
            surf_rh = c69_5k_I_rh
        elif surf == 'fsLR-32k':
            surf_lh = c69_32k_I_lh
            surf_rh = c69_32k_I_rh

        crange=(np.quantile(flair, 0.15), np.quantile(flair, 0.99))
        display = Display(visible=0, size=(900, 250))
        display.start()
        plot_hemispheres(surf_lh, surf_rh, array_name=flair, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=crange, cmap='afmhot', transparent_bg=False,
                         screenshot = True, offscreen=True, filename = flair_fig)
        display.stop()
        flair_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center;width:30%><b>{surf}</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{flair_fig}"></td></tr>'
        ).format(surf=surf,flair_fig=flair_fig)

    flair_table += (
        '</table>'
    )

    _static_block += flair_table

    return _static_block

## ---------------------------- PROC_DWI MODULE ---------------------------- ##
def qc_proc_dwi(proc_dwi_json=''):

    if 'acq-' in proc_dwi_json:
        tag_dwi = proc_dwi_json.split('module-proc_dwi_')[1].split('.json')[0]
        dwi_dir=f'dwi/{tag_dwi}'
        tag_dwi=f'_{tag_dwi}'
    else:
        tag_dwi=''
        dwi_dir='dwi'

    _static_block = qc_header()
    _static_block +=  report_module_header_template(module=f'proc_dwi{tag_dwi}')

    # QC summary
    _static_block += report_qc_summary_template(proc_dwi_json)

    if not check_json_complete(proc_dwi_json):
        return _static_block

    preproc_dwi_json = os.path.realpath(f'{subj_dir}/{dwi_dir}/{sbids}_desc-preproc_dwi.json')
    with open( preproc_dwi_json ) as f:
        preproc_dwi_json = json.load(f)

    # Inputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Inputs</b> </p>'
    )

    dwipeScan = preproc_dwi_json["DWIpe"]["fileName"].split()
    for i, s in enumerate(dwipeScan):

        acq = s.split("acq-")[1].split("_")[0]
        outPath = tmpDir + '/' + s.split('dwi/')[1].split('.nii.gz')[0] + '_mean.nii.gz'
        outName = "DWI pe %s (mean)"%(acq)
        figPath = "%s/dwipe_scan%s.png"%(tmpDir,i+1)
        _static_block += nifti_check(outName=outName, outPath=outPath, figPath=figPath)

    dwirpeScan = preproc_dwi_json["DWIrpe"]["fileName"].split()
    for i, s in enumerate(dwirpeScan):
        try:
            acq = s.split("acq-")[1].split("_")[0]
        except:
            acq = 'b0'
        outPath = tmpDir + '/' + s.split("dwi/")[1].split('.nii.gz')[0] + '_mean.nii.gz'
        outName = "DWI rpe scan %s (mean)"%(acq)
        figPath = "%s/dwirpe_scan%s.png"%(tmpDir,i+1)
        _static_block += nifti_check(outName=outName, outPath=outPath, figPath=figPath)

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
    )

    dwiflspreproc = preproc_dwi_json["dwiflspreproc"]
    dwiflspreproc_table = (
            '<table style="border:1px solid #666;width:100%">'
                '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>DWI processing options</b></td></tr>'
    )
    for k in dwiflspreproc:
        dwiflspreproc_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:left;width:15%><b>{left}</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:left>{right}</td></tr>'
        ).format(left=k, right=dwiflspreproc[k])
    dwiflspreproc_table += "</table>"

    _static_block += dwiflspreproc_table

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_desc-b0.nii.gz"
    figPath = "%s/%s_space-dwi_desc-b0.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="DWI b0 pre-processed", outPath=outPath, figPath=figPath)

    outPath = "%s/%s_space-dwi_desc-preproc_dwi_mean%s.nii.gz"%(tmpDir,sbids,tag_dwi)
    figPath = "%s/%s_space-dwi_desc-preproc_dwi_mean.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="DWI pre-processed", outPath=outPath, figPath=figPath)

    outPath = "%s/%s_space-dwi_model-CSD_map-FOD_desc-wmNorm%s.nii.gz"%(tmpDir,sbids,tag_dwi)
    figPath = "%s/%s_space-dwi_model-CSD_map-FOD_desc-wmNorm.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="White matter fibre orientation (FOD)", outPath=outPath, figPath=figPath)

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_model-DTI_map-ADC.nii.gz"
    figPath = "%s/%s_space-dwi_model-DTI_map-ADC.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="DTI apparent diffusion coefficient (ADC)", outPath=outPath, figPath=figPath)

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_model-DTI_map-FA.nii.gz"
    figPath = "%s/%s_space-dwi_model-DTI_map-FA.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="DTI fractional anisotropy (FA)", outPath=outPath, figPath=figPath)

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_desc-T1w_nativepro_SyN.nii.gz"
    if not os.path.isfile(outPath):
        outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_desc-T1w_nativepro_Affine.nii.gz"
    refPath = "%s/%s_space-dwi_desc-preproc_dwi_mean%s.nii.gz"%(tmpDir,sbids,tag_dwi)
    figPath = "%s/t1w_dwi_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: T1w in DWI space", outPath=outPath, refPath=refPath, figPath=figPath)

    outPath = "%s/%s_space-dwi_desc-5tt%s.nii.gz"%(tmpDir,sbids,tag_dwi)
    figPath = "%s/%s_space-dwi_desc-5tt.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="5 tissue segmentation (5TT) in DWI space", outPath=outPath, figPath=figPath)

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_desc-gmwmi-mask.nii.gz"
    figPath = "%s/%s_space-dwi_gmwmi-mask.png"%(tmpDir,sbids)
    _static_block += nifti_check(outName="Gray-white matter interface (GMWMI) in DWI space", outPath=outPath, figPath=figPath)

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_atlas-cerebellum.nii.gz"
    refPath = "%s/%s_space-dwi_desc-preproc_dwi_mean.nii.gz"%(tmpDir,sbids)
    figPath = "%s/DWI_cerebellum_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Cerebellum atlas in DWI space", outPath=outPath, refPath=refPath, figPath=figPath, roi=True)

    outPath = f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_atlas-subcortical.nii.gz"
    refPath = "%s/%s_space-dwi_desc-preproc_dwi_mean.nii.gz"%(tmpDir,sbids)
    figPath = "%s/DWI_subcortical_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Subcortical atlas in DWI space", outPath=outPath, refPath=refPath, figPath=figPath, roi=True)

    _static_block += (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b> DWI surface mapping of ADC and FA</b> </p>'
    )

    for measure in ['ADC', 'FA']:
        dti_surf_table = '<br>' if measure == 'FA' else ''

        dti_surf_table += (
            '<table style="border:1px solid #666;width:100%">'
                 '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>{measure} (fsLR-32k)</b></td></tr>'
         ).format(measure=measure)

        for surface in ['midthickness', 'white']:
            measure_c69_32k_lh = "%s/maps/%s_hemi-L_surf-fsLR-32k_label-%s_%s%s.func.gii"%(subj_dir,sbids,surface,measure,tag_dwi)
            measure_c69_32k_rh = "%s/maps/%s_hemi-R_surf-fsLR-32k_label-%s_%s%s.func.gii"%(subj_dir,sbids,surface,measure,tag_dwi)
            measure_c69_32k_png = "%s/%s_surf-fsLR-32k_label-%s_%s%s.png"%(tmpDir,sbids,surface,measure,tag_dwi)
            f = np.concatenate((nb.load(measure_c69_32k_lh).darrays[0].data, nb.load(measure_c69_32k_rh).darrays[0].data), axis=0)
            measure_crange=(np.quantile(f[mask_32k], 0.01), np.quantile(f[mask_32k], 0.98))
            # Replace values in f with NaN where mask_32k is False
            f[mask_32k == False] = np.nan
            display = Display(visible=0, size=(900, 250))
            display.start()
            plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                             nan_color=(0, 0, 0, 1), color_range=measure_crange, cmap='mako', transparent_bg=False,
                             screenshot = True, offscreen=True, filename = measure_c69_32k_png)
            display.stop()
            dti_surf_table += (
                     '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{surface}</b></td>'
                     '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{measure_c69_32k_png}"></td></tr>'
            ).format(surface=surface,
                measure_c69_32k_png=measure_c69_32k_png
            )
        _static_block += dti_surf_table

        _static_block += '</table>'

    return _static_block


## ---------------------------- PROC_FUNC MODULE --------------------------- ##
def qc_proc_func(proc_func_json=''):

    tag = proc_func_json.split('%s_module-proc_func-desc-'%(sbids))[1].split('.json')[0]
    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='proc_func (%s)'%(tag))

    # QC summary
    _static_block += report_qc_summary_template(proc_func_json)

    if not check_json_complete(proc_func_json):
        return _static_block

    func_clean_json = glob.glob("%s/func/desc-%s/volumetric/%s_space-func_desc*_preproc.json"%(subj_dir,tag,sbids))[0]
    with open( func_clean_json ) as f:
        func_clean_json = json.load(f)
    acquisition = func_clean_json["Acquisition"]

    # Inputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Inputs</b> </p>'
    )

    mainScan = func_clean_json["Preprocess"]["MainScan"].split()
    for i, s in enumerate(mainScan):
        mainScan_string = os.path.basename(s)
        outPath = tmpDir + '/' + mainScan_string.split('.nii.gz')[0] + '_mean.nii.gz'
        outName = "Main scan (mean)" if len(mainScan) == 1 else "Main scan - echo %s (mean)"%(i+1)
        figPath = "%s/fmri_mainScan%s.png"%(tmpDir,i+1)
        _static_block += nifti_check(outName=outName, outPath=outPath, figPath=figPath)

    mainPhaseScan = func_clean_json["Preprocess"]["MainPhaseScan"]
    if mainPhaseScan == 'DEFAULT':
        mainPhaseScan = os.getenv('default_mainPhase')
        print(mainPhaseScan)
        outPath = tmpDir + '/' + mainPhaseScan.split('nii.gz')[0] + '_mean.nii.gz'
    else:
        outPath = tmpDir + '/' + os.path.basename(mainPhaseScan).split('.nii.gz')[0] + '_mean.nii.gz'

    figPath = "%s/fmri_mainPhaseScan.png"%(tmpDir)
    _static_block += nifti_check(outName="Main phase scan (mean)", outPath=outPath, figPath=figPath)

    reversePhaseScan = func_clean_json["Preprocess"]["ReversePhaseScan"]
    if reversePhaseScan == 'DEFAULT':
        reversePhaseScan = os.getenv('default_reversePhase')
        outPath = tmpDir + '/' + reversePhaseScan.split('.nii.gz')[0] + '_mean.nii.gz'
    else:
        outPath = tmpDir + '/' + os.path.basename(reversePhaseScan).split('.nii.gz')[0] + '_mean.nii.gz'
    figPath = "%s/fmri_reveresePhaseScan.png"%(tmpDir)
    _static_block += nifti_check(outName="Reverse phase scan (mean)", outPath=outPath, figPath=figPath)

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
    )

    clean_json = os.path.realpath("%s/func/desc-%s/volumetric/%s_space-func_desc-%s_preproc.json"%(subj_dir,tag,sbids,acquisition))
    with open( clean_json ) as f:
        clean_json = json.load(f)

    fmripreproc = clean_json["Preprocess"]
    fmripreproc_table = (
        '<table style="border:1px solid #666;width:100%">'
        '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>fMRI processing options</b></td></tr>'
    )

    if acquisition == "me":
        fmripreproc_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:left;width:50%><b>Tedana</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:left;width:50%>YES</td></tr>'
        )

    preproc_info = ["TotalReadoutTime", "EchoTime", "Melodic", "FIX", "GlobalSignalRegression", "CSFWMSignalRegression", "dropTR"]
    for k in preproc_info:
        fmripreproc_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:left;width:50%><b>{left}</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:left;width:50%>{right}</td></tr>'
        ).format(left=k, right=fmripreproc[k])

    fmripreproc_table += "</table>"

    _static_block += fmripreproc_table

    outPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_brain.nii.gz"%(subj_dir,tag,sbids,acquisition)
    figPath = "%s/func_brain_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="fMRI brain", outPath=outPath, figPath=figPath)

    outPath = "%s/anat/%s_space-nativepro_desc-%s_mean.nii.gz"%(subj_dir,sbids,tag)
    refPath = "%s/anat/%s_space-nativepro_T1w.nii.gz"%(subj_dir,sbids)
    figPath = "%s/fmri_nativepro_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: fMRI in T1w nativepro space", outPath=outPath, refPath=refPath, figPath=figPath)

    outPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-T1w.nii.gz"%(subj_dir,tag,sbids)
    refPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_brain.nii.gz"%(subj_dir,tag,sbids,acquisition)
    figPath = "%s/nativepro_T1w_fmri_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Registration: T1w nativepro in fMRI space", outPath=outPath, refPath=refPath, figPath=figPath)

    outPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_cerebellum.nii.gz"%(subj_dir,tag,sbids,acquisition)
    refPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_brain.nii.gz"%(subj_dir,tag,sbids,acquisition)
    figPath = "%s/fMRI_cerebellum_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Cerebellum atlas in fMRI space", outPath=outPath, refPath=refPath, figPath=figPath, roi=True)

    outPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_subcortical.nii.gz"%(subj_dir,tag,sbids,acquisition)
    refPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_brain.nii.gz"%(subj_dir,tag,sbids,acquisition)
    figPath = "%s/fMRI_subcortical_screenshot.png"%(tmpDir)
    _static_block += nifti_check(outName="Subcortical atlas in fMRI space", outPath=outPath, refPath=refPath, figPath=figPath, roi=True)

    _static_block += '<div style="page-break-after: always;"></div>'

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Framewise displace: fMRI</b> </p>'
    )

    outPath = "%s/func/desc-%s/volumetric/%s_space-func_desc-%s_framewiseDisplacement.png"%(subj_dir,tag,sbids,acquisition)
    _static_block += report_module_output_template(outName='', outPath=outPath, figPath=outPath)

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Signal to Noise Ratio (tSNR)</b> </p>'
    )
    tSNR_L = f"{subj_dir}/func/desc-{tag}/surf/{sbids}_surf-fsnative_hemi-L_tSNR.shape.gii"
    tSNR_R = f"{subj_dir}/func/desc-{tag}/surf/{sbids}_surf-fsnative_hemi-R_tSNR.shape.gii"
    tSNR = np.concatenate((nb.load(tSNR_L).darrays[0].data, nb.load(tSNR_R).darrays[0].data), axis=0)
    tSNR = np.squeeze(tSNR)

    snr_fig = tmpDir + "/" + sbids + "_surf-native_tSNR.png"
    display = Display(visible=0, size=(900, 250))
    display.start()
    plot_hemispheres(inf_lh, inf_rh, array_name=tSNR, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), cmap='magma', color_range=(np.quantile(tSNR, 0.05), np.quantile(tSNR, 0.95)), transparent_bg=False, screenshot = True, offscreen=True, filename = snr_fig)
    display.stop()
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Vertex-wise (fsLR-5k) </b> </p>'
            '<center> <img style="width:500px%;margin-top:0px" src="{snr_fig}"> </center>'
    ).format(snr_fig=snr_fig)

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Functional connectomes</b> </p>'
    )

    fc_file = "%s/func/desc-%s/surf/%s_surf-fsLR-5k_desc-FC.shape.gii"%(subj_dir,tag,sbids)
    fc = nb.load(fc_file).darrays[0].data
    np.seterr(divide='ignore')
    fcz = np.arctanh(fc)
    fcz[~np.isfinite(fcz)] = 0
    fcz = np.triu(fcz,1)+fcz.T
    fc_pos = np.copy(fcz)
    fc_pos[0>fc_pos] = 0
    deg = np.sum(fc_pos,axis=1)
    deg_fig = tmpDir + "/" + sbids + "_surf-fsLR-5k_fc_degree.png"
    display = Display(visible=0, size=(900, 250))
    display.start()
    plot_hemispheres(c69_5k_I_lh, c69_5k_I_rh, array_name=deg, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), cmap='rocket', color_range=(np.quantile(deg, 0.001), np.quantile(deg, 0.95)), transparent_bg=False, screenshot = True, offscreen=True, filename = deg_fig)
    display.stop()
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Vertex-wise (fsLR-5k) </b> </p>'
            '<center> <img style="width:500px%;margin-top:0px" src="{deg_fig}"> </center>'
    ).format(deg_fig=deg_fig)

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Parcellated </b> </p>'
    )

    fc_connectome_table = (
        '<table style="border:1px solid #666;width:100%">'
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:center><b>Parcellation</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Connectome</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Degree</b></td></tr>'
    )

    label_dir = os.path.realpath("%s/parc/"%(subj_dir))
    files = os.listdir(label_dir)
    filtered_files = [file for file in files if "cerebellum" not in file and "subcortical" not in file and "fsLR-5k" not in file]
    atlas = sorted([file.split("atlas-")[1].split(".nii.gz")[0] for file in filtered_files])
    for annot in atlas:
        # fc connectomes
        fc_fig = tmpDir + "/" + sbids + "_surf-fsLR-32k_atlas-" + annot + "_fc.png"
        fc_file = "%s/func/desc-%s/surf/%s_atlas-%s_desc-FC.shape.gii"%(subj_dir,tag,sbids,annot)

        # Load shape.gii
        if os.path.isfile(fc_file):
            fc_mtx = nb.load(fc_file).darrays[0].data
            fc = fc_mtx[49:, 49:]
            fcz = np.arctanh(fc)
            fcz[~np.isfinite(fcz)] = 0
            fcz = np.triu(fcz,1)+fcz.T
            pltpy.imshow(fcz, cmap="Reds", aspect='auto')
            pltpy.savefig(fc_fig)

            # Degree
            deg_fig = tmpDir + "/" + sbids + "_atlas-" + annot + "_fc_degree.png"
            fc_pos = np.copy(fcz)
            fc_pos[0>fc_pos] = 0
            deg = np.sum(fc_pos,axis=1)
            annot_file = MICAPIPE + '/parcellations/' + annot + '_conte69.csv'
            if os.path.isfile(annot_file):
                labels_c69 = np.loadtxt(open(annot_file), dtype=int)
                mask_c69 = labels_c69 != 0

                deg_surf = map_to_labels(deg, labels_c69, fill=np.nan, mask=mask_c69)
                display = Display(visible=0, size=(900, 750))
                display.start()
                plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=deg_surf, size=(900, 750), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                                 nan_color=(0, 0, 0, 1), cmap='OrRd', layout_style='grid', transparent_bg=False,
                                 screenshot = True, offscreen=True, filename = deg_fig)
                display.stop()
                fc_connectome_table += (
                    '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{annot}</b></td>'
                    '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{fc_fig}"></td>'
                    '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{deg_fig}"></td></tr>'
                ).format(annot=annot,fc_fig=fc_fig,deg_fig=deg_fig)

    fc_connectome_table += "</table>"
    _static_block += fc_connectome_table

    # Yeon networks
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Yeo networks (schaefer-400) </b> </p>'
    )

    fc_file = "%s/func/desc-%s/surf/%s_atlas-schaefer-400_desc-FC.shape.gii"%(subj_dir,tag,sbids)
    fc_mtx = nb.load(fc_file).darrays[0].data
    fc = fc_mtx[49:, 49:]
    fcz = np.arctanh(fc)
    fcz[~np.isfinite(fcz)] = 0
    fcz = np.triu(fcz,1)+fcz.T
    fc_pos = np.copy(fcz)
    fc_pos[(0>fc_pos)] = 0

    lh_annot = nb.freesurfer.io.read_annot(MICAPIPE + '/parcellations/lh.schaefer-400_mics.annot')
    rh_annot = nb.freesurfer.io.read_annot(MICAPIPE + '/parcellations/rh.schaefer-400_mics.annot')
    annot = lh_annot[2][1:] + rh_annot[2][1:]

    labels_c69 = np.loadtxt(open(MICAPIPE + '/parcellations/schaefer-400_conte69.csv'))
    mask_c69 = labels_c69 != 0

    networks = ['Vis', 'SomMot', 'Default', 'SalVentAttn', 'DorsAttn', 'Cont', 'Limbic' ]

    yeo_table = (
        '<table style="border:1px solid #666;width:100%">'
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:center><b>Network</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Degree</b></td></tr>'
    )

    for i, n in enumerate(networks):
        idx = [n in str(x) for x in annot]
        net = np.sum(fc_pos[idx,:], axis=0)
        net_surf = map_to_labels(net, labels_c69, fill=np.nan, mask=mask_c69)
        net_fig = '%s/%s_atlas-schaefer400_desc-%s.png'%(tmpDir,sbids,n)
        display = Display(visible=0, size=(900, 250))
        display.start()
        plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=net_surf, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), cmap='RdGy_r', transparent_bg=False,
                         screenshot = True, offscreen=True, filename = net_fig)
        display.stop()
        yeo_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{n}</b></td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{net_fig}"></td>'
        ).format(n=n,net_fig=net_fig)

    yeo_table += (
        '</table>'
    )

    _static_block += yeo_table

    return _static_block

## ------------------------------- SC MODULE ------------------------------ ##
def qc_sc(sc_json=''):
    if 'acq-' in sc_json:
        tag_dwi = 'acq-'+sc_json.split('acq-')[1].split('.json')[0]
        dwi_dir=f'dwi/{tag_dwi}'
        tag_dwi=f'_{tag_dwi}'
        streamlines = sc_json.split('module-SC-')[1].split('.json')[0].split(tag_dwi)[0]
    else:
        tag_dwi=''
        dwi_dir='dwi'
        streamlines = sc_json.split(f'{sbids}_module-SC-')[1].split('.json')[0]

    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module=f"Structural connectomes ({streamlines} streamlines{tag_dwi})")

    # QC summary
    _static_block += report_qc_summary_template(sc_json)

    if not check_json_complete(sc_json):
        return _static_block

    tdi_json = os.path.realpath(f"{subj_dir}/{dwi_dir}/{sbids}_space-dwi_desc-iFOD2-{streamlines}_tractography.json")
    with open( tdi_json ) as f:
        tdi_json = json.load(f)

    # Inputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Inputs</b> </p>'
    )

    dti_fod = tdi_json["fileInfo"]["Name"]
    outPath = tmpDir + '/' + dti_fod.split("dwi/")[1].split('.mif')[0] + '.nii.gz'
    figPath = "%s/dti_FOD.png"%(tmpDir)
    _static_block += nifti_check(outName="Fiber orientation distribution", outPath=outPath, figPath=figPath)

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
    )

    tractography = tdi_json["Tractography"]
    tractography_table = (
            '<table style="border:1px solid #666;width:100%">'
                '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2">Tractography information</td></tr>'
    )
    for k in tractography:
        tractography_table += (
            '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:left;width:20%>{left}</td>'
            '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:left>{right}</td></tr>'
        ).format(left=k, right=tractography[k])
    tractography_table += "</table>"

    _static_block += tractography_table

    outPath = "%s/%s_space-dwi_desc-iFOD2-%s_tdi_mean%s.nii.gz"%(tmpDir,sbids,streamlines,tag_dwi)
    figPath = "%s/tdi_%s.png"%(tmpDir,streamlines)
    _static_block += nifti_check(outName="Track density imaging (%s tracks)"%(streamlines), outPath=outPath, figPath=figPath)

    label_dir = os.path.realpath("%s/parc/"%(subj_dir))
    files = os.listdir(label_dir)
    filtered_files = [file for file in files if "cerebellum" not in file and "subcortical" not in file and "fsLR-5k" not in file]
    atlas = sorted([file.split("atlas-")[1].split(".nii.gz")[0] for file in filtered_files])

    connectomes = ['full-connectome', 'full-edgeLengths'] if tractography["weighted_SC"] == "FALSE" else ['full-connectome', 'full-edgeLengths', 'full-weighted_connectome']

    for connectomeType in connectomes:
        c_file = f"{subj_dir}/{dwi_dir}/connectomes/{sbids}_surf-fsLR-5k_desc-iFOD2-{streamlines}-SIFT2_{connectomeType}.shape.gii"
        c = nb.load(c_file).darrays[0].data
        c = np.log(np.triu(c,1)+c.T)
        c[np.isneginf(c)] = 0
        c[c==0] = np.finfo(float).eps
        deg = np.sum(c,axis=1)

        deg_fig = tmpDir + "/" + sbids + "space-dwi_surf-fsLR-5k_desc-iFOD2-" + streamlines + "SIFT2_" + connectomeType + "_degree.png"
        display = Display(visible=0, size=(900, 250))
        display.start()
        plot_hemispheres(c69_5k_I_lh, c69_5k_I_rh, array_name=deg, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), cmap='BuPu', color_range=(np.quantile(deg,0.25), np.quantile(deg,0.99)), transparent_bg=False, screenshot = True, offscreen=True, filename = deg_fig)
        display.stop()
        vertex_wise = (
            '<center> <img style="width:500px%;margin-top:0px" src="{deg_fig}"> </center>'
        ).format(deg_fig=deg_fig)

        connectome_table = (
            '<table style="border:1px solid #666;width:100%">'
                '<tr><td style=padding-top:4px;padding-left:3px;text-align:center><b>Parcellation</b></td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Full connectomes</b></td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Degree</b></td></tr>'
        )

        for annot in atlas:
            if annot == 'aparc-a2009s':
                continue

            annot_lh_fs5= nb.freesurfer.read_annot(MICAPIPE + '/parcellations/lh.'+annot+'_mics.annot')
            Ndim = max(np.unique(annot_lh_fs5[0]))

            c_fig = tmpDir + "/" + sbids + "space-dwi_atlas-" + annot + "_desc-iFOD2-" + streamlines + "SIFT2_" + connectomeType + ".png"
            c_file = f"{subj_dir}/{dwi_dir}/connectomes/{sbids}_space-dwi_atlas-{annot}_desc-iFOD2-{streamlines}-SIFT2_{connectomeType}.shape.gii"
            c = nb.load(c_file).darrays[0].data
            c = np.log(np.triu(c,1)+c.T)
            c[np.isneginf(c)] = 0
            c[c==0] = np.finfo(float).eps
            pltpy.imshow(c, cmap="Purples", aspect='auto')
            pltpy.savefig(c_fig)

            deg_fig = tmpDir + "/" + sbids + "space-dwi_atlas-" + annot + "_desc-iFOD2-" + streamlines + "SIFT2_" + connectomeType + "_degree.png"
            c = c[49:,49:]
            c = np.delete(np.delete(c, Ndim, axis=0), Ndim, axis=1)
            deg = np.sum(c,axis=1)

            annot_file = MICAPIPE + '/parcellations/' + annot + '_conte69.csv'
            if os.path.isfile(annot_file):
                labels_c69 = np.loadtxt(open(annot_file), dtype=int)
                mask_c69 = labels_c69 != 0

                deg_surf = map_to_labels(deg, labels_c69, fill=np.nan, mask=mask_c69)
                display = Display(visible=0, size=(900, 750))
                display.start()
                plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=deg_surf, size=(900, 750), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                                 nan_color=(0, 0, 0, 1), cmap='BuPu', layout_style='grid', transparent_bg=False,
                                 screenshot = True, offscreen=True, filename = deg_fig)
                display.stop()
                connectome_table += (
                    '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center>{annot}</td>'
                    '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{c_fig}"></td>'
                    '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{deg_fig}"></td></tr>'
                ).format(annot=annot,c_fig=c_fig,deg_fig=deg_fig)

        connectome_table += "</table>"

        connectome_str = (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>{title} ({streamlines})</b> </p>'
                '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
                '<b> Vertex-wise (fsLR-5k) </b> </p>'
                '{vertex_wise}'
                '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
                '<b> Parcellated </b> </p>'
                '{connectome_table}'
        )

        if connectomeType == 'full-connectome':
            _static_block += connectome_str.format(title='Structural connectomes',
                    streamlines=streamlines,
                    vertex_wise=vertex_wise,
                    connectome_table=connectome_table
            )
        elif connectomeType == 'full-edgeLengths':
            _static_block += connectome_str.format(title='Edge length connectomes',
                    streamlines=streamlines,
                    vertex_wise=vertex_wise,
                    connectome_table=connectome_table
            )
        elif connectomeType == 'full-weighted_connectome':
            _static_block += connectome_str.format(title='Weighted structural connectomes',
                    streamlines=streamlines,
                    vertex_wise=vertex_wise,
                    connectome_table=connectome_table
            )

    return _static_block


## ------------------------------- MPC MODULE ------------------------------ ##
def qc_mpc(mpc_json=''):

    if 'MPC-SWM' in mpc_json:
        acq_str = 'MPC-SWM'
        mpc_dir = 'mpc-swm'
    else:
        acq_str = 'MPC'
        mpc_dir = 'mpc'
    acquisition = mpc_json.split(f'{sbids}_module-{acq_str}-')[1].split('.json')[0]

    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module=f'{acq_str} | Microstructural profile covariance ({acquisition})')

    # QC summary
    _static_block += report_qc_summary_template(mpc_json)

    if not check_json_complete(mpc_json):
        return _static_block

    # Inputs:
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main inputs:</b> </p>'
    )

    proc_mpc_json = os.path.realpath(f"{subj_dir}/{mpc_dir}/acq-{acquisition}/{sbids}_{acq_str}-{acquisition}.json")
    with open( proc_mpc_json ) as f:
        mpc_description = json.load(f)
    microstructural_img = mpc_description["microstructural_img"]
    microstructural_reg = mpc_description["microstructural_reg"]

    outPath = microstructural_img
    figPath = f"{tmpDir}/{acquisition}_microstructural_img.png"
    _static_block += nifti_check(outName="Microstructural image", outPath=outPath, figPath=figPath)

    outPath = microstructural_reg
    figPath = f"{tmpDir}/{acquisition}_microstructural_reg.png"
    _static_block += nifti_check(outName="Microstructural registration", outPath=outPath, figPath=figPath)

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
    )

    if 'MPC-SWM' in mpc_json:
        ref_space="nativepro"
        refPath = f"{subj_dir}/anat/{sbids}_space-nativepro_{acquisition}.nii.gz"
    else:
        ref_space="fsnative"
        refPath = f"{subj_dir}/anat/{sbids}_space-fsnative_{acquisition}.nii.gz"
    outPath = f"{subj_dir}/anat/{sbids}_space-{ref_space}_T1w.nii.gz"

    figPath = f"{tmpDir}/{acquisition}_fsnative_screenshot.png"
    _static_block += nifti_check(outName=f"Registration: {acquisition} in {ref_space} space", outPath=outPath, refPath=refPath, figPath=figPath)
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>MPC connectomes</b> </p>'
    )

    mpc_file = f"{subj_dir}/{mpc_dir}/acq-{acquisition}/{sbids}_surf-fsLR-5k_desc-MPC.shape.gii"
    mpc = nb.load(mpc_file).darrays[0].data
    mpc = np.triu(mpc,1)+mpc.T
    mpc[~np.isfinite(mpc)] = np.finfo(float).eps
    mpc[mpc==0] = np.finfo(float).eps
    deg = np.sum(mpc,axis=1)
    deg.shape
    deg_fig = tmpDir + "/" + sbids + "surf-fsLR-5k_desc-" + acquisition + "_mpc_degree.png"
    display = Display(visible=0, size=(900, 250))
    display.start()
    plot_hemispheres(c69_5k_I_lh, c69_5k_I_rh, array_name=deg, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), cmap='crest_r', color_range='sym', transparent_bg=False, screenshot = True, offscreen=True, filename = deg_fig)
    display.stop()

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Vertex-wise (fsLR-5k) </b> </p>'
            '<center> <img style="width:500px%;margin-top:0px" src="{deg_fig}"> </center>'
    ).format(deg_fig=deg_fig)

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Parcellated </b> </p>'
    )

    mpc_connectome_table = (
        '<table style="border:1px solid #666;width:100%">'
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:center><b>Parcellation</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Intensity profiles</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Connectomes</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Degree</b></td></tr>'
    )

    label_dir = os.path.realpath(f"{subj_dir}/parc/")
    files = os.listdir(label_dir)
    filtered_files = [file for file in files if "cerebellum" not in file and "subcortical" not in file and "fsLR-5k" not in file]
    atlas = sorted([file.split("atlas-")[1].split(".nii.gz")[0] for file in filtered_files])
    for annot in atlas:
        if annot == 'aparc-a2009s':
            continue

        # Intensity profiles
        ip_fig = tmpDir + "/" + sbids + "_atlas-" + annot + "_desc-" + acquisition + "_intensity_profiles.png"
        ip_file = f"{subj_dir}/{mpc_dir}/acq-{acquisition}/{sbids}_atlas-{annot}_desc-intensity_profiles.shape.gii"
        ip = nb.load(ip_file).darrays[0].data
        pltpy.imshow(ip, cmap="crest", aspect='auto')
        pltpy.savefig(ip_fig)

        # MPC connectomes
        annot_lh_fs5= nb.freesurfer.read_annot(MICAPIPE + '/parcellations/lh.'+annot+'_mics.annot')
        Ndim = max(np.unique(annot_lh_fs5[0]))

        mpc_fig = tmpDir + "/" + sbids + "_atlas-" + annot + "_desc-" + acquisition + "_mpc.png"
        mpc_file = f"{subj_dir}/{mpc_dir}/acq-{acquisition}/{sbids}_atlas-{annot}_desc-MPC.shape.gii"
        mpc = nb.load(mpc_file).darrays[0].data
        mpc = np.triu(mpc,1)+mpc.T
        mpc = np.delete(np.delete(mpc, 0, axis=0), 0, axis=1)
        mpc = np.delete(np.delete(mpc, Ndim, axis=0), Ndim, axis=1)

        mpc[~np.isfinite(mpc)] = np.finfo(float).eps
        mpc[mpc==0] = np.finfo(float).eps

        pltpy.imshow(mpc, cmap="crest", aspect='auto')
        pltpy.savefig(mpc_fig)

        # Degree
        deg_fig = tmpDir + "/" + sbids + "_atlas-" + annot + "_desc-" + acquisition + "_mpc_degree.png"
        deg = np.sum(mpc,axis=1)

        annot_file = MICAPIPE + '/parcellations/' + annot + '_conte69.csv'
        if os.path.isfile(annot_file):
            labels_c69 = np.loadtxt(open(annot_file), dtype=int)
            mask_c69 = labels_c69 != 0

            deg_surf = map_to_labels(deg, labels_c69, fill=np.nan, mask=mask_c69)
            display = Display(visible=0, size=(900, 750))
            display.start()
            plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=deg_surf, size=(900, 750), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                             nan_color=(0, 0, 0, 1), color_range='sym', cmap='mako', layout_style='grid', transparent_bg=False,
                             screenshot = True, offscreen=True, filename = deg_fig)
            display.stop()

            mpc_connectome_table += (
                '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center>{annot}</td>'
                '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{ip_fig}"></td>'
                '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{mpc_fig}"></td>'
                '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{deg_fig}"></td></tr>'
            ).format(annot=annot,ip_fig=ip_fig,mpc_fig=mpc_fig,deg_fig=deg_fig)

    mpc_connectome_table += "</table>"

    _static_block += mpc_connectome_table

    return _static_block


## -------------------------------- GD MODULE ------------------------------ ##
def qc_gd(gd_json=''):

    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='Geodesic distance')

    # QC summary
    _static_block += report_qc_summary_template(gd_json)

    if not check_json_complete(gd_json):
        return _static_block

    # Outputs
    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>Main outputs</b> </p>'
    )

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
            '<b>GD connectomes</b> </p>'
    )

    gd_file = "%s/dist/%s_surf-fsLR-5k_GD.shape.gii"%(subj_dir,sbids)
    gd = nb.load(gd_file).darrays[0].data
    deg = np.sum(gd,axis=1)
    deg_fig = tmpDir + "/" + sbids + "surf-fsLR-5k_GD_degree.png"

    display = Display(visible=0, size=(900, 250))
    display.start()
    plot_hemispheres(c69_5k_I_lh, c69_5k_I_rh, array_name=deg, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), cmap='Blues', transparent_bg=False, screenshot = True, offscreen=True, filename = deg_fig)
    display.stop()

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Vertex-wise (fsLR-5k) </b> </p>'
            '<center> <img style="width:500px%;margin-top:0px" src="{deg_fig}"> </center>'
    ).format(deg_fig=deg_fig)

    _static_block += (
            '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
            '<b> Parcellated </b> </p>'
    )

    gd_connectome_table = (
        '<table style="border:1px solid #666;width:100%">'
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:center><b>Parcellation</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Connectomes</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:center><b>Degree</b></td></tr>'
    )

    label_dir = os.path.realpath("%s/parc/"%(subj_dir))
    files = os.listdir(label_dir)
    filtered_files = [file for file in files if "cerebellum" not in file and "subcortical" not in file and "fsLR-5k" not in file]
    atlas = sorted([file.split("atlas-")[1].split(".nii.gz")[0] for file in filtered_files])
    for annot in atlas:
        if annot == 'aparc-a2009s':
            continue

        # gd connectomes
        annot_lh_fs5= nb.freesurfer.read_annot(MICAPIPE + '/parcellations/lh.'+annot+'_mics.annot')
        Ndim = max(np.unique(annot_lh_fs5[0]))

        gd_fig = tmpDir + "/" + sbids + "_atlas-" + annot + "_gd.png"
        gd_file = "%s/dist/%s_atlas-%s_GD.shape.gii"%(subj_dir,sbids,annot)
        gd = nb.load(gd_file).darrays[0].data
        pltpy.imshow(gd, cmap="Blues", aspect='auto')
        pltpy.savefig(gd_fig)

        # Degree
        gd = np.delete(np.delete(gd, 0, axis=0), 0, axis=1)
        gd = np.delete(np.delete(gd, Ndim, axis=0), Ndim, axis=1)

        deg_fig = tmpDir + "/" + sbids + "_atlas-" + annot + "_gd_degree.png"
        deg = np.sum(gd,axis=1)

        annot_file = MICAPIPE + '/parcellations/' + annot + '_conte69.csv'
        if os.path.isfile(annot_file):
            labels_c69 = np.loadtxt(open(annot_file), dtype=int)
            mask_c69 = labels_c69 != 0

            deg_surf = map_to_labels(deg, labels_c69, fill=np.nan, mask=mask_c69)
            display = Display(visible=0, size=(900, 750))
            display.start()
            plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=deg_surf, size=(900, 750), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                             nan_color=(1, 1, 1, 1), cmap='Blues', layout_style='grid', transparent_bg=False,
                             screenshot=True, filename=deg_fig)
            display.stop()

            gd_connectome_table += (
                '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center>{annot}</td>'
                '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{gd_fig}"></td>'
                '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{deg_fig}"></td></tr>'
            ).format(annot=annot,gd_fig=gd_fig,deg_fig=deg_fig)

    gd_connectome_table += "</table>"

    _static_block += gd_connectome_table

    return _static_block

#------------------------------------------------------------------------------#
# SWM generation and mapping QC
def qc_swm(swm_json=''):
    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='Superficial White Matter')

     # QC summary
    _static_block += report_qc_summary_template(swm_json)

    if not check_json_complete(swm_json):
        print('INCOMPLETE') # return(_static_block)

    # Outputs
    _static_block += (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>Main outputs</b> </p>')

    # SWM surfaces
    _static_block += (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>SWM surfaces</b> </p>'
                '<br />'
                '<table style="border:1px solid #666;width:100%">'
                '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>fsnative</b></td></tr>')

    def surf_table_row(Title, png_path):
        # Add a new row to the table
        surf_row = (
        '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{Title}</b></td>'
        '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{png_path}"></td></tr>').format(
        Title=Title, png_path=png_path)
        return(surf_row)

    # List all the Right and Left SWM surfaces
    surf_dir = f"{subj_dir}/surf"
    swm_L_files = sorted(glob.glob(f"{surf_dir}/{sbids}_hemi-L_surf-fsnative_label-swm*gii"))
    swm_R_files = sorted(glob.glob(f"{surf_dir}/{sbids}_hemi-R_surf-fsnative_label-swm*gii"))

    # Load each surface and plot them and save png
    for i,_ in enumerate(swm_L_files):
        swm_L = swm_L_files[i]
        swm_R = swm_R_files[i]
        # Set the label name
        swm_label = swm_L.replace(".surf.gii","").split('label-')[1]
        # Load the SWM surface
        lhSWM, rhSWM = load_surface(swm_L, swm_R, with_normals=True, join=False)

        # Plot the surfaces SWM
        dsize = (900, 250)
        display = Display(visible=0, size=dsize)
        display.start()
        swm_png=f"{tmpDir}/{sbids}_space-nativepro_surf-fsnative_label-{swm_label}.png"
        plot_hemispheres(lhSWM, rhSWM, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-1,1), transparent_bg=False,
                         screenshot = True, offscreen=True, filename = swm_png)
        _static_block += surf_table_row(swm_label, swm_png)
        display.stop()

    _static_block += ('</table>')

    #------------------------------------------------------------------------------#
    # SWM maps on surfs
    _static_block += (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>SWM maps</b> </p>')

    # List all the Right and Left SWM surfaces
    map_dir = f"{subj_dir}/maps"
    swm_files = sorted(glob.glob(f"{map_dir}/{sbids}_hemi-L_surf-fsLR-32k_label-swm*.func.gii"))

    # Get the unique maps IDs
    maps_str = list(set([file.split('mm_')[1][:-9] for file in swm_files]))

    # Get the unique surfaces
    surf_str = sorted(list(set([file.split('label-')[1].split('_')[0] for file in swm_files])))

    for measure in maps_str:
            swm_surf_table = '<br>'
            swm_surf_table += (
                '<table style="border:1px solid #666;width:100%">'
                     '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>{measure} (fsLR-32k)</b></td></tr>'
             ).format(measure=measure)

            for i, surface in enumerate(surf_str):
                measure_c69_32k_lh = f"{map_dir}/{sbids}_hemi-L_surf-fsLR-32k_label-{surface}_{measure}.func.gii"
                measure_c69_32k_rh = f"{map_dir}/{sbids}_hemi-R_surf-fsLR-32k_label-{surface}_{measure}.func.gii"
                measure_c69_32k_png = f"{tmpDir}/{sbids}_surf-fsLR-32k_label-{surface}_{measure}.png"
                f = np.concatenate((nb.load(measure_c69_32k_lh).darrays[0].data, nb.load(measure_c69_32k_rh).darrays[0].data), axis=0)

                if i == 0: measure_crange=(np.quantile(f[mask_32k], 0.01), np.quantile(f[mask_32k], 0.98))
                display = Display(visible=0, size=(900, 250))
                display.start()
                # Replace values in f with NaN where mask_32k is False
                f[mask_32k == False] = np.nan
                # Plot the values
                plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                                 nan_color=(0, 0, 0, 1), color_range=measure_crange, cmap='mako', transparent_bg=False,
                                 screenshot = True, offscreen=True, filename = measure_c69_32k_png)
                display.stop()
                swm_surf_table += (
                         '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{surface}</b></td>'
                         '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{measure_c69_32k_png}"></td></tr>'
                ).format(surface=surface,
                    measure_c69_32k_png=measure_c69_32k_png
                )
            _static_block += swm_surf_table

            _static_block += '</table>'

    return _static_block

# Utility function
def convert_html_to_pdf(source_html, output_filename):
    # open output file for writing (truncated binary)
    result_file = open(output_filename, "w+b")

    # convert HTML to PDFt
    pisa_status = pisa.CreatePDF(
            source_html,                # the HTML to convert
            dest=result_file)           # file handle to recieve result

    # close output file
    result_file.close()                 # close output file

    # return True on success and False on errors
    return pisa_status.err

# Generate PDF report of Micapipe QC
qc_module_function = {
   'modules':   ['SWM', 'proc_structural', 'proc_surf', 'post_structural', 'proc_dwi', 'proc_func', 'proc_flair', 'SC', 'MPC', 'GD'],
   'functions': [qc_swm, qc_proc_structural, qc_proc_surf, qc_post_structural, qc_proc_dwi, qc_proc_func, qc_proc_flair, qc_sc, qc_mpc, qc_gd]
}

for i, m in enumerate(qc_module_function['modules']):
    module_qc_json = glob.glob("%s/QC/%s_module-%s*.json"%(subj_dir,sbids,m))
    for j in module_qc_json:
        print('--------------------------------------------------')
        print(f"Running QC: {m}")
        if check_json_exist(j):
            try:
                static_report = qc_module_function['functions'][i](j)
                file_pdf=j.replace('.json','_qc-report.pdf')
                convert_html_to_pdf(static_report, file_pdf)
            except:
                print(f"Module QC failed: {m}")
