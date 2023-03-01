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

    tracts    : OPTIONAL number of streamlines, where 'M' stands for millions (default=40M)

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
import sys
import pandas as pd
import os
import argparse
import json
import glob
import nibabel as nb
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

parser.add_argument('-tracts',
                    dest='tracts',
                    type=str,
                    default='40M',
                    help='Number of streamlines'
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
tracts = args.tracts
tmpDir = args.tmpDir
quiet = args.nocleanup
version = args.version

# Optional inputs:
# Session
if ses == "":
    ses_number = "Not defined"
    sbids = sub
else:
    ses_number = ses
    ses = "ses-" + ses_number
    sbids = sub + "_" + ses

derivatives = out.split('/micapipe_v0.2.0')[0]

# Path to MICAPIPE
MICAPIPE=os.popen("echo $MICAPIPE").read()[:-1]


## ------------------------------------------------------------------------- ##
##                                                                           ##
##                      Helper functions to generate PDF                     ##
##                                                                           ##
## ------------------------------------------------------------------------- ##
def check_json_exist_complete(jsonPath=''):
    # Check if json file exists:
    json_exist = os.path.isfile(jsonPath)

    # Check if module is complete
    module_description = os.path.realpath(jsonPath)
    with open( module_description ) as f:
        module_description = json.load(f)
    json_complete = module_description["Status"] == "COMPLETED"

    return json_exist and json_complete


def report_header_template(sub='', ses_number='', dataset_name='', MICAPIPE=''):
    # Header
    report_header = (
        # Micapipe banner
        '<img id=\"top\" src=\"{MICAPIPE}/docs/figures/micapipe_long.png\" alt=\"micapipe\">'

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
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Status</td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{status}: {progress} steps completed</td></tr>'
            # Processing time
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Processing time</td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{time} minutes</td></tr>'
            # Thread
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Number of threads</td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{threads}</td></tr>'
            # Micapipe version
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Micapipe version</td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{micapipe_version}</td></tr>'
            # Date
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Date</td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{date}</td></tr>'
        '</table>'
    )
    return report_qc_summary.format(status=status, progress=progress, time=time, threads=threads, micapipe_version=micapipe_version, date=date)

# INPUT TEMPLATE
def report_module_input_template(inputs=''):
    # Module inputs
    report_module_input = (
        # Module inputs:
        '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0">'
        '<b>Inputs </b> </p>'
    )
    report_module_input.format(inputs=inputs)

    _input = ''
    for i in inputs.split():
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
            os.system("${MICAPIPE}/functions/nifti_capture.py -img %s %s %s -out %s"%(refPath, outPath, roi, figPath))
        else:
            os.system("${MICAPIPE}/functions/nifti_capture.py -img %s -out %s"%(outPath, figPath))

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
    cmap = plt.colors.ListedColormap(colors)
    return cmap

derivatives = out.split('/micapipe_v0.2.0')[0]

## ------------------------------------------------------------------------- ##
##                                                                           ##
##                              Build QC report                              ##
##                                                                           ##
## ------------------------------------------------------------------------- ##
static_report = ''


## ---------------------------- MICAPIPE header ---------------------------- ##
def qc_header()
    dataset_description = os.path.realpath("%s/dataset_description.json"%(bids))
    with open( dataset_description ) as f:
        dataset_description = json.load(f)
    dataset_name = dataset_description["Name"]
    _static_block = report_header_template(sub=sub, ses_number=ses_number, dataset_name=dataset_name, MICAPIPE=MICAPIPE)

    return _static_block


## ------------------------ PROC-STRUCTURAL MODULE ------------------------- ##

def qc_proc_structural(proc_structural_json=''):

    if check_json_exist_complete(proc_structural_json):

        # QC header
        _static_block = qc_header()
        _static_block +=  report_module_header_template(module='proc_structural')

        # QC summary
        _static_block += report_qc_summary_template(proc_structural_json)

        # Inputs
        nativepro_json = os.path.realpath("%s/%s/%s/anat/%s_space-nativepro_T1w.json"%(out,sub,ses,sbids))
        with open( nativepro_json ) as f:
            nativepro_json = json.load(f)
        inputs = nativepro_json["inputsRawdata"]
        _static_block += report_module_input_template(inputs=inputs)

        static_report += _static_block

        # Outputs
        _static_block = (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>Main outputs </b> </p>'
        )

        T1w_nativepro = "%s/%s/%s/anat/%s_space-nativepro_T1w.nii.gz"%(out,sub,ses,sbids)

        figPath = "%s/nativepro_T1w_screenshot.png"%(tmpDir)
        _static_block += nifti_check(outName="T1w nativepro", outPath=T1w_nativepro, figPath=figPath)

        outPath = "%s/%s/%s/anat/%s_space-nativepro_T1w_brain_mask.nii.gz"%(out,sub,ses,sbids)
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

        outPath =  "%s/%s/%s/anat/%s_space-nativepro_T1w_brain_pve_2.nii.gz"%(out,sub,ses,sbids)
        figPath = "%s/nativepro_T1w_brain_pve_2_screenshot.png"%(tmpDir)
        _static_block += nifti_check(outName="Partial volume: white matter", outPath=outPath, figPath=figPath)

        return _static_block


## --------------------------- PROC-SURF MODULE ---------------------------- ##
def qc_proc_surf(proc_surf_json=''):

    if check_json_exist_complete(proc_surf_json):

        # QC header
        _static_block = qc_header()

        surf_json = os.path.realpath("%s/%s/%s/anat/surf/%s_proc_surf.json"%(out,sub,ses,sbids))
        with open( os.path.realpath(surf_json) ) as f:
            surf_description = json.load(f)
        recon = surf_description["SurfRecon"]
        surfaceDir = surf_description["SurfaceDir"]
        _static_block +=  report_module_header_template(module="proc_surf (%s)"%(surf_description["SurfRecon"]))

        # QC summary
        _static_block += report_qc_summary_template(proc_surf_json)

        # Outputs
        _static_block = (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>Main outputs</b> </p>'
                '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
                '<b> Native surfaces </b> </p>'
        )

        # Native thickness
        th = np.concatenate((nb.freesurfer.read_morph_data(surfaceDir + sbids + '/surf/lh.thickness'), nb.freesurfer.read_morph_data(surfaceDir + sbids + '/surf/rh.thickness')), axis=0)
        plot_hemispheres(surf_lh, surf_rh, array_name=th, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(1.5, 4), cmap="inferno",transparent_bg=False,
                         screenshot = True, filename = tmpDir + sbids + '_space-fsnative_desc-surf_thickness.png')
        # Native curvature
        cv = np.concatenate((nb.freesurfer.read_morph_data(surfaceDir + sbids + '/surf/lh.curv'), nb.freesurfer.read_morph_data(surfaceDir + sbids + '/surf/rh.curv')), axis=0)
        plot_hemispheres(wm_lh, wm_rh, array_name=cv, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-0.2, 0.2), cmap=ColCurv,transparent_bg=False,
                         screenshot = True, filename = tmpDir + sbids + '_space-fsnative_desc-surf_curv.png')
        # Native sulcal depth
        sd = np.concatenate((nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/lh.sulc'), nb.freesurfer.read_morph_data(dir_fS + subBIDS + '/surf/rh.sulc')), axis=0)
        plot_hemispheres(wm_lh, wm_rh, array_name=sd, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-5, 5), cmap='cividis',transparent_bg=False,
                         screenshot = True, filename = tmpDir + sbids + '_space-fsnative_desc-surf_sulc.png')
        # Destrieux atlas (aparc.a2009s)
        parc = np.concatenate((nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/lh.aparc.a2009s.annot')[0], nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/rh.aparc.a2009s.annot')[0]), axis=0)
        plot_hemispheres(surf_lh, surf_rh, array_name=parc, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), cmap=cmap_gradient(len(np.unique(parc)), ['inferno', 'hsv', 'hsv', 'tab20b']),transparent_bg=False,
                         screenshot = True, filename = tmpDir + sbids + '_space-fsnative_desc-surf_a2009s.png')
        # Desikan-Killiany Atlas (aparc)
        parcDK = np.concatenate((nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/lh.aparc.annot')[0], nb.freesurfer.read_annot(dir_fS + subBIDS + '/label/rh.aparc.annot')[0]), axis=0)
        plot_hemispheres(surf_lh, surf_rh, array_name=parcDK, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), cmap=cmap_gradient(len(np.unique(parcDK)), ['inferno', 'hsv', 'hsv', 'tab20b']), transparent_bg=False,
                         screenshot = True, filename = tmpDir + sbids + '_space-fsnative_desc-surf_aparc.png')

        native_surface_table = (
            '<table style="border:1px solid #666;width:100%">'
                # Thickness
                '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Thickness</td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><img src="{thPath}"></td></tr>'
                # Curvature
                '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Curvature</td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><img src="{cvPath}"></td></tr>'
                # Sulcal depth
                '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Sulcal depth</td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><img src="{sdPath}"></td></tr>'
                # Destrieux Atlas (aparc.a2009s)
                '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Destrieuz Atlas (aparc.a2009s)</td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><img src="{parcPath}"></td></tr>'
                # Desikan-Killiany Atlas (aparc)
                '<tr><td style=padding-top:4px;padding-left:3px;text-align:left>Desikan-Killiany Atlas (aparc)</td>'
                '<td style=padding-top:4px;padding-left:3px;text-align:center><img src="{parcDKPath}"></td></tr>'
            '</table>'
        )

        _static_block += native_surface_table.format(thPath=tmpDir+sbids+'_space-fsnative_desc-surf_thickness.png',
            cvPath=tmpDir+sbids+'_space-fsnative_desc-surf_curv.png',
            sdPath=tmpDir+sbids+'_space-fsnative_desc-surf_sulc.png',
            parcPath=tmpDir+sbids+'_space-fsnative_desc-surf_a2009s.png'
            parcDKPath=tmpDir+sbids+'_space-fsnative_desc-surf_aparc.png'
        )

        return _static_block

## ------------------------- POST-STRUCTURAL MODULE ------------------------ ##

# QC summary =

# Utility function
def convert_html_to_pdf(source_html, output_filename):
    # open output file for writing (truncated binary)
    result_file = open(output_filename, "w+b")

    # convert HTML to PDF
    pisa_status = pisa.CreatePDF(
            source_html,                # the HTML to convert
            dest=result_file)           # file handle to recieve result

    # close output file
    result_file.close()                 # close output file

    # return True on success and False on errors
    return pisa_status.err

# Generate PDF report of Micapipe QC
convert_html_to_pdf(static_report, out+"/"+sub+'/'+ses+'/QC/'+sbids+'_test.pdf')