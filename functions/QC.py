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

<<<<<<< HEAD
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


# NIFTI_CHECK (generate qc images)
def nifti_check(outName='', outPath='', refPath='', roi=False, figPath=''):

    if os.path.exists(outPath):
        if os.path.exists(refPath):
            ROI = '-roi' if roi else ''
            print(ROI)
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

=======
derivatives = out.split('/micapipe_v0.2.0')[0]
>>>>>>> 74216e79748946c941edeb2bf290c2b5b0427c0e

## ------------------------------------------------------------------------- ##
##                                                                           ##
##                              Build QC report                              ##
##                                                                           ##
## ------------------------------------------------------------------------- ##
static_report = ''


## ---------------------------- MICAPIPE header ---------------------------- ##

# Dataset name
dataset_description = os.path.realpath("%s/dataset_description.json"%(bids))
with open( dataset_description ) as f:
    dataset_description = json.load(f)
dataset_name = dataset_description["Name"]
_static_block = report_header_template(sub=sub, ses_number=ses_number, dataset_name=dataset_name, MICAPIPE=MICAPIPE)

static_report += _static_block


## ------------------------ PROC-STRUCTURAL MODULE ------------------------- ##
_static_block =  report_module_header_template(module='proc_structural')

# QC summary
proc_structural_json = "%s/%s/%s/QC/%s_module-proc_structural.json"%(out,sub,ses,sbids)
_static_block += report_qc_summary_template(proc_structural_json)

# Inputs
nativepro_json = os.path.realpath("%s/%s/%s/anat/%s_space-nativepro_t1w.json"%(out,sub,ses,sbids))
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

T1w_nativepro = "%s/%s/%s/anat/%s_space-nativepro_t1w.nii.gz"%(out,sub,ses,sbids)

figPath = "%s/nativepro_t1w_screenshot.png"%(tmpDir)
_static_block += nifti_check(outName="T1w nativepro", outPath=T1w_nativepro, figPath=figPath)

outPath = "%s/%s/%s/anat/%s_space-nativepro_t1w_brain_mask.nii.gz"%(out,sub,ses,sbids)
figPath = "%s/nativepro_t1w_brain_mask_screenshot.png"%(tmpDir)
_static_block += nifti_check(outName="T1w nativepro brain mask", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

outPath =  "%s/%s/%s/xfm/%s_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz"%(out,sub,ses,sbids)
figPath = "%s/nativepro_t1w_brain_mni152_08_screenshot.png"%(tmpDir)
_static_block += nifti_check(outName="Registration: T1w nativepro in MNI152 0.8mm", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

outPath =  "%s/%s/%s/xfm/%s_from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN_1Warp.nii.gz"%(out,sub,ses,sbids)
figPath = "%s/nativepro_t1w_brain_mni152_08_screenshot.png"%(tmpDir)
_static_block += nifti_check(outName="Registration: T1w nativepro in MNI152 2mm", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

outPath =  "%s/%s/%s/xfm/%s_from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN_1Warp.nii.gz"%(out,sub,ses,sbids)
figPath = "%s/nativepro_t1w_brain_mni152_08_screenshot.png"%(tmpDir)
_static_block += nifti_check(outName="Registration: T1w nativepro in MNI152 2mm", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

outPath =  "%s/%s/%s/anat/%s_space-nativepro_t1w_brain_pve_2.nii.gz"%(out,sub,ses,sbids)
figPath = "%s/nativepro_t1w_brain_pve_2_screenshot.png"%(tmpDir)
_static_block += nifti_check(outName="Partial volume: white matter", outPath=outPath, refPath=T1w_nativepro, figPath=figPath)

# outName = "T1w nativepro white matter partial volume image (pve2)"
# outPath = "%s/%s/%s/xfm/%s_from-nativepro_brain_to-MNI152_2mm_mode-image_desc-SyN_1Warp.nii.gz"%(out,sub,ses,sbids)
# figPath = "%s/nativepro_t1w_brain_mni152_2_screenshot.png"%(tmpDir)
# if os.path.exists(outPath):
#     os.system("${MICAPIPE}/functions/nifti_capture.py -img %s %s -out %s -roi"%("%s/%s/%s/anat/%s_space-nativepro_t1w.nii.gz"%(out,sub,ses,sbids), outPath, figPath))
#     _static_block += report_module_output_template(outName=outName, outPath=outPath, figPath=figPath)
# else:
#     _static_block += ('<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
#         '<b> {outName} </b> </p>'
#         '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin:0px">'
#         'Filepath: does not exist </p>').format(outName=outName)
#
# outName = "T1w nativepro cerebellum atlas"
# outPath = "%s/%s/%s/anat/volumetric/%s_space-nativepro_t1w_atlas-cerebellum.nii.gz"%(out,sub,ses,sbids)
# figPath = "%s/nativepro_t1w_cerebllum_screenshot.png"%(tmpDir)
# if os.path.exists(outPath):
#     os.system("${MICAPIPE}/functions/nifti_capture.py -img %s %s -out %s -roi"%("%s/%s/%s/anat/%s_space-nativepro_t1w.nii.gz"%(out,sub,ses,sbids), outPath, figPath))
#     _static_block += report_module_output_template(outName=outName, outPath=outPath, figPath=figPath)
# else:
#     _static_block += ('<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
#         '<b> {outName} </b> </p>'
#         '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin:0px">'
#         'Filepath: does not exist </p>').format(outName=outName)
#
# outName = "T1w nativepro subcortical atlas"
# outPath = "%s/%s/%s/anat/volumetric/%s_space-nativepro_t1w_atlas-subcortical.nii.gz"%(out,sub,ses,sbids)
# figPath = "%s/nativepro_t1w_subcortical_screenshot.png"%(tmpDir)
# if os.path.exists(outPath):
#     os.system("${MICAPIPE}/functions/nifti_capture.py -img %s %s -out %s -roi"%("%s/%s/%s/anat/%s_space-nativepro_t1w.nii.gz"%(out,sub,ses,sbids), outPath, figPath))
#     _static_block += report_module_output_template(outName=outName, outPath=outPath, figPath=figPath)
# else:
#     _static_block += ('<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
#         '<b> {outName} </b> </p>'
#         '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin:0px">'
#         'Filepath: does not exist </p>').format(outName=outName)
#
# ################# FIX outPath ####################################
# outName = "T1w nativepro 5 tissue segmentation (5TT)"
# outPath = "%s/%s/%s/anat/%s_space-nativepro_t1w_TT.nii.gz"%(out,sub,ses,sbids)
# figPath = "%s/nativepro_t1w_5TT_screenshot.png"%(tmpDir)
# if os.path.exists(outPath):
#     os.system("${MICAPIPE}/functions/nifti_capture.py -img %s %s -out %s"%("%s/%s/%s/anat/%s_space-nativepro_t1w.nii.gz"%(out,sub,ses,sbids), outPath, figPath))
#     _static_block += report_module_output_template(outName=outName, outPath=outPath, figPath=figPath)
# else:
#     _static_block += ('<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
#         '<b> {outName} </b> </p>'
#         '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin:0px">'
#         'Filepath: does not exist </p>').format(outName=outName)

static_report += _static_block


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
