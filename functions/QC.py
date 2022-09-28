#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MICA pipe Quality Check script:

Generates a pdf file for QC of the processing

    Parameters
    ----------

    subj      : str, subject identification

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

# Arguments
parser = argparse.ArgumentParser()

# Required
parser.add_argument('-subj',
                    dest='subj',
                    type=str,
                    help='Subject identification',
                    nargs='1'
                    required=True
                    )

parser.add_argument('-out',
                    dest='out'
                    type=str,
                    help='Output directory for the processed files <derivatives>'
                    required=True
                    )

parser.add_argument('-bids',
                    dest='bid',
                    type=str,
                    help='Path to BIDS Directory'
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
subj = args.subj
out = os.path.realpath(args.out)
bids = os.path.realpath(args.bids)
ses = args.ses
tracts = args.tracts
tmpDir = args.tmpDir
quiet = args.nocleanup
version = args.version

derivatives = out.split('/micapipe_v1.0.0')[0]


# Dataset name
dataset_description = os.path.realpath("%s/dataset_description.json"%(bids))
with open( dataset_description ) as f:
    dataset_description = json.load(f)
dataset_name = dataset_description["Name"]

# Opitonal inputs:
# Session
if ses == "":
    ses = 'Not defined'

# QC module summary
status =
progress =
time =
threads =
micapipe_version =
module =

# QC png
qc_png = os.path.realpth("")

report_block_header = (
    # ============================================================================================================= #
    # =============================================== P A G E   # 1 =============================================== #
    # ============================================================================================================= #

    # Micapipe banner
    '<img id=\"top\" src=\"${MICAPIPE}/docs/figures/micapipe_long.png\" style=\"width:100%\"  alt=\"micapipe\">'

    # Dataset name
    '<h1 style="color:#343434;font-family:Helvetica, sans-serif !important;text-align:center;margn-bottom:0">'
    '{dataset_name} <h1>'

    # Subject's ID | Session
    '<h3 style="color:#343434;font-family:Helvetica, sans-serif;text-align:center;margin-bottom:0">'
    '<b>Subject</b>: {subj} &nbsp | &nbsp <b>Session</b>: {ses} </h3>'
)

report_module = (
    # ===================================== Structural Prossing - proc_struct ===================================== #
    # Module header:
    '<div class="boxed" style="border:2px solid #666;padding:10px;background-color:#eee;font-family:Helvetica, '
    'sans-serif;font-size:14px">'
    '<b>Module: {module}</b> </div>'

    # QC summary:
    '<h4 style="font-family:Helvetica, sans-serif;text-align:Left;margin-bottom:10px">'
    'Module-specific QC summary </h4>'

    '<style> td {padding:5px;text-align:left} </style>'
    '<table border="1px solid" style="width:100%;border-collapse:collapse">'
        '<tr> <td>Status</td> <td>{status}: {progress} steps completed </td> </tr>'
        '<tr> <td>Processing time</td> <td>{time} minutes</td> </tr>'
        '<tr> <td>Number of threads</td> <td>{threads}</td> </tr>'
        '<tr> <td>Micapipe version</td> <td>{micapipe_version</td> </tr> </table>'

    # Module inputs:
    '<h4 style="font-family:Helvetica, sans-serif;text-align:Left;margin-bottom:10px">'
    'Module-specific input files </h4>'

    '<style> td {padding:5px;text-align:left} </style>'
    '<table border="1px solid" style="width:100%;border-collapse:collapse">'
        '<tr> <td>Status</td> <td>{status}</td> </tr>'
        '<tr> <td>Progress</td> <td>{progress} steps completed</td> </tr>'
        '<tr> <td>Number of threads</td> <td>{threads}</td> </tr>'
        '<tr> <td>Processing time</td> <td>{time} minutes</td> </tr> </table>'

    # Module Outputs
    # T1w nativepro
    '<h4 style="font-family:Helvetica, sans-serif;text-align:Left;margin-bottom:10px">'
    'T1w nativepro </h4>'

    '<img src="{derivatives}/micapipe/{subj}/{ses}"{}'
    '<br>\n'
)


convert_html_to_pdf(, )
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
