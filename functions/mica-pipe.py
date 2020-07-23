#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=================================
BIDS: MICA structural processing
=================================

Preprocessing workflow for structural T1w.

This workflow makes use of ANTS, FSL and Freesurfer

Atlas an templates are avaliable from:

https://github.com/MICA-MNI/micaopen/templates


Created on Mon Jun  1 17:40:23 2020

@author: rcruces

"""

from os.path import join as opj
import os
import json
from nipype.interfaces import fsl

from nipype.interfaces.spm import Smooth
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink
from nipype.algorithms.rapidart import ArtifactDetect
from nipype import Workflow, Node

import copy
from nipype.interfaces.ants import N4BiasFieldCorrection
import nipype.interfaces.mrtrix3 as mrt

gen5tt = mrt.Generate5tt()
n4 = N4BiasFieldCorrection()


# Variables USE parse argument
bids_dir = '~/tmp/BIDS'  #experiment_dir
output_dir = '~/tmp/derivatives'
working_dir = '/tmp/'
participant_label='sub-HC10'

# ----------------------------------------------------------------
#           T1w processing
# check https://github.com/llevitis/APPIAN

# Initiate a node to Reorient to RPI with AFNI
reorient = Node(afni.)

# Initiate a node to merge and average
    # note: ONLY IF multiple T1w were provided
flirt
fslmaths

# Initiate a node: Intensity Non-uniform correction
N4BiasFieldCorrection

# Initiate a node: ImageMath
ImageMath


# Initiate a node: brain extraction bet
btr = fsl.BET()

# Initiate a node: run first_all
first = fsl.FIRST()

# Initiate a node: Register to MNI $mmTemplates
aw = fsl.ApplyWarp()
applyxfm = fsl.preprocess.ApplyXFM()

# Initiate a node: 5ttgen

# Initiate a node: freesurfer
from nipype.interfaces.freesurfer import ReconAll
reconall = ReconAll()


# ----------------------------------------------------------------
#           DWI processing
