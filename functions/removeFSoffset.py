#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Handling surfaces
-----------------

"""
import os
import sys
import numpy as np
import nibabel as nb


def normalize_surfs(in_file, fname):

    img = nb.load(in_file)
    pointset = img.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0]
    coords = pointset.data
    c_ras_keys = ('VolGeomC_R', 'VolGeomC_A', 'VolGeomC_S')
    ras = np.array([[float(pointset.metadata[key])]
                    for key in c_ras_keys])
    # Apply C_RAS translation to coordinates
    pointset.data = coords + ras.T

    # Remove C_RAS translation from metadata to avoid double-dipping in FreeSurfer
    for c in c_ras_keys:
        pointset.meta[f'{c}'] = '0.000000'

    # save
    img.to_filename(fname)
    return os.path.abspath(fname)

normalize_surfs(sys.argv[1], sys.argv[2])
