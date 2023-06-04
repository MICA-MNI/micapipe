#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generates PNG images from a nifti for Quality Check in micapipe


    Parameters
    ----------

    conn     : str, path to connectivity matrix

    lut1     : str, path to LUT

    lut2     : str, path to LUT

    mica   : str, MICAPIPE path

    Usage
    -----
    connectome_slicer.py --conn outname.png --lut1 --lut2 --mica

Created on June 2023 (the year of light)

@author: rcruces
"""

import argparse
import pandas as pd
import nibabel as nib
import numpy as np
import os

# Function save as gifti
def save_gii(data_array, file_name):
    # Initialize gifti: NIFTI_INTENT_SHAPE - 2005, FLOAT32 - 16
    gifti_data = nib.gifti.GiftiDataArray(data=data_array, intent=2005, datatype=16)

    # this is the GiftiImage class
    gifti_img = nib.gifti.GiftiImage(meta=None, darrays=[gifti_data])

    # Save the new GIFTI file
    nib.save(img=gifti_img, filename=file_name)

# Argument parsing
parser = argparse.ArgumentParser(description="micapipe connectome slicer")
parser.add_argument("--conn", help="path to connectivity matrix", required=True)
parser.add_argument("--lut1", help="path to LUT", required=True)
parser.add_argument("--lut2", help="path to LUT", required=True)
parser.add_argument("--mica", help="MICAPIPE path", required=True)
args = parser.parse_args()

# Read LUT files
lut1 = pd.read_csv(args.lut1)
lut2 = pd.read_csv(args.lut2)
indx = sorted(list(set(lut1['mics']).union(set(lut2['mics']))))

# Read connectivity matrix
M = pd.read_csv(args.conn, sep=" ", header=None).values
print("INFO.... Connectome dimensions: {} x {}".format(M.shape[0], M.shape[1]))

# Check if connectivity matrix contains "fsLR-5k"
if "fsLR-5k" in args.conn:
    print("Connectome is 'fsLR-5k'. Saving as GIFTI without slicing.")
    output_file = args.conn.replace("txt", "shape.gii")
    save_gii(M, output_file)
    print("GIFTI file saved as:", output_file)
else:
    if M.shape[0] == len(indx):
        print("Connectome is already sliced!")
    else:
        print("INFO.... Connectome new dimensions: {} x {}".format(len(indx), len(indx)))
        indx = [i - 1 for i in indx]
        M = M[np.ix_(indx, indx)]
        # Save the GIFTI data to a file
        output_file = args.conn.replace("txt", "shape.gii")
        save_gii(M, output_file)
        print("GIFTI file saved as:", output_file)

# Remove the txt file
os.remove(args.conn)
