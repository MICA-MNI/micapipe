#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 19:41:17 2022

@author: Vladimir Fonov, modified by rcruces
   
"""

import nibabel as nib
import argparse
import numpy as np
from time import gmtime, strftime
from copy import deepcopy
import sys


def format_history(argv):
    stamp=strftime("%a %b %d %T %Y>>>", gmtime())
    return stamp+(' '.join(argv))

def mp2rage_robust_combination(filename_uni, filename_inv1, filename_inv2, filename_output = None, multiplying_factor=5):
    """ This is a python version of the 'RobustCombination.mat' function
    described by O'Brien et al. and originally implemented by Jose Marques.
    
    Reference: 
        Marques, J. P., Kober, T., Krueger, G., van der Zwaag, W., Van de Moortele, P. F., & Gruetter, R. (2010). 
        MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. 
        Neuroimage, 49(2), 1271-1281. https://doi.org/10.1016/j.neuroimage.2009.10.002

    Parameters
    ----------
    filename_uni : string
        path to the uniform T1-image (UNI)
    filename_inv1 : string
        path to the first inversion image (INV1)
    filename_inv2 : string
        path to the second inversion image (INV2)
    filename_output : string  (optional)
        path to output image
    multiplying_factor : int  (optional)
        if it's too noisy, give it a bigger value
    """
    # define relevant functions
    mp2rage_robustfunc  =  lambda inv1, inv2, beta: (inv1.conj() * inv2 - beta) / (np.square(inv1) + np.square(inv2) + 2*beta)

    rootsquares_pos  = lambda a,b,c: (-b+np.sqrt(np.square(b) -4 *a*c))/(2*a)
    rootsquares_neg  = lambda a,b,c: (-b-np.sqrt(np.square(b) -4 *a*c))/(2*a)

    # load data
    image_uni  = nib.load(filename_uni)
    image_inv1 = nib.load(filename_inv1)
    image_inv2 = nib.load(filename_inv2)

    image_uni_fdata = image_uni.get_fdata()
    image_inv1_fdata = image_inv1.get_fdata()
    image_inv2_fdata  = image_inv2.get_fdata()
    
    # Maximun intensity of uni image
    uni_max = np.amax(image_uni_fdata)
    
    # scale UNI image values
    if (np.amin(image_uni_fdata) >=0) and (np.amax(image_uni_fdata >= 0.51)):
        scale = lambda x: (x - np.amax(image_uni_fdata)/2) / np.amax(image_uni_fdata)
        scaleback = lambda x: (x*uni_max+uni_max/2)
        image_uni_fdata = scale(image_uni_fdata)

    # correct polarity for INV1
    image_inv1_fdata = np.sign(image_uni_fdata) * image_inv1_fdata

    # MP2RAGEimg is a phase sensitive coil combination.. some more maths has to
    # be performed to get a better INV1 estimate which here is done by assuming
    # both INV2 is closer to a real phase sensitive combination
    inv1_pos = rootsquares_pos(-image_uni_fdata, image_inv2_fdata, -np.square(image_inv2_fdata) * image_uni_fdata)
    inv1_neg = rootsquares_neg(-image_uni_fdata, image_inv2_fdata, -np.square(image_inv2_fdata) * image_uni_fdata)

    image_inv1_final_fdata = deepcopy(image_inv1_fdata)

    image_inv1_final_fdata[np.abs(image_inv1_fdata - inv1_pos) >  np.abs(image_inv1_fdata - inv1_neg)] = inv1_neg[np.abs(image_inv1_fdata - inv1_pos) >  np.abs(image_inv1_fdata - inv1_neg)]
    image_inv1_final_fdata[np.abs(image_inv1_fdata - inv1_pos) <= np.abs(image_inv1_fdata - inv1_neg)] = inv1_pos[np.abs(image_inv1_fdata - inv1_pos) <= np.abs(image_inv1_fdata - inv1_neg)]

    noiselevel = multiplying_factor * np.mean(np.mean(np.mean(image_inv2_fdata[1:,-10:,-10:])))

    image_output = mp2rage_robustfunc(image_inv1_final_fdata, image_inv2_fdata, np.square(noiselevel))

    _history = format_history(sys.argv)
    
    clean_img = nib.Nifti1Image(scaleback(image_output), image_uni.affine, image_uni.header)

    # Save output 
    print(_history)
    if filename_output == None:
        filename_output=filename_uni.replace('.gz','').replace('.nii','') + '_denoised_mf-'+str(multiplying_factor)+'.nii.gz'
    nib.save(clean_img, filename_output)

def parse_options():

    parser = argparse.ArgumentParser(description=""" Remove background noise from mp2rage.\n
                                     This is a python version of the 'RobustCombination.mat' function
    described by O'Brien et al. and originally implemented by Jose Marques.
    Reference: \n
        Marques, J. P., Kober, T., Krueger, G., van der Zwaag, W., Van de Moortele, P. F., & Gruetter, R. (2010). 
        MP2RAGE, a self bias-field corrected sequence for improved segmentation and T1-mapping at high field. 
        Neuroimage, 49(2), 1271-1281. https://doi.org/10.1016/j.neuroimage.2009.10.002,""",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("uni", type=str, default=None,
                    help="UNI NIFTI/NIFTI.GZ file")
    
    parser.add_argument("inv1", type=str,
                        help="INV1 NIFTI/NIFTI.GZ file")

    parser.add_argument("inv2", type=str, nargs='?',
                        help="INV2 NIFTI/NIFTI.GZ file")

    parser.add_argument("output", type=str, nargs='?',
                        help="Output NIFTI.GZ file")

    parser.add_argument("--mf", type=int, default=5,
                        help="Multiplying Factor - Greater value removes more noise")

    params = parser.parse_args()

    params.instance_norm=False

    return params

if __name__ == '__main__':
    params = parse_options()
    mp2rage_robust_combination(params.uni,params.inv1,params.inv2,filename_output=params.output,multiplying_factor=params.mf)
