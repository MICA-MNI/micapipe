#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generates PNG images from a nifti for Quality Check in micapipe


    Parameters
    ----------

    img     : str, Path to 3D nifti image(s)

    out     : str, Out name

    roi     : bool, optional
              Plots as ROI if True or overlay if False. Default is False.

    title   : str, optional
              Main title of the image. Default is ''.

    Usage
    -----
    nifti_capture.py -out outname.png -img img1.nii.gz img2nii.gz -roi -title 'My NIFTI'

Created on Tue June 8  2021 (the second year of the pademic)

@author: rcruces
"""


import os
import argparse
from nilearn import plotting
import matplotlib as plt

# Arguments
parser = argparse.ArgumentParser()

parser.add_argument('-img',
                    dest='img',
                    type=str,
                    help='Path to input image(s)',
                    nargs='+'
                    )

parser.add_argument('-roi',
                    action='store_true',
                    help='Is the overlay ROI? {True, False}'
                    )

parser.add_argument('-out',
                    dest='out',
                    type=str,
                    help='Output png name and path'
                    )

parser.add_argument('-title',
                    type=str,
                    dest='title',
                    default='',
                    help='Main title of the image'
                    )

args = parser.parse_args()

# Arguments
out=os.path.realpath(args.out)
Nimg=len(args.img)

# Print the arguments
# print('%s = %d' % (' + '.join([str(i) for i in args.img]) ))
print("\nCheck inputs:")
print(f'  -img   : "{args.img}"')
print(f'   Number of img: "{Nimg}"')
print(f'  -out   : "{out}"')
print(f'  -roi   : "{args.roi}"')
print(f'  -title :  "{args.title}"')


# Mandatory inputs
if Nimg < 1 :
    print("\nERROR at least one image is mandatory\n")
elif Nimg > 2 :
    print("\nWARING the scripts only will use the first 2 images\n")

# Create the PNG
if Nimg == 1 :
    print("\n[INFO]... Creating single image png")
    display = plotting.plot_img(args.img[0], colorbar=True, display_mode='ortho', draw_cross=False,
                                title=args.title, cmap=plt.cm.gray, black_bg=True, output_file=out)
elif Nimg >=2 :
    img=os.path.realpath(args.img[0])
    roi=os.path.realpath(args.img[1])
    if args.roi == False :
        print("\n[INFO]... Creating image overlaying NIFTI png")
        display = plotting.plot_img(roi, colorbar=False, display_mode='ortho', draw_cross=False,
                                    title=args.title, cmap=plt.cm.gist_heat, black_bg=True)
        display.add_overlay(img, threshold=0.1, alpha=0.6, cmap=plt.cm.viridis)
        display.savefig(out)
        display.close()
    else:
        print("\n[INFO]... Creating image overlaying ROI png")
        plotting.plot_roi(roi, img, threshold=0.5, alpha=0.5, draw_cross=False, title=args.title, output_file=args.out)
