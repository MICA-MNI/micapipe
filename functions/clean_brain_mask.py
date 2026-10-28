#!/usr/bin/env python3
"""Keep the largest connected piece of a binary brain mask and fill its interior holes.

FastScrubber can leave detached non-brain pieces (75 mL on a 3T scan), which make
FreeSurfer's mri_normalize fail, and holes in the ventricles.

    clean_brain_mask.py IN.mgz OUT.mgz
"""
import sys

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

src, dst = sys.argv[1:3]
img = nib.load(src)
mask = np.asarray(img.dataobj) > 0
labels, n = ndi.label(mask)
if n == 0:
    sys.exit(f"{src}: the mask is empty")
largest = labels == np.argmax(np.bincount(labels.ravel())[1:]) + 1
clean = ndi.binary_fill_holes(largest)
ml = float(np.prod(img.header.get_zooms()[:3])) / 1000.0
print(f"{src}: {n} piece(s); removed {(mask & ~largest).sum() * ml:.1f} mL outside the largest, "
      f"filled {(clean & ~largest).sum() * ml:.1f} mL of interior holes")
nib.save(nib.MGHImage(clean.astype(np.uint8), img.affine, img.header), dst)
