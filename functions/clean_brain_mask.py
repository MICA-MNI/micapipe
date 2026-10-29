#!/usr/bin/env python3
"""Make a FastScrubber mask safe for FastSurfer's surface reconstruction.

    clean_brain_mask.py IN.mgz OUT.mgz [--keep ASEG.mgz]

1. With --keep, add back every voxel the segmentation labels as deep brain
   (everything except cerebral and cerebellar cortex and CSF). FastScrubber
   decides the outer boundary only: a few white matter voxels left out of the
   mask turn mri_normalize's white matter peak into NaN and the surfaces
   collapse (BIDS_CI sub-mri7T: 1280 voxels, cortex 2.15 -> 0.76 mm).
2. Keep the largest connected piece: on a 3T scan FastScrubber left a detached
   75 mL non-brain piece, which makes mri_normalize fail.
3. Fill interior holes, such as the ventricles.
"""
import argparse

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

# FreeSurfer labels left to FastScrubber: cerebral and cerebellar cortex, CSF
OUTER_LABELS = (3, 42, 8, 47, 24)

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("src")
parser.add_argument("dst")
parser.add_argument("--keep", help="segmentation whose deep brain labels must stay in the mask")
args = parser.parse_args()

img = nib.load(args.src)
mask = np.asarray(img.dataobj) > 0
ml = float(np.prod(img.header.get_zooms()[:3])) / 1000.0
notes = []
if args.keep:
    seg = np.asarray(nib.load(args.keep).dataobj)
    if seg.shape != mask.shape:
        raise SystemExit(f"{args.keep} {seg.shape} does not match {args.src} {mask.shape}")
    deep = (seg > 0) & ~np.isin(seg, OUTER_LABELS)
    notes.append(f"added {(deep & ~mask).sum() * ml:.1f} mL of deep brain labels")
    mask |= deep
labels, n = ndi.label(mask)
if n == 0:
    raise SystemExit(f"{args.src}: the mask is empty")
largest = labels == np.argmax(np.bincount(labels.ravel())[1:]) + 1
clean = ndi.binary_fill_holes(largest)
notes.append(f"{n} piece(s), removed {(mask & ~largest).sum() * ml:.1f} mL outside the largest")
notes.append(f"filled {(clean & ~largest).sum() * ml:.1f} mL of interior holes")
print(f"{args.src}: " + "; ".join(notes))
nib.save(nib.MGHImage(clean.astype(np.uint8), img.affine, img.header), args.dst)
