#!/usr/bin/env python3
"""Make a FastScrubber mask safe for FastSurfer's surface reconstruction.

    clean_brain_mask.py IN.mgz OUT.mgz [--keep ASEG.mgz] [--trim SEG.mgz ...]

1. With --keep, add back every voxel the segmentation labels as deep brain
   (everything except cerebral and cerebellar cortex and CSF). FastScrubber
   decides the outer boundary only, it does not cut into the white matter.
2. Keep the largest connected piece: on a 3T scan FastScrubber left a detached
   75 mL non-brain piece, which makes mri_normalize fail.
3. Fill interior holes, such as the ventricles.
4. With --trim, set every label outside the mask to 0 in these segmentations
   (in place), as FreeSurfer's aseg never leaves the brain mask. mri_normalize
   takes the aseg white matter as control points; where norm.mgz is 0 they
   turn its white matter peak into NaN and the surfaces collapse (BIDS_CI:
   FastSurfer 2.6.0 labels a detached piece 110 mm below the brain as white
   matter on sub-mri3T, and 252 white matter voxels outside the mask were
   enough on sub-mri7T).
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
parser.add_argument("--trim", nargs="+", default=[], help="segmentations to clear outside the mask, in place")
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

for path in args.trim:
    seg_img = nib.load(path)
    seg = np.asarray(seg_img.dataobj)
    if seg.shape != clean.shape:
        raise SystemExit(f"{path} {seg.shape} does not match {args.src} {clean.shape}")
    outside = (seg != 0) & ~clean
    print(f"{path}: cleared {outside.sum() * ml:.1f} mL of labels outside the mask "
          f"({np.isin(seg[outside], (2, 41)).sum() * ml:.1f} mL white matter)")
    if outside.any():
        seg = seg.copy()
        seg[outside] = 0
        nib.save(nib.MGHImage(seg, seg_img.affine, seg_img.header), path)
