import sys
import numpy as np
import nibabel as nib
import pandas as pd

subject = sys.argv[1]  # subject name
labelDir = sys.argv[2] # labelDir='/host/fladgate/local_raid/MICA-MTL/HC12/scan_session_01/proc_struct/surfaces/HC12/label/'
annotName = sys.argv[3] # schaefer_100

# Read annotation files
path_lh = labelDir + 'lh.' + annotName + '_mics.annot'
[labels_lh, ctab_lh, names_lh] = nib.freesurfer.io.read_annot(path_lh, orig_ids=True)

path_rh = labelDir + 'rh.' + annotName + '_mics.annot'
[labels_rh, ctab_rh, names_rh] = nib.freesurfer.io.read_annot(path_rh, orig_ids=True)

# Fill vector to output csv
nativeLength = len(labels_lh)+len(labels_rh)
native_parc = [0] * nativeLength

for x in range(len(labels_lh)):
	native_parc[x] = np.where(ctab_lh[:,4] == labels_lh[x])[0][0]

for x in range(len(labels_rh)):
	native_parc[x + len(labels_lh)] = np.where(ctab_rh[:,4] == labels_rh[x])[0][0] + len(ctab_lh)

outputFile = labelDir + subject + '_' + annotName + '_native.csv'
pd.DataFrame(native_parc).to_csv(outputFile, header=None, index=None)
