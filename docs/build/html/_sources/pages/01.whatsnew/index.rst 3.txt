.. _whatsnew:

.. title:: What's new?

What's new?
================================================

v0.1.2
------------------------
**Fixed**

-  Added missing semicolon to SC line 74

-  Issue with transformations management in SC when only AFFINE was applied

- 'str_dwi_affine' variable name corrected dwi in -SC

-  Typo in notification of completition (a missing 's') in proc_struc

-  Variable changed from th to TH in qc_surf.py

-  Issue with the reo file when the original T1w was nii not nii.gz and single run

**Enhaced**

-  Reorder c69 surfaces, plot T1onDWI either Affine or SyN on `micapipe_qc`

-  Added version to group QC table in `utilities.sh`

-  Added full path to `nifti_capture.py` in proc_dwi

-  Created a conte69 dir with `mkdir -p` in freesurfer directory

-  When using antsApplyTransforms replace transformations in tissue_series with variable ($transformsInv)

-  Added option to apply only an affine registration in proc_dwi

-  Added conte69 surfaces to the freesurfer dir and pipeline

-  Increase compatibility if rpe and pe have different size in proc_dwi

-  Update function from app.add_stylesheet to app.add_css_file (conf.py)

-  Erase all MP-PCA and deGibbs files with micapipe_cleanup

-  Manages single session T1w.nii: compression to NIFTI_GZ

-  So not append invidivual QC log

**Documentation**

-  Updated python libraries in README

-  Added gradients tutorial single subject

-  Update native sphere visualization in R with fsbrain

-  Updated surfaces, updated FAQ and references

-  Surface visualization (python)

-  Update doi and reference of SUDMEX dataset

-  Added ipynb and R files, organized surfaces and gradients

-  Added tutorial 'Matrices', made draft of Surface visualization and Gradients


v0.1.1
------------------------

- Documentation update

- Added a missing string in the output names of 02_proc-dwi.sh "*space-dwi_from-dwi_to-nativepro_mode-image_desc-affine_*"

- Added umask to micapipe_cleanup

- Fix a typo in the mica-pipe help (distortion)

- Rename flags and variables in *micapipe_anonymize* from *refacePro* to *warpface*

- Update print version in mica-pipe

- micapipe_qc: added print info for Sankey diagram and Surfaces visualization

- micapipe_qc: added full path to nifti_capture.py

- Fixed an error in *02_proc-rsfmri.sh*, wrong assignation of fmri_pe!

v0.1.0 (Roadrunner)
------------------------

- We are currently on the initial release version of the **micapipe**

- From now on, we'll keep track of the major changes here

- Start keeping changelog 👾🤓👾
