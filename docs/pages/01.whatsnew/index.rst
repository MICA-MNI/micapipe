.. _whatsnew:

.. title:: What's new?

What's new?
================================================

v0.2.+ 'Northern Flicker'
-------------------------------------

- All versions 0.2.0 and above are compatible with each other and have the same software versions.
- Each version update has code improvement and enhancements based on the comments of the users.
- We highly recommend to use the latest version.

v0.2.3, Jan 18, 2024
-------------------------------------

**General**

- Updated RtD
- Added tSNR volumetric and mapping, matched names on all scripts, and updated QC_subj accordingly
- Center-aligned images before affine registration.
- Rename transformation json before it was being overwritten by multiple func acquisitions (e.g. task)
- Quiet percentage of completition of MRtrix
- Center align b0 and T1 for the first affine registration
- New flag -upsample_dwi for proc_dwi
- Fixed MPC cleanup
- Upgraded flair preprocessing with WM-peak, GM mode normalization
- Fixed SINGLE session issues (e.g proc_func and QC)
- Fixed R libraries issues for fls-fix in new Dockerfile

**New Modules**

- Both modules include: json files, subject QC and cleanup integration.
- SWM: superficial white matter generation and mapping in nativepro space.
- SWM_MPC: Superficial white matter MPC in native qMRI space.


v0.2.2, Aug 9, 2023
-------------------------------------

**General**

- FIX issue with multiple acquisitions when `-proc_dwi` was used with `-dwi_acq`. Including the following modules `SC`, `QC` and `QC_subj`.
- FIX instance when multi echo processing has no spike regression file now calculates the nuisance regression correctly.
- FIX an issue with the dimensions of the tSNR file in `proc_func`.

**QC group**

- Set higher thresholds for GD and SC in the similarity subjects matrices for better outlier detection.
- NEW vertex-wise and ROI variance per feature.
- NEW histograms of features and variances per ROI and vertex-wise.

**QC_subj**

- NEW a surface plot of the tSNR to the `proc_func` module in `QC_subj`.


v0.2.1, Jun 26, 2023
-------------------------------------
- We fixed the usage of `-QC` and `-QC_subj` for those datasets with single session.
- Update a conditional in `proc_surf` that search for the `ribbon.mgz`
- Add a missing variable in `-SC` when using `-autotrack`


v0.2.0 'Northern Flicker', Jun 16, 2023
-------------------------------------

- 👥 Compatible with multiple acquisitions (func, dwi, MPC). multiple quantitative maps, DWI acquisitions and func acquisitions including task
- 🧠 Improved registrations between modalities. We improved the registration between multiple acquisitions to optimize the performance.
- 🌐 Surface mapping from native space in all modules.
- 🔬 Optimized for high resolution processing (tested on 0.5-1mm3, e.g 7T)
- 💦 MP2RAGE processing workflow (3T and 7T) for structural images, including a cleaning of the salt and pepper background noise.
- ✨ Optimize to process ultra high resolution 7T datasets, with new CNN based processing tools.
- ⚡️ NEW FLAIR processing module
- 🔍 New QC subject and group level more informative with a lot of new features and pdf (new json files for QC)
- 👾 New algorithms  (brain mask, surface reconstructions - fastsurfer)
- 🛥️ Docker container
- 🐛 A lot of bugs fixed


v0.1.5, Jun 10, 2023
-------------------------------------

**General changes**

- Update print help
- Removed unused flags
- Update the documentation (partially)
- UPDATE jsons to increase `pyBIDS` compatibility
- FIX bug on MPC output directory

**`proc_func` and `FC`**

- Rename flags on `proc_func`  `-regress_WM_CSF` to `-NSR`
- FIX bug on glm regression matrix of cofounding (higher precision)

**`proc_dwi` and `SC` upgrade**

- Manage multiple acquisitions `-dwi_acq`
- `micapipe_cleanup` is up to date to erase multiple acquisitions
- `-b0thr` allows the user to set a threshold to determine the volumes that correspond to b=0, default=61
- `-no_bvalue_scaling` disable the diffusion b-values scaling of the corresponding DWI norm (see Mrtrix3 for further info)


v0.1.4 'Roadrunner', Nov 3, 2022
-------------------------------------

- Update FSL from 6.0.0 to 6.0.3
- `mica-pipe` command goes to `micapipe`
- To print help should specify it with the flag `micapipe -h` or `micapipe -help`
- `post_structural` will always run schefer-400 by default if is not included in the `-atlas` list
- `proc_rsfmri` is deprecated and replaced by `proc_func`
- `proc_func` handles more than one functional acquisition (e.g. tasks), as well as multi echo data (tedana https://github.com/ME-ICA/tedana).
- `proc_func` Added option to drop the first 5 TRs -dropTR (by default is not dropped)
- `proc_func` Added option to not run the functional connectomes -noFC (only func surface data)
- `proc_func` Added 6 parameters of motion to the regression of -regress_WM_CSF (func\~spikes+6motion+wm+csf)
- `proc_func` Added 6 parameters of motion to the regression of -GSR (func\~spikes+6motion+wm+csf+gs)
- `MPC` can processes more than one quantitative map (at a time) with the flag `-mpc_acq <qMRI_name>`
- `micapipe_cleanup` can be called from `micapipe` command: `micapipe -cleanup`
- `micapipe -cleanup` uses the string `-acqStr` to erase multiple acquisitions of `-proc_func` and `-MPC`
- Improved comments and print logs of `-proc_func` and `FC.py`
- `proc_func` generates new jsons files of each acquisition with metadata about processing and completion status
- `proc_func` exit status when Melodic/FIX fail
- `-QC_subj` is not compatible with `proc_func` yet.... or with MPC multiple acquisitions but it is still with the old `proc_rsfmri` outputs
- NOTE: the read the docs is not updated yet.


v0.1.2
-------------------------------------

**Fixed**

- Added missing semicolon to SC line 74
- Issue with transformations management in SC when only AFFINE was applied
- 'str_dwi_affine' variable name corrected dwi in -SC
- Typo in notification of completition (a missing 's') in proc_struc
- Variable changed from th to TH in qc_surf.py
- Issue with the reo file when the original T1w was nii not nii.gz and single run

**Enhaced**

- Reorder c69 surfaces, plot T1onDWI either Affine or SyN on `micapipe_qc`
- Added version to group QC table in `utilities.sh`
- Added full path to `nifti_capture.py` in proc_dwi
- Created a conte69 dir with `mkdir -p` in freesurfer directory
- When using antsApplyTransforms replace transformations in tissue_series with variable ($transformsInv)
- Added option to apply only an affine registration in proc_dwi
- Added conte69 surfaces to the freesurfer dir and pipeline
- Increase compatibility if rpe and pe have different size in proc_dwi
- Update function from app.add_stylesheet to app.add_css_file (conf.py)
- Erase all MP-PCA and deGibbs files with micapipe_cleanup
- Manages single session T1w.nii: compression to NIFTI_GZ
- So not append invidivual QC log

**Documentation**

- Updated python libraries in README
- Added gradients tutorial single subject
- Update native sphere visualization in R with fsbrain
- Updated surfaces, updated FAQ and references
- Surface visualization (python)
- Update doi and reference of SUDMEX dataset
- Added ipynb and R files, organized surfaces and gradients
- Added tutorial 'Matrices', made draft of Surface visualization and Gradients

v0.1.1
-------------------------------------

- Documentation update
- Added a missing string in the output names of 02_proc-dwi.sh "*space-dwi_from-dwi_to-nativepro_mode-image_desc-affine_*"
- Added umask to micapipe_cleanup
- Fix a typo in the mica-pipe help (distortion)
- Rename flags and variables in *micapipe_anonymize* from *refacePro* to *warpface*
- Update print version in mica-pipe
- micapipe_qc: added print info for Sankey diagram and Surfaces visualization
- micapipe_qc: added full path to nifti_capture.py
- Fixed an error in *02_proc-rsfmri.sh*, wrong assignation of fmri_pe!


v0.1.0 'Wobbly'
-------------------------------------

- We are currently on the initial release version of the **micapipe**
- From now on, we'll keep track of the major changes here
- Start keeping changelog 👾🤓👾
