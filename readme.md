# BIDS to MICA processing pipeline  

## Software versions
Freesurfer  6.0   (https://surfer.nmr.mgh.harvard.edu/)
FSL         6.0   (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
AFNI        1.10  (https://afni.nimh.nih.gov/download)
MRtrix3     3.0.0 (https://www.mrtrix.org)
ANTs        2.3.3 (https://github.com/ANTsX/ANTs)

## Define global Variables
export MICAPIPE=<github_Directory>/micapipe
export CORES=20
export TMP=<path to temporal directory>

## Anonymize your dataset
We strongly recomend to anonymize the dataset before running the pipelines

## Run `mica-pipe`
