# BIDS to MICA processing pipeline  

## Software versions
dicom2niix  
Freesurfer  6.0   (https://surfer.nmr.mgh.harvard.edu/)  
FSL         6.0   (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)  
AFNI        1.10  (https://afni.nimh.nih.gov/download)  
MRtrix3     3.0.1 (https://www.mrtrix.org)  
ANTs        2.3.3 (https://github.com/ANTsX/ANTs)  
workbench   1.3.2

## Define global Variables
```bash
export MICAPIPE=<github_Directory>/micapipe  
export PATH=$PATH:${MICAPIPE}:${MICAPIPE}/functions  
export CORES=20  
export tmp=<path to temporal directory>  
export OMP_NUM_THREADS=4
```

## Anonymize your dataset
We strongly recommend to anonymize the dataset before running the pipelines. This pipeline includes scripts to do so using AFNI's reface plus tool (https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/refacer/refacer_run.html) to anonymize T1w as well as quantitative T1 images.

## Run `mica-pipe`
