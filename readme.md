# micapipe  

## Requirements
| Software   |     Version   | Further info |
|------------|---------------|--------------|  
| dcm2niix   | v1.0.20190902 | https://github.com/rordenlab/dcm2niix |
| Freesurfer | 6.0.0         | https://surfer.nmr.mgh.harvard.edu/ |
| FSl        | 6.0           | https://fsl.fmrib.ox.ac.uk/fsl/fslwiki |
| AFNI       | 20.2.06       | https://afni.nimh.nih.gov/download |
| MRtrix3    | 3.0.1         | https://www.mrtrix.org |
| ANTs       | 2.3.4         | https://github.com/ANTsX/ANTs |
| workbench  | 1.3.2         | https://www.humanconnectome.org/software/connectome-workbench |
| FIX        | 1.06          | https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX |
| R          | 3.6.3         | https://www.r-project.org |
| python     | 3.7.6         | https://www.python.org/downloads/ |

 > The FIX package (FMRIB's ICA-based Xnoiseifier) requires FSL, R and one of MATLAB Runtime Component, full MATLAB or Octave. We recommend the use of the MATLAB Runtime Component. Additionally, it requires the following R libraries:  'kernlab','ROCR','class','party','e1071','randomForest'

### python packages
- ENIGMA toolbox (https://github.com/MICA-MNI/ENIGMA)  

## R libraries  
- 'kernlab'
- 'ROCR'  
- 'class'  
- 'party'  
- 'e1071'  
- 'randomForest'

## Define global Variables
```bash
export MICAPIPE=<github_Directory>/micapipe  
export PATH=$PATH:${MICAPIPE}:${MICAPIPE}/functions  
export CORES=20  
export OMP_NUM_THREADS=4  
export tmp=<path to temporal directory>  
```

## DWI processing (`-proc_dwi`)  
This step requires the `json` files in the BIDS database to contain the following variables:  
```bash
"PhaseEncodingDirection": "j-",
"TotalReadoutTime": 0.047,
```

## Anonymize your dataset
We strongly recommend to anonymize the dataset before running the pipelines. This pipeline includes scripts to do so using AFNI's reface plus tool (https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/refacer/refacer_run.html) to anonymize T1w as well as quantitative T1 images.

## Run `mica-pipe`
proc_rsfmri does not run time slicing correction because our scans where acquired with multiband accelerated fmri.
