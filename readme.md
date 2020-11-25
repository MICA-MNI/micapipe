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
|     Package     |  Version  |
|:---------------:|:---------:|
| brainspace      | 0.1.1     |
| certifi         | 2020.6.20 |
| cycler          | 0.10.0    |
| enigmatoolbox   | 0.0.1     |
| joblib          | 0.16.0    |
| kiwisolver      | 1.2.0     |
| matplotlib      | 3.3.1     |
| nibabel         | 3.1.1     |
| numpy           | 1.19.1    |
| packaging       | 20.4      |
| pandas          | 1.1.1     |
| Pillow          | 7.2.0     |
| pyparsing       | 2.4.7     |
| python-dateutil | 2.8.1     |
| pytz            | 2020.1    |
| scikit-learn    | 0.23.2    |
| scipy           | 1.5.2     |
| six             | 1.15.0    |
| threadpoolctl   | 2.1.0     |
| vtk             | 9.0.1     |

### R libraries  
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
  
## Anonymize your dataset
We strongly recommend to anonymize the dataset before running the pipelines. This pipeline includes a script to anonymize the T1w as well as quantitative T1 images.

# `mica-pipe` command and modules  

## MRI structural processing (`-proc_structural`)  

## fMRI processing (`-proc_rsfmri`)  
proc_rsfmri does not run time slicing correction.

## DWI processing (`-proc_dwi`)  
This step requires the `json` files in the BIDS database to contain the following variables:  
```bash
"PhaseEncodingDirection": "j-",
"TotalReadoutTime": 0.047,
```
