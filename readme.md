![micapipe logo](/docs/figures/micapipe_long.png)

## Multimodal connectome processing with the `micapipe` ##

[![License: GPL v3](https://img.shields.io/github/license/MICA-MNI/micapipe)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/micapipe/badge/?version=latest)](https://micapipe.readthedocs.io/en/latest/?badge=latest)
[![CircleCI](https://circleci.com/gh/MICA-MNI/micapipe/tree/master.svg?style=shield)](https://circleci.com/gh/rcruces/MICA-MNI/tree/master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7da8a9a3524745bc9616fd465a17f41b)](https://app.codacy.com/gh/rcruces/micapipe?utm_source=github.com&utm_medium=referral&utm_content=rcruces/micapipe&utm_campaign=Badge_Grade)
[![GitHub stars](https://img.shields.io/github/stars/MICA-MNI/micapipe)](https://github.com/MICA-MNI/micapipe/stargazers)
[![GitHub issues](https://img.shields.io/github/issues/MICA-MNI/micapipe)](https://github.com/MICA-MNI/micapipe/issues)
[![version](https://img.shields.io/github/v/tag/MICA-MNI/micapipe)](https://github.com/MICA-MNI/micapipe)
[![Docker Pulls](https://img.shields.io/docker/pulls/micalab/micapipe)](https://hub.docker.com/r/micalab/micapipe)
![Docker Image Version (latest semver)](https://img.shields.io/docker/v/micalab/micapipe?color=orange&label=docker%20version)

[`micapipe`](micapipe.readthedocs.io) is developed by [MICA-lab](https://mica-mni.github.io) at McGill University for use at [the Neuro](https://www.mcgill.ca/neuro/), McConnell Brain Imaging Center ([BIC](https://www.mcgill.ca/bic/)).  
> The main goal of this pipeline is to provide a semi-flexible and robust framework to process MRI images and generate ready to use modality based connectomes.    
> The `micapipe` utilizes a set of known software dependencies, different brain atlases, and software developed in our laboratory.
> The basic cutting edge processing of our pipelines aims the *T1 weighted images*, *resting state fMRI* and *Diffusion weighted images*.

## Documentation ##
You can find the documentation in [micapipe.readthedocs.io](http://micapipe.readthedocs.io/en/latest/)

## Container ##
You can find the latest version of the container in [Docker](https://hub.docker.com/r/micalab/micapipe/)

## Advantages ##
-   Microstructure Profile Covariance ([Paquola C et al. Plos Biology 2019](https://doi.org/10.1371/journal.pbio.3000284)).  
-   Multiple parcellations (18 x 3).  
-   Includes cerebellum and subcortical areas.  
-   Surface based analysis.  
-   Latest version of software dependencies.  
-   Ready to use outputs.  
-   Easy to use.  
-   Standardized format (BIDS).  

## How to cite micapipe ##
> Raúl R. Cruces, Jessica Royer, Peer Herholz, Sara Larivière, Reinder Vos de Wael, Casey Paquola, Oualid Benkarim, Bo-yong Park, Janie Degré-Pelletier, Mark Nelson, Jordan DeKraker, Christine Tardif, Jean-Baptiste Poline, Luis Concha, Boris C. Bernhardt. (2022). *Micapipe: A Pipeline for Multimodal Neuroimaging and Connectome Analysis*. bioRxiv 2022.01.31.478189. doi: https://doi.org/10.1101/2022.01.31.478189

## Dependencies ##
| *Software*   |     *Version*   | *Further info* |
|------------|---------------|--------------|  
| dcm2niix   | v1.0.20190902 | https://github.com/rordenlab/dcm2niix |
| Freesurfer | 6.0.0         | https://surfer.nmr.mgh.harvard.edu/ |
| FSl        | 6.0.3           | https://fsl.fmrib.ox.ac.uk/fsl/fslwiki |
| AFNI       | 20.3.03       | https://afni.nimh.nih.gov/download |
| MRtrix3    | 3.0.1         | https://www.mrtrix.org |
| ANTs       | 2.3.4         | https://github.com/ANTsX/ANTs |
| workbench  | 1.4.2         | https://www.humanconnectome.org/software/connectome-workbench |
| FIX        | 1.06          | https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX |
| R          | 3.6.3         | https://www.r-project.org |
| python     | 3.7.6         | https://www.python.org/downloads/ |

 > The FIX package (FMRIB's ICA-based Xnoiseifier) requires FSL, R and one of MATLAB Runtime Component, full MATLAB or Octave. We recommend the use of the MATLAB Runtime Component. Additionally, it requires the following R libraries:  'kernlab','ROCR','class','party','e1071','randomForest'

### `python` packages ###
|     *Package*     |  *Version*  |
|:---------------:|:---------:|
| brainspace      | 0.1.1     |
| certifi         | 2020.6.20 |
| cycler          | 0.10.0    |
| argparse        | 0.1.1     |
| joblib          | 0.16.0    |
| kiwisolver      | 1.2.0     |
| matplotlib      | 3.3.1     |
| nibabel         | 3.1.1     |
| nilearn         | 0.6.2     |
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

### `R` libraries ###
| *Core   library* |  *version* | *Dependency library* |  *version* | *Dependency library* |   *version* |
|:--------------:|:--------:|:------------------:|:--------:|:------------------:|:---------:|
|         plotly |  4.9.2.1 |         tidyselect |  1.1.0   |           parallel |  3.6.3    |
|        viridis |  0.5.1   |               coin |  1.3-1   |            TH.data |  1.0-10   |
|    viridisLite |  0.3.0   |              purrr |  0.3.4   |               Rcpp |  1.0.5    |
|          tidyr |  1.1.2   |            splines |  3.6.3   |           jsonlite |  1.6.1    |
|        ggplot2 |  3.3.2   |            lattice |  0.20-41 |          gridExtra | 2.3       |
|         scales |  1.1.1   |         colorspace |  1.4-1   |             digest |  0.6.27   |
|   randomForest |  4.6-14  |              vctrs |  0.3.5   |              dplyr |  1.0.2    |
|          e1071 |  1.7-4   |           generics |  0.0.2   |              tools |  3.6.3    |
|          party |  1.3-5   |          htmltools |  0.4.0   |           magrittr |  2.0.1    |
|    strucchange |  1.5-2   |           survival |  3.1-12  |           lazyeval |  0.2.2    |
|       sandwich |  2.5-1   |              rlang |  0.4.9   |             tibble |  3.0.4    |
|            zoo |  1.8-7   |             pillar |  1.4.3   |             crayon |  1.3.4    |
|     modeltools |  0.2-23  |               glue |  1.4.2   |          pkgconfig |  2.0.3    |
|        mvtnorm |  1.1-1   |              withr |  2.3.0   |               MASS |  7.3-51.5 |
|          class |  7.3-17  |        matrixStats |  0.56.0  |           ellipsis |  0.3.0    |
|           ROCR |  1.0-11  |           multcomp |  1.4-13  |            libcoin |  1.0-6    |
|        kernlab |  0.9-29  |          lifecycle |  0.2.0   |             Matrix |  1.2-18   |
|      networkD3 |    0.4   |            munsell |  0.5.0   |         data.table |  1.12.8   |
|                |          |             gtable |  0.3.0   |               httr |  1.4.1    |
|                |          |        htmlwidgets |  1.5.1   |                 R6 |  2.4.1    |
|                |          |          codetools |  0.2-16  |           compiler |  3.6.3    |
