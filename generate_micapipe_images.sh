#!/bin/bash

set -e

# Generate Dockerfile.
# change ANTs version manually to 2.3.4 in the generated Dockerfile (from 2.3.1)
# change Freesurfer version manually to 7.4.0 in the generated Dockerfile (from 6.0.0)
# Manually erased the next lines from the Dockerfile see commit :
# -    && echo "Installing FSL conda environment ..." \
# -    && bash /opt/fsl-6.0.2/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.2

generate() {
  docker run --rm repronim/neurodocker:0.7.0 generate "$1" \
    --base=ubuntu:bionic-20201119 \
    --pkg-manager=apt \
    --install "gcc g++ lsb-core bsdtar jq libopenblas-dev tree openjdk-8-jdk libstdc++6" \
    --dcm2niix version=v1.0.20190902 method=source\
    --fsl version=6.0.2 \
    --run-bash 'bash /opt/fsl-6.0.2/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.2' \
    --freesurfer version=6.0.0 \
    --matlabmcr version=2017b\
    --afni version=latest\
    --ants version=2.3.1 \
    --run-bash "apt-get update && apt-get install -y gnupg2 && wget -O- http://neuro.debian.net/lists/xenial.de-fzj.full | tee /etc/apt/sources.list.d/neurodebian.sources.list && apt-key adv --recv-keys --keyserver hkps://keyserver.ubuntu.com 0xA5D32F012649A5A9 && apt-get update && apt-get install -y connectome-workbench=1.3.2-2~nd16.04+1" \
    --run-bash "cd /opt/ && wget http://www.fmrib.ox.ac.uk/~steve/ftp/fix1.068.tar.gz && tar xvfz fix1.068.tar.gz && rm fix1.068.tar.gz" \
    --user=mica \
    --miniconda \
      conda_install="python=3.9.13 certifi==2020.6.20
                     cycler==0.10.0 joblib==0.16.0
                     kiwisolver==1.2.0 matplotlib==3.4.3 nibabel==4.0.2
                     numpy==1.21.5 packaging==20.4 pandas==1.4.4
                     pillow==7.2.0 pyparsing==2.4.7 python-dateutil==2.8.1
                     pytz==2020.1 scikit-learn==1.0.2 scipy==1.9.1
                     six==1.15.0 threadpoolctl==2.1.0 vtk==9.2.2"\
      pip_install='brainspace==0.1.4 seaborn==0.11.2 tedana==0.0.12 pygeodesic==0.1.8 xhtml2pdf==0.2.9' \
      pip_install='sys os glob argparse==1.1 copy time re warnings subprocess' \
      create_env="micapipe" \
      activate=true \
    --run-bash 'source activate micapipe && conda install -c mrtrix3 mrtrix3==3.0.1 && pip install git+https://github.com/MICA-MNI/ENIGMA.git' \
    --run-bash 'git clone https://github.com/Deep-MI/FastSurfer.git && mv FastSurfer /opt/ && conda env create -f /opt/FastSurfer/fastsurfer_env_cpu.yml' \
    --run-bash 'source activate fastsurfer_cpu && python /opt/FastSurfer/FastSurferCNN/download_checkpoints.py --all && source deactivate' \
    --user=root\
    --run "set -uex; \
           LD_LIBRARY_PATH=/lib64/:$PATH; \
           apt install -y software-properties-common apt-transport-https; \
           apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9; \
           add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'; \
           apt update; \
           apt install -y r-base libblas-dev liblapack-dev gfortran g++ libgl1-mesa-glx; \
           rm -rf /var/lib/apt/lists/*;" \
    --run-bash "wget https://www.dropbox.com/s/47lu1nojrderls1/install_R_env.sh?dl=0 -O /opt/install_R_env.sh &&
                bash /opt/install_R_env.sh && cd /opt/afni-latest && rPkgsInstall -pkgs ALL" \
    --copy . /opt/micapipe \
    --run-bash "cd /opt/micapipe && mv fix_settings.sh /opt/fix1.068/settings.sh && mv fsl_conf/* /opt/fsl-6.0.2/etc/flirtsch/" \
    --run-bash "mv /opt/micapipe/surfaces/fsaverage5 /opt/freesurfer-6.0.0/subjects" \
    --workdir='/home/mica' \
    --env MICAPIPE='/opt/micapipe'\
    --env PROC='container-micapipe v0.2.0' \
    --add-to-entrypoint "source /opt/freesurfer-6.0.0/SetUpFreeSurfer.sh" \
    --add-to-entrypoint "export FIXPATH=/opt/fix && export PATH="${FIXPATH}:${PATH}"" \
    --entrypoint "/neurodocker/startup.sh /opt/micapipe/micapipe"
  }


generate docker > Dockerfile
