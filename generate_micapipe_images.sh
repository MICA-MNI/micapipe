#!/bin/bash

set -e

# Generate Dockerfile.
# change ANTs version manually to 2.3.4 in the generated Dockerfile (from 2.3.1)
# change Freesurfer version manually to 7.4.0 in the generated Dockerfile (from 6.0.0)
# Manually erased the next lines from the Dockerfile see commit :
# -    && echo "Installing FSL conda environment ..." \
# -    && bash /opt/fsl-6.0.3/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.3
# workbench    --run-bash "apt-get update && apt-get install -y gnupg2 && wget -O- http://neuro.debian.net/lists/xenial.de-fzj.full | tee /etc/apt/sources.list.d/neurodebian.sources.list && apt-key adv --recv-keys --keyserver hkps://keyserver.ubuntu.com 0xA5D32F012649A5A9 && apt-get update && apt-get install -y connectome-workbench=1.3.2-2~nd16.04+1" \

generate() {
  docker run --rm repronim/neurodocker:0.7.0 generate "$1" \
    --base=ubuntu:bionic-20201119 \
    --pkg-manager=apt \
    --install "gcc g++ lsb-core bsdtar jq libopenblas-dev tree openjdk-8-jdk libstdc++6" \
    --dcm2niix version=v1.0.20190902 method=source\
    --fsl version=6.0.3 \
    --run-bash 'bash /opt/fsl-6.0.3/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.3' \
    --freesurfer version=6.0.0 \
    --matlabmcr version=2017b\
    --afni version=latest\
    --ants version=2.3.1 \
    --install connectome-workbench \
    --run-bash "cd /opt/ && wget http://www.fmrib.ox.ac.uk/~steve/ftp/fix1.068.tar.gz && tar xvfz fix1.068.tar.gz && rm fix1.068.tar.gz" \
    --user=mica \
    --miniconda \
      version=22.11.1 \
      conda_install="python=3.9.13 aiohttp==3.8.4 aiosignal==1.3.1
                    asn1crypto==1.5.1 async-timeout==4.0.2 attrs==22.2.0
                    bokeh==2.2.3  cffi==1.15.1 charset-normalizer==3.1.0 click==8.1.3 contourpy==1.0.7
                    cryptography==39.0.2 cycler==0.11.0
                    fonttools==4.39.2 frozenlist==1.3.3 html5lib==1.1
                    idna==3.4 importlib-resources==5.12.0 jinja2==3.0.1
                    joblib==1.2.0 kiwisolver==1.4.4 lxml==4.9.2
                    markupsafe==2.1.2 matplotlib==3.4.3
                    multidict==6.0.4 nibabel==4.0.2 nilearn==0.10.0
                    numpy==1.21.5 packaging==23.0
                    pandas==1.4.4 pillow==9.4.0 pycparser==2.21
                    pyhanko-certvalidator==0.20.1
                    pyparsing==3.0.9 pypdf==3.6.0 pypng==0.20220715.0
                    python-bidi==0.4.2 python-dateutil==2.8.2 pytz==2022.7.1
                    pytz-deprecation-shim==0.1.0.post0 pyyaml==6.0 qrcode==7.4.2
                    reportlab==3.6.12 requests==2.28.2 scikit-learn==1.0.2
                    scipy==1.9.1 seaborn==0.11.2 six==1.16.0
                    svglib==1.5.1 threadpoolctl==3.1.0
                    tinycss2==1.2.1 tornado==6.2 typing-extensions==4.5.0
                    tzlocal==4.3 uritools==4.0.1 urllib3==1.26.15
                    vtk==9.2.2 webencodings==0.5.1 wslink==1.10.1 yarl==1.8.2 zipp==3.15.0 " \
      pip_install='argparse==1.1 brainspace==0.1.4 tedana==0.0.12 pyhanko==0.17.2 mapca==0.0.3
                    xhtml2pdf==0.2.9 oscrypto==1.3.0 tzdata==2022.7 arabic-reshaper==3.0.0
                    cssselect2==0.7.0 pygeodesic==0.1.8' \
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
    --run-bash "wget https://sourceforge.net/projects/c3d/files/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz/download -O itksnap.tar.gz &&
                tar -xfv itksnap.tar.gz -C /opt/" \
    --env PATH="/opt/itksnap/bin/:${PATH}" \
    --run-bash "wget https://www.dropbox.com/s/47lu1nojrderls1/install_R_env.sh?dl=0 -O /opt/install_R_env.sh &&
                bash /opt/install_R_env.sh && cd /opt/afni-latest && rPkgsInstall -pkgs ALL" \
    --copy . /opt/micapipe \
    --run-bash "cd /opt/micapipe && mv fix_settings.sh /opt/fix1.068/settings.sh && mv fsl_conf/* /opt/fsl-6.0.3/etc/flirtsch/" \
    --run-bash "mv /opt/micapipe/surfaces/fsaverage5 /opt/freesurfer-6.0.0/subjects" \
    --workdir='/home/mica' \
    --env MICAPIPE='/opt/micapipe'\
    --env PROC='container-micapipe v0.2.0' \
    --add-to-entrypoint "export FIXPATH=/opt/fix && export PATH="${FIXPATH}:${PATH}"" \
    --entrypoint "/neurodocker/startup.sh /opt/micapipe/micapipe"
  }


generate docker > Dockerfile

echo -e "###########################################################################################\n
NOTES:
> change ANTs version manually to 2.3.4 in the generated Dockerfile (from 2.3.1)

> change Freesurfer version manually to 7.4.0 in the generated Dockerfile (from 6.0.0)
    REPLACE: surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.0/freesurfer-Linux-centos6_x86_64-stable-pub-v7.4.0.tar.gz
    with surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.0/freesurfer-linux-ubuntu18_amd64-7.4.0.tar.gz

> Manually erased the next lines from the Dockerfile see commit :
-  && echo Installing FSL conda environment ... \
-  && bash /opt/fsl-6.0.3/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.3

-  RUN bash -c 'bash /opt/fsl-6.0.3/etc/fslconf/fslpython_install.sh -f /opt/fsl-6.0.3'

> libxp6 issue
-  replaCE : http://mirrors.kernel.org/debian/pool/main/libx/libxp/libxp6_1.0.2-2_amd64.deb
-      WITH: http://launchpadlibrarian.net/160108232/libxp6_1.0.2-1ubuntu1_amd64.deb

> Miniconda
- https://repo.continuum.io/miniconda/Miniconda3-22.11.1-Linux-x86_64.sh
- https://repo.anaconda.com/miniconda/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh

###########################################################################################\n"
