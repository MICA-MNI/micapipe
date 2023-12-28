#!/bin/bash
#
#apt-get install -y apt-transport-https software-properties-common build-essential
#sh -c 'echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
#gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
#gpg -a --export E084DAB9 | apt-key add -
#add-apt-repository -r 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' -y
#apt-get install -y apt-transport-https software-properties-common build-essential
# apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
# add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
# apt-get update -y
# apt-get -y install r-base libapparmor1 libcurl4-gnutls-dev libxml2-dev libssl-dev gdebi-core libcairo2-dev libxt-dev
R -e 'install.packages("/opt/randomForest_4.6-14.tar.gz", repos = NULL, type = "source")'
R -e 'install.packages("/opt/MASS_7.3-54.tar.gz", repos = NULL, type = "source")'
R -e 'install.packages("/opt/class_7.3-19.tar.gz", repos = NULL, type = "source")'
R -e 'install.packages("/opt/kernlab_0.9-32.tar.gz", repos = NULL, type = "source")'

R -e "install.packages('versions')"
R -e "install.packages('Matrix', version = '1.2-18')"
R -e "install.packages('scales', version = '1.1.1')"
R -e "install.packages('strucchange', version = '1.5-2')"
R -e "install.packages('sandwich', version = '2.5-1')"
R -e "install.packages('zoo', version = '1.8-7')"
R -e "install.packages('modeltools', version = '0.2-23')"
R -e "install.packages('mvtnorm', version = '1.1-1')"
R -e "install.packages('coin', version = '1.3-1')"
R -e "install.packages('pkgconfig', version = '2.0.3')"
R -e "install.packages('libcoin')"

# Packages for fsl-fix
R -e "install.packages('ROCR', version = '1.0-11')"
R -e "install.packages('party', version = '1.3-13')"
R -e "install.packages('e1071', version = '1.7-4')"
