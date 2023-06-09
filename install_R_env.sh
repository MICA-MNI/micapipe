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

su - -c "R -e \"install.packages('versions')\""
su - -c "R -e \"require(versions) | install.versions('scales', '1.1.1')\""
su - -c "R -e \"require(versions) | install.versions('randomForest', '4.6-14')\""
su - -c "R -e \"require(versions) | install.versions('e1071', '1.7-4')\""
su - -c "R -e \"require(versions) | install.versions('party', '1.3-5')\""
su - -c "R -e \"require(versions) | install.versions('strucchange', '1.5-2')\""
su - -c "R -e \"require(versions) | install.versions('sandwich', '2.5-1')\""
su - -c "R -e \"require(versions) | install.versions('zoo', '1.8-7')\""
su - -c "R -e \"require(versions) | install.versions('modeltools', '0.2-23')\""
su - -c "R -e \"require(versions) | install.versions('mvtnorm', '1.1-1')\""
su - -c "R -e \"require(versions) | install.versions('class', '7.3-17')\""
su - -c "R -e \"require(versions) | install.versions('ROCR', '1.0-11')\""
su - -c "R -e \"require(versions) | install.versions('kernlab', '0.9-29')\""
su - -c "R -e \"require(versions) | install.versions('coin', '1.3-1')\""
su - -c "R -e \"require(versions) | install.versions('pkgconfig', '2.0.3')\""
su - -c "R -e \"require(versions) | install.versions('MASS', '7.3-51.5')\""
su - -c "R -e \"require(versions) | install.versions('libcoin', 'libcoin')\""
su - -c "R -e \"require(versions) | install.versions('Matrix', '1.2-18')\""
