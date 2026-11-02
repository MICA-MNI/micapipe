FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash bc binutils build-essential bzip2 ca-certificates csh curl dc file \
    git libarchive-tools libeigen3-dev libfftw3-dev libfontconfig1 \
    libgl1 libglib2.0-0 libglu1-mesa libgomp1 libgsl27 libjpeg62 \
    liblapack3 libopengl0 libosmesa6 libpcre3 libquadmath0 libsm6 \
    libvtk9.1 libxext6 libxmu6 libxrender1 libxt6 locales perl \
    python3 r-base-core r-cran-e1071 r-cran-kernlab r-cran-party \
    r-cran-randomforest r-cran-rocr octave octave-signal octave-statistics \
    tar tcsh unzip wget xvfb zlib1g-dev \
    && locale-gen en_US.UTF-8 && rm -rf /var/lib/apt/lists/*
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
# Retained AFNI binaries need these historical ABI libraries.
COPY vendor/ /tmp/legacy-libs/
RUN mkdir -p /tmp/legacy-root &&     for deb in /tmp/legacy-libs/*.deb; do dpkg-deb -x "$deb" /tmp/legacy-root; done &&     cp -a /tmp/legacy-root/usr/lib/x86_64-linux-gnu/. /usr/lib/x86_64-linux-gnu/ &&     cp -a /tmp/legacy-root/lib/x86_64-linux-gnu/. /usr/lib/x86_64-linux-gnu/ &&     ldconfig && rm -rf /tmp/legacy-libs /tmp/legacy-root
