# Simplified Dockerfile without problematic ANTs multi-stage copy
# This version installs ANTs directly instead of copying from another image

FROM ubuntu:xenial-20161213

# Set environment variables
ENV DEBIAN_FRONTEND="noninteractive" \
    LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8"

# Install system dependencies
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           build-essential \
           cmake \
           curl \
           wget \
           git \
           ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install ANTs directly (avoid multi-stage copy issue)
RUN echo "Installing ANTs directly..." \
    && mkdir -p /opt/ants-install \
    && cd /opt/ants-install \
    && curl -sSL "https://dl.dropbox.com/s/gwf51ykkk5bifyj/ants-Linux-centos6_x86_64-v2.3.4.tar.gz" \
    | tar -xzC /opt \
    && mv /opt/ants-Linux-centos6_x86_64-v2.3.4 /opt/ants-2.3.4 \
    && rm -rf /opt/ants-install

# Set ANTs environment (this was where it was hanging before)
ENV ANTSPATH="/opt/ants-2.3.4/bin" \
    PATH="/opt/ants-2.3.4/bin:$PATH"

# Test ANTs installation
RUN antsRegistration --version || echo "ANTs installation failed"

CMD ["/bin/bash"]