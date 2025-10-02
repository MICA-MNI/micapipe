FROM ubuntu:20.04

# Test the COPY instruction
COPY downloads/ /downloads/

# Debug what got copied
RUN echo "Testing COPY instruction..." \
    && echo "Contents of /downloads:" \
    && ls -la /downloads/ \
    && echo "Checking for FSL file:" \
    && test -f /downloads/fsl-6.0.2-centos6_64.tar.gz && echo "FSL file found" || echo "FSL file NOT found" \
    && echo "Checking for FreeSurfer file:" \
    && test -f /downloads/freesurfer-linux-centos6_x86_64-7.4.1.tar.gz && echo "FreeSurfer file found" || echo "FreeSurfer file NOT found"