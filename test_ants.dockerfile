# Test Dockerfile to isolate ANTs installation issue
FROM ubuntu:xenial-20161213

# Test just the ANTs copy and environment setup
COPY --from=kaczmarj/ants:2.3.4 /opt/ants /opt/ants-2.3.4

# This is where it's hanging - test this step
ENV ANTSPATH="/opt/ants-2.3.4/bin" \
    PATH="/opt/ants-2.3.4/bin:$PATH"

# Verify installation
RUN ls -la /opt/ants-2.3.4/bin/ || echo "ANTs not found"

CMD ["/bin/bash"]