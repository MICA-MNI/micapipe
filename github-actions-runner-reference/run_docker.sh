export DOCKER_CONTENT_TRUST=0

docker run -d --name micapipe-runner --restart always --privileged  \
    -v /data/mica1/03_projects/enning/BIDS_CI/sing_cache:/data/mica1/03_projects/enning/BIDS_CI/sing_cache \
    -v /export02/local/singularity_tmp/:/export02/local/singularity_tmp/ \
    -v /data/mica3/BIDS_CI/rawdata:/data/mica3/BIDS_CI/rawdata \
    -v /data_/mica3/BIDS_CI:/data_/mica3/BIDS_CI \
    -v /data/mica1/03_projects/enning/BIDS_CI:/data/mica1/03_projects/enning/BIDS_CI \
    -v /tmp:/tmp \
    -v /var/run/docker.sock:/var/run/docker.sock \
    micapipe-runner