rule proc_structural:
    input:
        t1w=lambda wildcards: f"{config['bids_dir']}/{wildcards.subject}/{wildcards.session}/anat/{wildcards.subject}_{wildcards.session}_{config['parameters']['proc_structural']['T1wStr']}"
    output:
        volumetric = f"{config['output_dir']}/{{subject}}/{{session}}/anat/processed_volumetric.nii.gz"
    params:
        threads = config["threads"],
        UNI = config["parameters"]["proc_structural"]["UNI"],
        MF = config["parameters"]["proc_structural"]["MF"]
    shell:
        """
        echo Your bash command is being constructed
        """
