rule proc_dwi:
    input:
        # DWI processing requires structural output as dependency
        structural_output=lambda w: f"{output_dir}/sub-{w.subject}/ses-{w.session}/anat/processed_volumetric.nii.gz"
    output:
        processed_dwi=f"{output_dir}/sub-{{subject}}/ses-{{session}}/dwi/processed_dwi.mif"
    params:
        tmpDir="tmp",
        dwi_main=config["parameters"]["proc_dwi"].get("dwi_main", "DEFAULT"),
        dwi_rpe=config["parameters"]["proc_dwi"].get("dwi_rpe", "DEFAULT"),
        b0thr=config["parameters"]["proc_dwi"].get("b0thr", 61),
        sub=lambda w: w.subject,
        ses=lambda w: w.session
    threads: config.get("threads", 4),
    shell:
        """
        bash {script_dir}/02_proc-dwi.sh \
            {bids_dir} {params.sub} {output_dir} {params.ses} \
            --threads {threads} --tmpDir {params.tmpDir} \
            --dwi_main {params.dwi_main} --dwi_rpe {params.dwi_rpe} \
            --b0thr {params.b0thr}
        """

rule sc:
    input:
        dwi_output=lambda w: f"{output_dir}/sub-{w.subject}/ses-{w.session}/dwi/processed_dwi.mif",
        post_structural=lambda w: f"{output_dir}/sub-{w.subject}/ses-{w.session}/anat/post_structural.nii.gz"
    output:
        sc_output=f"{output_dir}/sub-{{subject}}/ses-{{session}}/connectome/sc.csv"
    params:
        tmpDir="tmp",
        sub=lambda w: w.subject,
        ses=lambda w: w.session
    threads: config.get("threads", 4),
    shell:
        """
        bash {script_dir}/03_SC.sh \
            {bids_dir} {params.sub} {output_dir} {params.ses} \
            --threads {threads} --tmpDir {params.tmpDir}
        """
