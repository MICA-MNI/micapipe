# Rule for structural processing
rule proc_structural:
    input:
        t1w=lambda wildcards: f"{bids_dir}/sub-{wildcards.subject}/ses-{wildcards.session}/anat/sub-{wildcards.subject}_ses-{wildcards.session}_{config['parameters']['proc_structural']['T1wStr']}.nii.gz",
    output:
        processed_volumetric=f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/processed_volumetric.nii.gz",
    params:
        tmpDir="tmp",
        T1wStr=config["parameters"]["proc_structural"].get("T1wStr", "T1w.nii"),
        UNI=config["parameters"]["proc_structural"].get("UNI", "FALSE"),
        MF=config["parameters"]["proc_structural"].get("MF", 3),
        subject_short=lambda wildcards: f"{wildcards.subject}",
        session_short=lambda wildcards: f"ses-{wildcards.session}",
    threads: config.get("threads", 4),
    shell:
        """
        echo "Running structural processing with subject={wildcards.subject} session={wildcards.session}, full_subject={params.subject_short}, full_session={params.session_short}"
        bash {script_dir}/01_proc-structural.sh \
            {bids_dir} {params.subject_short} {output_dir} -ses {params.session_short} \
            --threads {threads} --tmpDir {params.tmpDir} --T1wStr {params.T1wStr} --uni {params.UNI} --mf {params.MF}
        """

# Rule for cortical surface reconstruction
rule proc_surf:
    input:
        structural_output=f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/processed_volumetric.nii.gz",
    output:
        processed_surf=f"{output_dir}/sub-{{subject}}/ses-{{session}}/surf/processed_surf.gii",
    params:
        surf_dir="surf",
        freesurfer=config["parameters"]["proc_surf"].get("freesurfer", "FALSE"),
        fs_licence=config["parameters"]["proc_surf"].get("fs_licence", "path/to/fs_licence.txt"),
        sub=lambda wildcards: wildcards.subject,
        ses=lambda wildcards: wildcards.session,
    threads: config.get("threads", 4),
    shell:
        """
        bash {script_dir}/01_proc-surf.sh \
            {bids_dir} {params.sub} {output_dir} {params.ses} \
            --threads {threads} --surf_dir {params.surf_dir} --freesurfer {params.freesurfer} --fs_licence {params.fs_licence}
        """

# Rule for post structural processing
rule post_structural:
    input:
        structural_output=f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/processed_volumetric.nii.gz",
        surf_output=f"{output_dir}/sub-{{subject}}/ses-{{session}}/surf/processed_surf.gii",
    output:
        post_structural=f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/post_structural.nii.gz",
    params:
        atlas=config["parameters"]["post_structural"].get("atlas", "default"),
        freesurfer=config["parameters"]["post_structural"].get("freesurfer", "FALSE"),
        sub=lambda wildcards: wildcards.subject,
        ses=lambda wildcards: wildcards.session,
    threads: config.get("threads", 4),
    shell:
        """
        bash {script_dir}/02_post-structural.sh \
            {bids_dir} {params.sub} {output_dir} {params.ses} \
            --threads {threads} --atlas {params.atlas} --freesurfer {params.freesurfer}
        """

# Rule for geodesic distance
rule proc_geodesic_distance:
    input:
        post_structural_output=f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/post_structural.nii.gz",
    output:
        geodesic_distance=f"{output_dir}/sub-{{subject}}/ses-{{session}}/maps/geodesic_distance.nii.gz",
    params:
        sub=lambda wildcards: wildcards.subject,
        ses=lambda wildcards: wildcards.session,
    threads: config.get("threads", 4),
    shell:
        """
        bash {script_dir}/03_GD.sh \
            {bids_dir} {params.sub} {output_dir} {params.ses} \
            --threads {threads}
        """
