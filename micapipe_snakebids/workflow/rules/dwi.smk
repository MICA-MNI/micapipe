def get_dwi_outputs(inputs, output_dir):
    return bids(
        root=f"{output_dir}/micapipe_v0.2.0",
        datatype="dwi",
        space="dwi",
        desc="preproc",
        suffix="dwi.mif",
        subject="{subject}",
        session="{session}"
    ) 

def get_sc_outputs(inputs, output_dir):
    tracts = config['parameters']['SC']['tracts'].replace("M", "000000")
    return bids(
        root=f"{output_dir}/micapipe_v0.2.0",
        datatype="dwi",
        space="dwi",
        desc=f"iFOD2-{tracts}", #TODO: check this
        suffix="tdi.nii.gz",
        subject="{subject}",
        session="{session}"
    ) 

# rule for diffusion processing
rule proc_dwi:
    input:
        # DWI processing requires structural output as dependency
        inputs['t1w'].expand(
            get_structural_outputs(inputs, output_dir)
        ),
        inputs['dwi'].expand()
    output:
        processed_dwi=get_dwi_outputs(inputs, output_dir)
    params:
        dwi_main=process_multi_inputs(config["parameters"]["proc_dwi"]["dwi_main"]),
        dwi_rpe=process_multi_inputs(config["parameters"]["proc_dwi"]["dwi_rpe"]),
        dwi_processed=config["parameters"]["proc_dwi"]["dwi_processed"],
        rpe_all=config["parameters"]["proc_dwi"]["rpe_all"],
        b0thr=config["parameters"]["proc_dwi"]["b0thr"],
        # the following are just flags
        dwi_acq=process_optional_flags(config["parameters"]["proc_dwi"]["dwi_acq"], "dwi_acq"),
        no_bvalue_scaling=process_flags(config["parameters"]["proc_dwi"]["no_bvalue_scaling"], "no_bvalue_scaling"),
        dwi_upsample=process_flags(config["parameters"]["proc_dwi"]["dwi_upsample"], "dwi_upsample"),
    threads: config.get("threads", 4),
    shell:
        """
        {command} -sub sub-{wildcards.subject} -out {output_args} -bids {bids_args} -proc_dwi \
            -threads {threads} -ses {wildcards.session} -dwi_main {params.dwi_main} -dwi_rpe {params.dwi_rpe} \
            -dwi_processed {params.dwi_processed} -rpe_all {params.rpe_all} \
            -b0thr {params.b0thr} {params.dwi_acq} {params.no_bvalue_scaling} {params.dwi_upsample}
        """

rule sc:
    input:
        # inputs['t1w'].expand(
        #     get_structural_outputs(inputs, output_dir)
        # ),
        inputs['t1w'].expand(
            get_surf_outputs(inputs, output_dir)
        ),
        inputs['t1w'].expand(
            get_post_structural_outputs(inputs, output_dir)
        ),
        inputs['dwi'].expand(
            get_dwi_outputs(inputs, output_dir)
        ),
    
    output:
        get_sc_outputs(inputs, output_dir)
    params:
        tracts=config["parameters"]["SC"]["tracts"],
        keep_tck=config["parameters"]["SC"]["keep_tck"],
        autoTract=config["parameters"]["SC"]["autoTract"],
        dwi_acq=process_optional_flags(config["parameters"]["proc_dwi"]["dwi_acq"], "dwi_acq"),
        weighted_SC=config["parameters"]["SC"]["weighted_SC"],
        #flags
        tck=process_flags(config["parameters"]["SC"]["tck"], "tck " + config["parameters"]["SC"]["tck"]),
    threads: config.get("threads", 4),
    shell:
        """
        {command} -sub sub-{wildcards.subject} -out {output_args} -bids {bids_args} -SC \
            -threads {threads} -ses {wildcards.session} -tracts {params.tracts} -keep_tck {params.keep_tck} \
            -autoTract {params.autoTract} -weighted_SC {params.weighted_SC} {params.dwi_acq} {params.tck}
        """
