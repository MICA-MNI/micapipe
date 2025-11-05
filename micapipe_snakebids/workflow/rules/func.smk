def get_func_outputs(inputs, output_dir):
    return bids(
        root=f"{output_dir}/micapipe_v0.2.0",
        datatype="func/desc-se_task-{task}_acq-{acq}_bold/volumetric",
        space="func",
        desc="se",
        suffix="preproc.nii.gz",
        subject="{subject}",
        session="{session}"
    )

def get_func_inputs(inputs, output_dir):
    return bids(
        root=f"{output_dir}/micapipe_v0.2.0",
        datatype="func",
        suffix="bold.nii.gz",
        **inputs["func"].wildcards
    )

rule proc_func:
    input:
        inputs['func'].expand(),
        inputs['t1w'].expand(
            get_structural_outputs(inputs, output_dir)
        ),
        inputs['t1w'].expand(
            get_surf_outputs(inputs, output_dir)
        ),
        inputs['t1w'].expand(
            get_post_structural_outputs(inputs, output_dir)
        )
    output:
        processed_func=get_func_outputs(inputs, output_dir)
    params:
        mainScanStr=config["parameters"]["proc_func"]["mainScanStr"],
        func_pe=config["parameters"]["proc_func"]["func_pe"],
        func_rpe=config["parameters"]["proc_func"]["func_rpe"],

        # the following are just optional flags, default is False otherwise the expected value
        mainScanRun = process_optional_flags(config["parameters"]["proc_func"]["mainScanRun"], "mainScanRun"),
        phaseReversalRun=process_optional_flags(config["parameters"]["proc_func"]["phaseReversalRun"], "phaseReversalRun"),
        topupConfig=process_optional_flags(config["parameters"]["proc_func"]["topupConfig"], "topupConfig"),
        icafixTraining=process_optional_flags(config["parameters"]["proc_func"]["icafixTraining"], "icafixTraining"),

        # the following are just flags
        # Boolean flags that only appear if set to TRUE
        smoothWithWB=process_flags(config["parameters"]["proc_func"]["smoothWithWB"], "smoothWithWB"),
        NSR=process_flags(config["parameters"]["proc_func"]["NSR"], "NSR"),
        GSR=process_flags(config["parameters"]["proc_func"]["GSR"], "GSR"),
        noFIX=process_flags(config["parameters"]["proc_func"]["noFIX"], "noFIX"),
        dropTR=process_flags(config["parameters"]["proc_func"]["dropTR"], "dropTR"),
        noFC=process_flags(config["parameters"]["proc_func"]["noFC"], "noFC"),

    threads: config.get("threads", 4)
    shell:
        """
        {command} -sub sub-{wildcards.subject} -out {output_args} -bids {bids_args} -proc_func \
            -threads {threads} -ses {wildcards.session} -mainScanStr {params.mainScanStr} -func_pe {params.func_pe} \
            -func_rpe {params.func_rpe} {params.mainScanRun} {params.phaseReversalRun} {params.topupConfig} \
            {params.icafixTraining} {params.smoothWithWB} {params.NSR} {params.GSR} {params.noFIX} \
            {params.dropTR} {params.noFC} 
        """