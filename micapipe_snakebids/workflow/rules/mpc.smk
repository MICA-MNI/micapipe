
def get_mpc_outputs(inputs, output_dir):
    return bids(
        root=f"{output_dir}/micapipe_v0.2.0",
        datatype="xfm",     
        subject="{subject}",
        session="{session}",
        suffix="from-T1map_to-fsnative_0GenericAffine.mat"     # TODO: check how to get acq?
    )



rule proc_mpc:
    # A) Inputs
    input:
        inputs["t1w"].expand(
            get_post_structural_outputs(inputs, output_dir)
        ),
    output:
        processed_mpc = get_mpc_outputs(inputs, output_dir)
    params:
        microstructural_img = process_optional_flags(
            config["parameters"]["proc_mpc"]["microstructural_img"],
            "microstructural_img"
        ),
        mpc_acq = process_optional_flags(
            config["parameters"]["proc_mpc"]["mpc_acq"],
            "mpc_acq"
        ),

        # Boolean flags that only appear if set to TRUE
    threads: config.get("threads", 4)
    shell:
        """
        {command} -sub sub-{wildcards.subject} -out {output_args} -bids {bids_args} -MPC \
            -threads {threads} -ses {wildcards.session} \
            {params.microstructural_img} \
            {params.mpc_acq}
        """

