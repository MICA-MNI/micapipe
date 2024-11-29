configfile: "config.yaml"

include: "rules/structural.smk"
# include: "rules/diffusion.smk"

rule all:
    input:
        expand("{output_dir}/{subject}/{session}/anat/processed_volumetric.nii.gz",
               output_dir=config["output_dir"],
               subject=config["subjects"],
               session=config["sessions"])
