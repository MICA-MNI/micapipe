# configfile: "sample_run.yaml"
configfile: "config.yaml"

# Define common directories and variables
script_dir = config["script_dir"]
bids_dir = config["bids_dir"]
output_dir = config["output_dir"]

# Include any rules needed first
include: "rules/structural.smk"
include: "rules/dwi.smk"

# Define the targets for structural processing once at the top
structural_targets = [
    # processed volumetric
    expand(
        f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/processed_volumetric.nii.gz",
        subject=config["subjects"],
        session=config["sessions"]
    ),
    # processed_surf
    expand(
        "{output_dir}/sub-{subject}/ses-{session}/surf/processed_surf.gii",
        output_dir=output_dir,
        subject=config["subjects"],
        session=config["sessions"]
    ),
    # post_structural
    expand(
        "{output_dir}/sub-{subject}/ses-{session}/anat/post_structural.nii.gz",
        output_dir=output_dir,
        subject=config["subjects"],
        session=config["sessions"]
    ),
    # geodesic_distance
    expand(
        "{output_dir}/sub-{subject}/ses-{session}/maps/geodesic_distance.nii.gz",
        output_dir=output_dir,
        subject=config["subjects"],
        session=config["sessions"]
    )
]

dwi_targets =  [
       expand(f"{output_dir}/sub-{{subject}}/ses-{{session}}/dwi/processed_dwi.mif",
               subject=config["subjects"],
               session=config["sessions"]),
        expand(f"{output_dir}/sub-{{subject}}/ses-{{session}}/connectome/sc.csv",
               subject=config["subjects"],
               session=config["sessions"])
]

# The master "all" rule at the top
rule all:
    input:
       structural_targets,
       dwi_targets

# The structural rule depends on the same targets
rule structural:
    input:
       structural_targets
       
rule dwi:
    input:
       dwi_targets
