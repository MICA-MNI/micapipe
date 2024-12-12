# configfile: "config.yaml"

# include: "rules/structural.smk"
#
# rule all:
#     input:
#         expand("{output_dir}/sub-{subject}/ses-{session}/anat/processed_volumetric.nii.gz",
#                output_dir=config["output_dir"],
#                subject=config["subjects"],
#                session=config["sessions"]),
#         expand("{output_dir}/sub-{subject}/ses-{session}/surf/processed_surf.gii",
#                output_dir=config["output_dir"],
#                subject=config["subjects"],
#                session=config["sessions"]),
#         expand("{output_dir}/sub-{subject}/ses-{session}/anat/post_structural.nii.gz",
#                output_dir=config["output_dir"],
#                subject=config["subjects"],
#                session=config["sessions"]),
#         expand("{output_dir}/sub-{subject}/ses-{session}/maps/geodesic_distance.nii.gz",
#                output_dir=config["output_dir"],
#                subject=config["subjects"],
#                session=config["sessions"]),
configfile: "sample_run.yaml"

# Define common directories for the workflow
script_dir = config["script_dir"]
bids_dir = config["bids_dir"]
output_dir = config["output_dir"]

# Include the structural rules
include: "rules/structural.smk"

# Default target rule
rule all:
    input:
        # Expand processed volumetric outputs for all subjects and sessions
        expand(
            f"{output_dir}/sub-{{subject}}/ses-{{session}}/anat/processed_volumetric.nii.gz",
            subject=config["subjects"],
            session=config["sessions"]
        ),
        expand("{output_dir}/sub-{subject}/ses-{session}/surf/processed_surf.gii",
               output_dir=config["output_dir"],
               subject=config["subjects"],
               session=config["sessions"]),
        expand("{output_dir}/sub-{subject}/ses-{session}/anat/post_structural.nii.gz",
               output_dir=config["output_dir"],
               subject=config["subjects"],
               session=config["sessions"]),
        expand("{output_dir}/sub-{subject}/ses-{session}/maps/geodesic_distance.nii.gz",
               output_dir=config["output_dir"],
               subject=config["subjects"],
               session=config["sessions"]),
