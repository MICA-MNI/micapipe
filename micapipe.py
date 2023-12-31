import argparse
import subprocess
import shlex
import sys


# Get all arguments in a BIDS-Apps like manner
def parse_arguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Python script to translate BIDS App inputs to MICAPIPE inputs.')

    # Define arguments similar to BIDS-Apps, mandatory ones
    parser.add_argument('input_dir', help='BIDS input directory')
    parser.add_argument('output_dir', help='BIDS output directory')
    parser.add_argument('analysis_level', help='Analysis level')
    parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label corresponds to sub-<participant_label> from'
                                                    'the BIDS spec (so it does not include "sub-"). If this parameter is not provided all subjects should be analyzed.'
                                                    'Multiple participants can be specified with a space separated list.')
    parser.add_argument('--session', help='The label(s) of the session(s) that should be analyzed. The label corresponds to ses-<participant_label> from'
                                          'the BIDS spec (so it does not include "ses-"). If this parameter is not provided all sessions should be analyzed.'
                                          'Multiple sessions can be specified with a space separated list.')

    # Define arguments similar to BIDS-Apps, optional ones
    parser.add_argument('-proc_structural', action='store_true', help='Basic volumetric processing on T1-weighted data')
    parser.add_argument('-proc_surf', action='store_true', help="Run freesurfer's recon-all pipeline on T1-weighted data")
    parser.add_argument('-post_structural', action='store_true', help='Further structural processing after quality control')
    parser.add_argument('-GD', action='store_true', help='Generate geodesic distance matrices')
    parser.add_argument('-MPC', action='store_true', help='Equivolumetric surface mapping and computation of MPC matrices')
    parser.add_argument('-MPC_SWM', action='store_true', help='Superficial white matter mapping and computation of MPC matrices')
    parser.add_argument('-SWM', action='store_true', help='Superficial white matter processing')
    parser.add_argument('-proc_flair', action='store_true', help='T2/FLAIR processing')
    parser.add_argument('-proc_dwi', action='store_true', help='Basic diffusion-weighted imaging processing')
    parser.add_argument('-SC', action='store_true', help='Diffusion tractography and structural connectomes generation')
    parser.add_argument('-proc_func', action='store_true', help='Resting-state functional processing and functional connectomes generation')
    parser.add_argument('-QC_subj', action='store_true', help='Creates an individual report of the different modules')
    parser.add_argument('-QC', action='store_true', help='Creates a group-level pdf of the subjects with QC')

    # Define arguments similar to BIDS-Apps, additional optional ones
    parser.add_argument('--h', '-help', action='store_true', help='Print help message')
    parser.add_argument('--v', '-version', action='store_true', help='Print software version')
    parser.add_argument('-force', action='store_true', help='Overwrite existing data')
    parser.add_argument('-quiet', action='store_true', help='Do not print comments and warnings')
    parser.add_argument('-nocleanup', action='store_true', help='Prevent deletion of temporary directory')
    parser.add_argument('-threads', type=int, default=6, help='Number of threads (default=6)')
    parser.add_argument('-tmpDir', help='Custom location for temporary directory')
    parser.add_argument('-regSynth', action='store_true', help='Perform registration based on synthseg')

    # Parse arguments
    args = parser.parse_args()
    return args


def check_mandatory_arguments(args):
    missing_args = []
    if not args.input_dir:
        missing_args.append('input')
    if not args.output_dir:
        missing_args.append('output')
    if not args.analysis_level:
        missing_args.append('analysis level')
    if args.analysis_level == 'participant':
        if args.participant_label is None:
            print("Analysis level is set to 'participant', but no participant label is provided.")
            missing_args.append('participant_label')
    if not args.session:
        missing_args.append('session')

    if missing_args:
        print(f"Missing mandatory argument(s): {', '.join(missing_args)}. Please set these arguments.")
        sys.exit(1)


def construct_micapipe_command(args):
    # Start constructing the command
    command = "micapipe"

    # Add all arguments to the command
    for arg, value in vars(args).items():
        if value and arg not in ['input_dir', 'output_dir', 'analysis_level', 'participant_label', 'session', 'h', 'help', 'v', 'version']:
            command += f" -bids {args.input_dir} -out {args.output_dir} -sub {args.participant_label} -ses {args.session}"
            if isinstance(value, bool):
                command += f" -{arg}"
            else:
                command += f" -{arg} {value}"

    # Print the command
    print(command)

    # Return the command
    return command


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Check if all mandatory arguments are provided
    check_mandatory_arguments(args)

    # Construct the bash command
    bash_command = construct_micapipe_command(args)

    # Safely parse the command into a format suitable for subprocess
    parsed_command = shlex.split(bash_command)

    # Execute the bash command
    subprocess.run(parsed_command)


if __name__ == "__main__":
    main()
