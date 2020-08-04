####################################################################################################

# Translated from matlab: https://github.com/MICA-MNI/micaopen/blob/master/MPC/scripts/03_surf2mpc.m
# Original script by Casey Paquola
# Translated to python by Jessica Royer, with tiny changes for MICs dataset structure (e.g. paths, parc_name)

# Description from original matlab script:
# "This script can be used as a wrapper to run the function build_mpc, and 
# thus construct microstructure profile covariance (MPC) matrices for a 
# group of individiuals. The following four variable need to be updated for
# your individual system/study. It will automatically save the intensity
# profiles and MPC matrices as a text files in the subject's BIDS folder, 
# alongside the preconstructed equivolumetric surfaces."

# INPUT
# dataDir 		BIDS derivatives directory
# sub 			subject id
# num_surf		surface solution number (default is 12)
# parc_name		name of parcellation in annotation file (default is vosdewael 200)

# EXAMPLE INPUTS FOR MICS
# dataDir = '/data_/mica3/BIDS_MIC/derivatives/' 
# sub = 'HC12'
# num_surf = 12 (Note that 14 surfaces were generated but we discard outermost and innermost surfaces)
# parc_name = 'vosdewael_200_native.csv'
####################################################################################################

# Import packages
import sys
import os
import numpy as np
import nibabel as nib
from build_mpc import build_mpc

# Define input arguments
dataDir = sys.argv[1] # = /host/fladgate/local_raid/jessica/
sub = sys.argv[2]
num_surf = sys.argv[3]
parc_name = sys.argv[4]

# Define default inpute if none given
if len(sys.argv) < 4:
    parc_name = 'vosdewael_200'

if len(sys.argv) < 3:
    num_surf = 12

# setting output directory
OPATH = "{dataDir}/sub-{sub}/ses-pre/proc_struct/surfaces/equivSurfs/{num_surf:d}surfs_out/".format(dataDir=dataDir, sub=sub, num_surf=num_surf+2)

if os.path.exists(OPATH):
    try:
        # Get data for specified hemisphere and surface number
        def get_hemisphere(surface_number, hemi):
            thisname_mgh = "{OPATH}{hemi}h.{surface_number:d}.mgh".format(OPATH=OPATH, hemi=hemi, surface_number=surface_number+1)
            img = nib.load(thisname_mgh)
            data = img.get_fdata()
            return data.reshape((1,-1))
        BBl = np.concatenate(
            [get_hemisphere(ii, 'l') for ii in range(int(num_surf))],
            axis=0
        )
        BBr = np.concatenate(
            [get_hemisphere(ii, 'r') for ii in range(int(num_surf))],
            axis=0
        )
        # Concatenate hemispheres and flip so pial surface is at the top
        BB = np.flipud(np.concatenate((BBl, BBr), axis = 1))
        # Load parcellation in native surface space
        pathToParc = "{dataDir}/{sub}/surfaces/fsspace_native_derivatives/{sub}_{parc_name}.csv".format(dataDir=dataDir, sub=sub, parc_name=parc_name)
        labeling = np.loadtxt(pathToParc, dtype=np.int)
        # Create MPC matrix (and nodal intensity profiles if parcellating)
        print("")
        print("--------------------------")
        print("Ello! Let's build the MPC!")
        print("--------------------------")
        print("")
        (MPC, I, problemNodes) = build_mpc(BB, labeling)
        # Check success of MPC and save output
        if np.isnan(np.sum(MPC)):
            print("")
            print("-------------------------------------")
            print("MPC building failed for subject {sub}".format(sub=sub))
            print("-------------------------------------")
            print("")
            sys.exit(0)
        else:
            np.savetxt("{output}/mpc_{parc_name}.txt".format(output=OPATH, parc_name=parc_name), MPC)
            np.savetxt("{output}/intensity_profiles_{parc_name}.txt".format(output=OPATH, parc_name=parc_name), I)
            print("")
            print("-------------------------------------")
            print("MPC building successful for subject {sub}".format(sub=sub))
            print("Check it out in {output}".format(output=OPATH))
            print("-------------------------------------")
            print("")
            sys.exit(1)
    except Exception as e:
        print("")
        print("---------------------------------------------------------------------")
        print("Something went wrong in loading or processing files for subject {sub}".format(sub=sub))
        print(e)
        print("---------------------------------------------------------------------")
        print("")
        sys.exit(-1)
