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
# num_surf		surface solution number (default is 14)
# parc_name		name of parcellation in annotation file (default is vosdewael 200)

# EXAMPLE INPUTS FOR MICS
# dataDir = '/data_/mica3/BIDS_MIC/derivatives/' 
# sub = 'HC012'
# num_surf = 14
# parc_name = 'vosdewael-200_mics.annot'
####################################################################################################

# Import packages
import sys
import os
import numpy as np
import nibabel as nib
from build_mpc import build_mpc

# Define input arguments
dataDir = sys.argv[1] 
sub = sys.argv[2]
num_surf = sys.argv[3]
parc_name = sys.argv[4]
ses_num = sys.argv[5]

# Define default inpute if none given
if len(sys.argv) < 4:
    parc_name = 'vosdewael-200_mics.annot'

if len(sys.argv) < 3:
    num_surf = 14

# setting output directory
OPATH = "{dataDir}/sub-{sub}/ses-{ses}/proc_struct/surfaces/micro_profiles/".format(dataDir=dataDir, sub=sub, ses=ses_num)

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
        pathToParc = "{dataDir}/sub-{sub}/ses-{ses}/proc_struct/surfaces/{sub}/label/".format(dataDir=dataDir, sub=sub, ses=ses_num)
        # Load annot files
        fname_lh = 'lh.' + parc_name
        ipth_lh = os.path.join(pathToParc, fname_lh)
        [labels_lh, ctab_lh, names_lh] = nib.freesurfer.io.read_annot(ipth_lh, orig_ids=True)
        fname_rh = 'rh.' + parc_name
        ipth_rh = os.path.join(pathToParc, fname_rh)
        [labels_rh, ctab_rh, names_rh] = nib.freesurfer.io.read_annot(ipth_rh, orig_ids=True)
        # Join hemispheres
        parcLength = len(labels_lh)+len(labels_rh)
        parc = np.zeros((parcLength))
        for x in range(len(labels_lh)):
            parc[x] = np.where(ctab_lh[:,4] == labels_lh[x])[0][0]
        for x in range(len(labels_rh)):
            parc[x + len(labels_lh)] = np.where(ctab_rh[:,4] == labels_rh[x])[0][0] + len(ctab_lh)
        
        # Create MPC matrix (and nodal intensity profiles if parcellating)
        print("")
        print("--------------------------")
        print("Ello! Let's build the MPC!")
        print("--------------------------")
        print("")
        (MPC, I, problemNodes) = build_mpc(BB, parc)
        
        # Check success of MPC and save output
        if np.isnan(np.sum(MPC)):
            print("")
            print("-------------------------------------")
            print("MPC building failed for subject {sub}".format(sub=sub))
            print("-------------------------------------")
            print("")
            sys.exit(1)
        else:
            parc_str = parc_name.replace('.annot', "")
            np.savetxt("{output}/mpc_{parc_str}.txt".format(output=OPATH, parc_str=parc_str), MPC)
            np.savetxt("{output}/intensity_profiles_{parc_str}.txt".format(output=OPATH, parc_str=parc_str), I)
            print("")
            print("-------------------------------------")
            print("MPC building successful for subject {sub}".format(sub=sub))
            print("Check it out in {output}".format(output=OPATH))
            print("-------------------------------------")
            print("")
            sys.exit(0)
    except Exception as e:
        print("")
        print("---------------------------------------------------------------------")
        print("Something went wrong in loading or processing files for subject {sub}".format(sub=sub))
        print(e)
        print("---------------------------------------------------------------------")
        print("")
        sys.exit(-1)
