####################################################################################################

# Translated from matlab: https://github.com/MICA-MNI/micaopen/blob/master/MPC/scripts/build_mpc.m
# Original script by Casey Paquola
# Translated to python by Jessica Royer

# build_mpc   Construct a microstructure profile covariance matrix
  
# INPUT
# data          surfaces x vertices matrix containing intensity values
# parc          1 x vertices vector with unique integrers correpsonding to 
#               node assignment. Leave empty, [], for a vertex-wise MPC
#               Make sure that you provide different parcel # for each hemisphere
#
# OUTPUT
# MPC           microstructural profile covariance matrix
# I             microstructure intensity profiles (node-wise if parc is given)
# problemNodes  which nodes in intensity profiles are made up of NaNs

####################################################################################################

# Import packages
import sys
import numpy as np
import pandas as pd
from pingouin import partial_corr
import scipy
from scipy import special
import enigmatoolbox.datasets
from enigmatoolbox.utils.parcellation import map_to_labels, reduce_by_labels, relabel_consecutive

def build_mpc(data, parc=None):
    
    
    # If no parcellation is provided, MPC will be computed vertexwise
    if parc is None:
        downsample = 0
    else:
        downsample = 1


    # Parcellate input data according to parcellation scheme provided by user
    if downsample == 1:
        uparcel = np.unique(parc)
        I = np.zeros([data.shape[0], len(uparcel)])
        # Parcellate data by averaging profiles within nodes
        print("")
        print("----------------------------------------------------------")
        print("Parcellating microstructural profiles with provided labels")
        print("----------------------------------------------------------")
        print("")
        for ii in range(len(uparcel)):
            # Get vertices within parcel
            thisparcel = uparcel[ii]
            tmpData = data[:, parc == thisparcel]
            tmpData[:,np.mean(tmpData) == 0] = 0
            # Define function to find outliers: Return index of values above three scaled median absolute deviations of input data
            # https://www.mathworks.com/help/matlab/ref/isoutlier.html
            def find_outliers(data_vector):
                c = -1 / (np.sqrt(2) * scipy.special.erfcinv(3/2))
                scaled_MAD = c * np.median(np.abs(data_vector - np.median(data_vector)))
                is_outlier = np.greater(data_vector, (3 * scaled_MAD) + np.median(data_vector))
                idx_outlier = [i for i, x in enumerate(is_outlier) if x]
                return idx_outlier
            # Find if there are any outliers in vertex-wise average profile within given parcel
            idx = find_outliers(np.mean(tmpData, axis = 0))
            if len(idx) > 0: 
                tmpData[:,idx] = np.nan
            # Average profiles within parcels
            I[:,ii] = np.nanmean(tmpData, axis = 1)
        # Get matrix sizes
        szI = I.shape
        szZ = [len(uparcel), len(uparcel)]
    else:
        I = data
        szI = data.shape
        szZ = np.empty((data.shape[1], data.shape[1]))
    

    # Build MPC
    if np.isnan(np.sum(I)):
        # Find where are the NaNs
        is_nan = np.isnan(I[1,:])
        problemNodes = [i for i, x in enumerate(is_nan) if x]
        print("")
        print("---------------------------------------------------------------------------------")
        print("There seems to be an issue with the input data or parcellation. MPC will be NaNs!")
        print("---------------------------------------------------------------------------------")
        print("")
        # Fill matrices with NaN for return
        I = np.zeros(szI)
        I[I == 0] = np.nan
        MPC = np.zeros(szZ)
        MPC[MPC == 0] = np.nan
    else:
        problemNodes = 0
        # Have to convert profiles np array into pandas dataframe to use partial_corr function
        col_name = np.arange(szI[1]).astype(str)
        I_df = pd.DataFrame.from_records(I, columns = col_name)
        # Add a column at the end for the mean (also needed for partial_corr function)
        I_M = np.mean(I, axis = 1).reshape(-1,1)
        I_M_df = pd.DataFrame.from_records(I_M).rename(columns={0: "Mean"})
        I_df = pd.concat([I_df, I_M_df], axis = 1)
        # Loop through columns of I dataframe and correlate them, controlling for column "Mean"
        print("")
        print("--------------------------------------------------------")
        print("Running partial correlations, this could take a while...")
        print("--------------------------------------------------------")
        print("")
        R = np.zeros(szZ)
        for col in range(len(uparcel)):
            print("\rProgress: {}/{}".format(col + 1, len(uparcel)),end="")
            if col == len(uparcel)//2:
                print("\nSit tight we're halfway there!")
            for col_corr in range(len(uparcel)):
                if col != col_corr:
                    this_corr = partial_corr(data = I_df, x = col_name[col], y = col_name[col_corr], covar = 'Mean')
                    R[col,col_corr] = this_corr.r.pearson
                else:
                    R[col,col_corr] = 0        
        # Log transform
        MPC = 0.5 * np.log( np.divide(1 + R, 1 - R) )
        MPC[np.isnan(MPC)] = 0
        MPC[np.isinf(MPC)] = 0
        # Convert I back to numpy array
        I = I_df.to_numpy()
        I = I[:,:-1]

    # Output MPC, microstructural profiles, and problem nodes
    return (MPC, I, problemNodes)
    print("Goodbye!")
