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
import numpy as np
import scipy.special
import scipy.stats

def build_mpc(data, parc=None, idxExclude=None):
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
        for (ii, _) in enumerate(uparcel):

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

        # Calculate mean across columns, excluding mask and any excluded labels input
        I_mask = I
        for i in idxExclude:
            I_mask[:, i] = np.nan
        I_M = np.nanmean(I_mask, axis = 1)

        # Get residuals of all columns (controlling for mean)
        I_resid = np.zeros(I.shape)
        for c in range(I.shape[1]):
            y = I[:,c]
            x = I_M
            slope, intercept, _, _, _ = scipy.stats.linregress(x,y)
            y_pred = intercept + slope*x
            I_resid[:,c] = y - y_pred

        R = np.corrcoef(I_resid, rowvar=False)

        # Log transform
        MPC = 0.5 * np.log( np.divide(1 + R, 1 - R) )
        MPC[np.isnan(MPC)] = 0
        MPC[np.isinf(MPC)] = 0

        # CLEANUP: correct diagonal and round values to reduce file size
        # Replace all values in diagonal by zeros to account for floating point error
        for i in range(0,MPC.shape[0]):
                MPC[i,i] = 0
        # Replace lower triangle by zeros
        MPC = np.triu(MPC)

    # Output MPC, microstructural profiles, and problem nodes
    return (MPC, I, problemNodes)
