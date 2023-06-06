#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun May 14 21:17:57 2023

This script will create a pdf with the QC at group level of the outputs of the micapipe, generating:
    a) A table of the modules by cases processed
    b) A graph of the each module progress (in percertange)
    c) The processing times displayed in graph and table by subject per module
    
    STRUCTURAL
    d) The mean group cortical thickness (vertex-wise)
    e) cortical thickness similarity by subject
    e.1) gradient of the mean GD (scaheffer-400)
    
    DWI (schaefer-400)
    f) IF dwi has been process, mean ADC map
    g) ADC similarity by subject
    g.1) gradient of the mean SC
    
    FUNC (one per module. e.g. rest, task etc. schaefer-400)
    h) Mean functional connectome
    i) similarity matrix of func data between subjects
    i.1) gradient of the mean FC
    
    MAPS (schaefer-400)
    j) mean qMRI map
    k) similarity matrix of the qMRI data
    k.1) gradient of the mean MPC
    
    flair/T2 (vertex-wise)
    l) mean flair image
    m) similarity matrix of all flair/T2 data
    
    Structural connectome (scheffer-400)
    n) mean strength
    o) similarity matrix
    p) gradient (scheffer-400)

    Parameters
    ----------

    out       : Output directory for the processed files <derivatives>

    bids      : Path to BIDS Directory

    tmpDir    : OPTIONAL specifcation of temporary directory location <path> (Default is /tmp)

    version   : OPTIONAL print software version


    USAGE
    --------
    >>> QC_group.py -out <outputDirectory> -bids <BIDS-directory>

  McGill University, MNI, MICA-lab, Created 13 September 2022
  https://github.com/MICA-MNI/micapipe
  http://mica-mni.github.io/

  author: alexander-ngo    
    
@author: rcruces
"""
from xhtml2pdf import pisa
import pandas as pd
import numpy as np
import os
import glob
import json
import nibabel as nb
import argparse
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting import plot_hemispheres
from brainspace.utils.parcellation import map_to_labels
import seaborn as sns
import math
from matplotlib.cm import Spectral
import matplotlib.pyplot as plt
import cmocean
cmaps = cmocean.cm.cmap_d

# Arguments
parser = argparse.ArgumentParser()

# Required
parser.add_argument('-out',
                    dest='out',
                    type=str,
                    help='Output directory for the processed files <derivatives>',
                    required=True
                    )

parser.add_argument('-bids',
                    dest='bids',
                    type=str,
                    help='Path to BIDS Directory',
                    )

parser.add_argument('-tmpDir',
                    dest='tmpDir',
                    type=str,
                    default='/tmp',
                    help='Path to temporary directory location'
                    )

parser.add_argument('-version',
                    action='store_true',
                    help='Print software version used to run the dataset'
                    )

args = parser.parse_args()

# Arguments
out = os.path.realpath(args.out)
version = args.version

# Set the paths
derivatives = out.split('/micapipe_v0.2.0')[0]
out=derivatives+'/micapipe_v0.2.0'

# Check pipeline status
def check_json_exist(jsonPath=''):
    json_exist = os.path.isfile(jsonPath)
    return json_exist

data_desc = out+'/dataset_description.json'
if check_json_exist(data_desc):
    with open( data_desc ) as f:
        data_json = json.load(f)
    generated_by = data_json["GeneratedBy"][0]
else:
    print('ERROR... dataset_description.json does not exit, or this dataset has not been processed with micapipe')
    exit()

# Print version
if version == True:
    print(generated_by["Name"] + ' ' + generated_by["Version"])
    exit()

# ----------------------------------------------------------- #
#bids = os.path.realpath(args.bids)
# Check if tmpDir exists, otherwise set it to '/tmp' as default
tmpDir = args.tmpDir
if not os.path.exists(tmpDir):
    tmpDir = '/tmp'

tmpDir = tmpDir+'/micapipe_QC-group_'+str(np.random.randint(low=1000, high=10000))
os.mkdir(tmpDir)

# ----------------------------------------------------------- #
# Path to MICAPIPE
micapipe=os.popen("echo $MICAPIPE").read()[:-1]

## ------------------------------------------------------------------------- ##
##                                                                           ##
##                              Helper functions                             ##
##                                                                           ##
## ------------------------------------------------------------------------- ##
# Set dataset PNI as working directory
os.chdir(out)

# Load native mid surface
inf_lh = read_surface(micapipe + '/surfaces/fsLR-32k.L.inflated.surf.gii', itype='gii')
inf_rh = read_surface(micapipe + '/surfaces/fsLR-32k.R.inflated.surf.gii', itype='gii')
i5_lh = read_surface(micapipe + '/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
i5_rh = read_surface(micapipe + '/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')

def plot_connectome(mtx, Title='matrix plot', xlab='X', ylab='Y', col='rocket', vmin=None, vmax=None,
                   xticklabels='auto', yticklabels='auto',xrot=90, yrot=0, save_path=None):
    
    '''
    This optional function, only plots a connectome as a heatmap
    Parameters
    ----------
    mtx : np.array
    Returns
    -------
    f : plot
    '''
    f, ax = plt.subplots(figsize=(15,10))
    g = sns.heatmap(mtx, ax=ax, cmap=col, vmin=vmin, vmax=vmax, xticklabels=xticklabels, yticklabels=yticklabels)
    g.set_xlabel(xlab)
    g.set_ylabel(ylab)
    g.set_title(Title)
    # Rotate the x-axis labels
    # rotate tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=xrot, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=yrot, ha='right')
    
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

def load_mpc(File, Ndim):
    """Loads and process a MPC"""
    
    # load the matrix
    mtx_mpc = nb.load(File).darrays[0].data
    
    # Mirror the matrix
    MPC = np.triu(mtx_mpc,1)+mtx_mpc.T
    
    # Remove the medial wall
    MPC = np.delete(np.delete(MPC, 0, axis=0), 0, axis=1)
    MPC = np.delete(np.delete(MPC, Ndim, axis=0), Ndim, axis=1)
    
    return MPC

def load_gd(File, Ndim):
    """Loads and process a GD"""
    
    # load the matrix
    mtx_gd = nb.load(File).darrays[0].data
    
    # Remove the Mediall Wall
    mtx_gd = np.delete(np.delete(mtx_gd, 0, axis=0), 0, axis=1)
    GD = np.delete(np.delete(mtx_gd, Ndim, axis=0), Ndim, axis=1)
    
    return GD

def load_fc(File, Ndim, parc=''):
    """Loads and process a functional connectome"""
    
    # load the matrix
    mtx_fs = nb.load(File).darrays[0].data
    
    # slice the matrix remove subcortical nodes and cerebellum
    FC = mtx_fs[49:, 49:]
    
    # Fisher transform
    FCz = np.arctanh(FC)

    # replace inf with 0
    FCz[~np.isfinite(FCz)] = 0
    
    # Mirror the matrix
    FCz = np.triu(FCz,1)+FCz.T
    return FCz

def load_sc(File, Ndim):
    """Loads and process a structura connectome"""
    
    # load the matrix
    mtx_sc = nb.load(File).darrays[0].data
    
    # Mirror the matrix
    mtx_sc = np.log(np.triu(mtx_sc,1)+mtx_sc.T)
    mtx_sc[np.isneginf(mtx_sc)] = 0
    
    # slice the matrix remove subcortical nodes and cerebellum
    SC = mtx_sc[49:, 49:]
    SC = np.delete(np.delete(SC, Ndim, axis=0), Ndim, axis=1)
    
    # replace 0 values with almost 0
    SC[SC==0] = np.finfo(float).eps
    
    return SC

# Load files
def load_connectomes(files, Ndim, func):
    # Load all the matrices
    M=np.empty([Ndim*2, Ndim*2, len(files)], dtype=float)
    for i, f in enumerate(files):
        M[:,:,i] = func(f, Ndim)
        
    return M

# Load annotation file in fsaverage5
atlas='schaefer-400'
annot_lh_fs5= nb.freesurfer.read_annot(micapipe + '/parcellations/lh.' + atlas + '_mics.annot')
Ndim = max(np.unique(annot_lh_fs5[0]))

# Read label for conte69
labels_c69 = np.loadtxt(open(micapipe + '/parcellations/'+atlas+'_conte69.csv'), dtype=int)

# mask of the medial wall
mask_c69 = labels_c69 != 0

## ---------------------------- MICAPIPE header ---------------------------- ##


# --------------------------------------------------------------------
# json csv and table
# --------------------------------------------------------------------
# Get the list of json files and sessions
dir_list = glob.glob(out+'/sub*/*')
dir_ses = [entry.split('/')[-1] for entry in dir_list]
dir_ses = list(set(dir_ses))

if any('ses-' in string for string in dir_ses):
    dir_str='sub*/ses*'
else:
    dir_str = 'sub*'
jsons = sorted(glob.glob(out+'/'+dir_str+'/QC/*json'))   

# Sort the JSON files by creation time (newest at the bottom)
jsons = sorted(jsons, key=lambda x: os.path.getctime(x))
                   
# Define the keys to extract from each JSON file
keys_to_extract = ['Subject', 'Session', 'Module', 'Status', 'Progress', 'User', 'Workstation', 'Date', 'Processing.time', 'Processing', 'micapipeVersion']

# Initialize an empty list to hold the individual dataframes
dataframes = []

# Loop through each JSON file in the specified directory
for file_name in jsons:

    # Read in the JSON data from the file
    with open(file_name, 'r') as f:
        json_data = json.load(f)

    # Extract the desired keys and create a dataframe
    data = {k: [] for k in keys_to_extract}
    
    for k in keys_to_extract:
        data[k].append(json_data[k])

    df = pd.DataFrame(data)

    # Append the dataframe to the list
    dataframes.append(df)

# Concatenate all the individual dataframes into a single dataframe
result = pd.concat(dataframes, ignore_index=True)

# Rename Version column
result.rename(columns={'micapipeVersion': 'Version'})

# Save the resulting dataframe as a CSV file
result.to_csv(out+'/micapipe_processed_sub.csv', index=False)

# --------------------------------------------------------
# Create an empty dictionary to hold the data
data = {}

# Loop through each row in the concatenated dataframe
for i in range(len(result)):

    # Extract thesubject and session and module and status from the row
    subject = result['Subject'][i]
    session = result['Session'][i]
    module = result['Module'][i]
    status = result['Status'][i]

    # Convert the status to a binary value
    if status == 'INCOMPLETE':
        value = 0
    elif status == 'COMPLETED':
        value = 1
    else:
        value = np.nan

    # Add the value to the data dictionary
    key = f"{subject}-{session}"
    if key not in data:
        data[key] = {}
    data[key][f"{module}"] = value

# Create a dataframe from the data dictionary
df = pd.DataFrame(data).T

# Sort the columns by module and status
df = df.reindex(sorted(df.columns), axis=1)

# Sort the rows by subject and session
df = df.sort_index()


# Check if the number of rows in df is greater than 50
if len(df) > 50:
    # Calculate the number of subplots needed
    num_subplots = math.ceil(len(df) / 50)
    # Loop through the subplots and create each one
    for i in range(num_subplots):
        start_idx = i * 50
        end_idx = (i+1) * 50
        subset_df = df.iloc[start_idx:end_idx]
        # Determine if there are both 0 and 1 values in the data
        has_zeros = any(0 in row.values for _, row in subset_df.iterrows())
        has_ones = any(1 in row.values for _, row in subset_df.iterrows())

        # Choose the appropriate color map
        if has_zeros and has_ones:
            cmap = 'vlag_r'  # Reverse color map if both 0 and 1 values exist
        else:
            cmap = 'vlag'  # Use original color map if only 0 or 1 values exist
        # Set the size of the plot according to the number of rows in the subset DataFrame
        fig, ax = plt.subplots(figsize=(8, len(subset_df)*0.3))
        sns.heatmap(subset_df, cmap=cmap, annot=True, cbar=False, fmt='.0f', linewidths=1, linecolor='white', annot_kws={'fontsize':12, 'color':'white'})
        # Set the x and y axis labels to white
        plt.xlabel('Module-Status', fontsize=14, color='white')
        plt.ylabel('Subject-Session', color='white')
        # Set the plot title and axis labels
        plt.title(f'Module Status by Subject-Session ({start_idx+1}-{end_idx})')
        plt.xlabel('Module-Status')
        plt.ylabel('Subject-Session')
        # Rotate the x-axis labels for better visibility
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        # Save the plot as a PNG file
        plt.savefig(f'{tmpDir}/micapipe_qc_module_status_{start_idx+1}_{end_idx}.png', dpi=300, bbox_inches='tight')
        # Close the plot to release resources
        plt.close()
else:
    # Determine if there are both 0 and 1 values in the data
    has_zeros = any(0 in row.values for _, row in df.iterrows())
    has_ones = any(1 in row.values for _, row in df.iterrows())

    # Choose the appropriate color map
    if has_zeros and has_ones:
        cmap = 'vlag_r'  # Reverse color map if both 0 and 1 values exist
    else:
        cmap = 'vlag'  # Use original color map if only 0 or 1 values exist
    # Set the size of the plot according to the number of rows in df
    fig, ax = plt.subplots(figsize=(8, len(df)*0.3))
    sns.heatmap(df, cmap=cmap, annot=True, cbar=False, fmt='.0f', linewidths=1, linecolor='white', annot_kws={'fontsize':12, 'color':'white'})
    # Set the x and y axis labels to white
    plt.xlabel('Module-Status', fontsize=14, color='white')
    plt.ylabel('Subject-Session', color='white')
    # Set the plot title and axis labels
    plt.title('Module Status by Subject-Session')
    plt.xlabel('Module-Status')
    # Save the plot as a PNG file
    plt.savefig(f'{tmpDir}/micapipe_qc_module_status.png', dpi=300, bbox_inches='tight')
    # Close the plot to release resources
    plt.close()

# --------------------------------------------------------------------
# Progress barplot
# --------------------------------------------------------------------
# Get progress percentage
module_progress = np.sum(df, axis=0)/len(df)

# sort the series by values
module_progress_sorted = module_progress.sort_values(ascending=False)

# create the horizontal barplot
fig, ax = plt.subplots(figsize=(8, 6))
module_progress_sorted.plot.barh(ax=ax, color=Spectral(module_progress_sorted))

# set the x and y axis labels and title
ax.set_xlabel('Progress percentage')
ax.set_title('Module Progress')

# Save the figure as a PNG file
plt.savefig(f'{tmpDir}/micapipe_qc_module_progress_plot.png', dpi=300, bbox_inches='tight')

# Close the plot to release resources
plt.close()

# --------------------------------------------------------------------
# Processing time plot
# --------------------------------------------------------------------
# Convert the processing.time to numeric
result['Processing.time'] = pd.to_numeric(result['Processing.time'])

# Calculate the mean and standard deviation by module and sort by mean processing time
stats_df = result.groupby('Module')['Processing.time'].agg(['mean', 'std'])
stats_df = stats_df.sort_values('mean')
stats_df['mean_std'] = stats_df.apply(lambda row: '{:.2f} +/- {:.2f}'.format(row['mean'], row['std']), axis=1)
stats_df = stats_df.reset_index()

# Create the box plot with custom box and whisker colors and sorted order
sns.boxplot(x='Module', y='Processing.time', data=result, order=stats_df['Module'], color='thistle')

# Add the individual data points as dots
sns.stripplot(x='Module', y='Processing.time', data=result, order=stats_df['Module'], color='purple')

# Set the title of the plot
plt.title('Processing Time')

# Rotate the x-axis labels
plt.xticks(rotation=90)
# labels
plt.xlabel('')
plt.ylabel('minutes')
plt.title('Processing time')

# Save the plot as a PNG file
plt.savefig(f'{tmpDir}/micapipe_qc_time.png', dpi=300, bbox_inches='tight')
# Close the plot to release resources
plt.close()

# --------------------------------------------------------------------
# Processing time table
# --------------------------------------------------------------------
# Calculate the mean and standard deviation by module and sort by mean processing time
stats_df = result.groupby('Module')['Processing.time'].agg(['mean', 'std'])
stats_df = stats_df.sort_values('mean')
stats_df['mean_std'] = stats_df.apply(lambda row: '{:.1f} ± {:.1f}'.format(row['mean'], row['std']), axis=1)

# Generate the HTML table with pandas Styler
html_table = stats_df[['mean_std']].style.set_table_styles([
    {'selector': 'th', 'props': [('text-align', 'left'), ('background-color', '#f2f2f2'), ('padding', '1px')]},
    {'selector': 'td', 'props': [('border', '1px solid #ddd'), ('padding', '1px'), ('white-space', 'nowrap')]},
])
# Render the styled table
styled_table = html_table.to_html()



# --------------------------------------------------------------------
# Write pdf
# --------------------------------------------------------------------
def report_header_template(dataset_name='', micapipe=''):
    # Header
    report_header = (
        # Micapipe banner
        '<img id=\"top\" src=\"{micapipe}/micapipe_long.png\" alt=\"micapipe\">'

        # Dataset name
        '<h1 style="color:#343434;font-family:Helvetica, sans-serif !important;text-align:center;margn-bottom:0">'
        '{dataset_name} <h1>'

        # Subject's ID | Session
        '<h3 style="color:#343434;font-family:Helvetica, sans-serif;text-align:center;margin-bottom:25px">'
        '<b>Group level QC</b>'
    )
    return report_header.format(dataset_name=dataset_name, micapipe=micapipe)

def qc_header():
    dataset_description = os.path.realpath("%s/dataset_description.json"%(out))
    with open( dataset_description ) as f:
        dataset_description = json.load(f)
    dataset_name = dataset_description["Name"]
    _static_block = report_header_template(dataset_name=dataset_name, micapipe=micapipe)

    return _static_block

# Header template
def report_module_header_template(module=''):
    # Module header:
    report_module_header = (
        '<p style="border:2px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica, '
        'sans-serif;font-size:14px;text-align:center;">'
        '<b>{module}</b> <p>'
    )
    return report_module_header.format(module=module)

# OUTPUT TEMPLATE
def report_module_output_figure(outName='', figPath=''):
    # Module Outputs
    report_module_output = (
        '<p style="font-family:Helvetica, sans-serif;font-size:10px;text-align:Left;margin-bottom:0px">'
        '<b> {outName} </b> </p>'

        '<center> <img style="width:500px%;margin-top:0px" src="{figPath}"> </center>'
    )
    return report_module_output.format(outName=outName, figPath=figPath)

def report_titleh2(Title='', Size='12'):
    # Module Outputs
    report_module_output = (
        '<h2 style="color:#343434;font-family:Helvetica, sans-serif !important;text-align:Left;margin-bottom:0">'
        '{Title} </h2>'
    )
    return report_module_output.format(Title=Title)

# --------------------------------------------------------------------
# Surface and similarity matrix
# --------------------------------------------------------------------
def report_surface_similarity(out, lh_str, rh_str, out_png, cmap, quantile=(0.01, 0.95)):
    # change working directory
    os.chdir(out)
    # Load and plot the mean maps on the surface
    # Load all the cortical thickness
    lh_files=sorted(glob.glob(lh_str))
    rh_files=sorted(glob.glob(rh_str))
    
    if len(lh_files) > 0:
        # get the BIDS ids
        bids_ids = [elem.split('/')[0] + '_' + elem.split('/')[1] for elem in lh_files]
        
        # Load all the thickness data
        Nth=np.concatenate((nb.load(lh_files[0]).darrays[0].data, nb.load(rh_files[0]).darrays[0].data), axis=0).shape[0]
        
        surf_map=np.empty([len(lh_files), Nth], dtype=float)
        for i, f in enumerate(lh_files):
            #print(f)
            surf_map[i,:] = np.hstack(np.concatenate((nb.load(lh_files[i]).darrays[0].data, nb.load(rh_files[i]).darrays[0].data), axis=0))
        # Mean matrix across the x axis (vertices)
        map_mean = np.mean(surf_map, axis=0)
        
        # Plot the mean thickness 10mm on conte69 surface
        Range=(np.quantile(map_mean, quantile[0]), np.quantile(map_mean, quantile[1]))
        plot_hemispheres(i5_lh, i5_rh, array_name=map_mean, cmap=cmap, nan_color=(0, 0, 0, 1),
                              zoom=1.3, size=(900, 750), embed_nb=True,
                              color_bar='right', layout_style='grid', color_range=Range,
                              label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                              screenshot=True, filename=f'{tmpDir}/micapipe_qc_{out_png}.png')
        # Calculate the similarity matrix
        ## correlation matrix
        corr = np.corrcoef(surf_map)
        plot_connectome(corr, f'Subject similarity - {out_png}', xlab=None, ylab=None, col='Spectral_r', vmin=0.1, vmax=1,
                       yticklabels=bids_ids, xticklabels=bids_ids, save_path=f'{tmpDir}/micapipe_qc_{out_png}_matrix.png')
        plt.close()
        # plot a barplot if the subjects mean similarity is below a threshold 0.2???
        # Exclude diagonal values
        np.fill_diagonal(corr, np.nan)
        
        # Calculate column means
        colmean = np.nanmean(corr, axis=0)
        
        # Get the position indices where colmean < 0.2
        indices = np.where(colmean < 0.25)[0]
        
        # Write html code
        _static_block = '<div style="page-break-after: always;"></div>'
        _static_block +=  report_module_header_template(module=f'{out_png} | vertex-wise')
        _static_block += report_titleh2(f'{out_png} | vertex-wise Group Mean')
        _static_block += report_module_output_figure('', f'{tmpDir}/micapipe_qc_{out_png}.png')
        _static_block += report_titleh2(f'{out_png} | between Subjects Similarity')
        _static_block += report_module_output_figure('', f'{tmpDir}/micapipe_qc_{out_png}_matrix.png')
        
        # Print the position indices
        if len(indices) > 0:
            # Create a dictionary with the data for the table
            table_data = {
                "Subject": [bids_ids[i] for i in indices],
                "Mean similarity Value": colmean[indices]
            }
            
            # Create a DataFrame from the table data
            table_df = pd.DataFrame(table_data)
            
            # Apply the table styles
            table_style = [
                {'selector': 'th', 'props': [('text-align', 'left'), ('background-color', '#f2f2f2'), ('padding', '1px')]},
                {'selector': 'td', 'props': [('border', '1px solid #ddd'), ('padding', '1px'), ('white-space', 'nowrap')]}
            ]
            styled_table = table_df.style.set_table_styles(table_style)
            
            # Display the styled table
            _static_block += report_titleh2(f'{out_png} | subjects with low within group similarity')
            _static_block += styled_table.to_html()
    else:
        _static_block = ''
        
    return(_static_block)

def report_roi_similarity(out, file_str, out_png, cmap, load_cnn):
    # change working directory
    os.chdir(out)
    # Loads and plots all the connectomes
    files=sorted(glob.glob(file_str))
    
    if len(files) > 0:
        # get the BIDS ids
        bids_ids = [elem.split('/')[0] + '_' + elem.split('/')[1] for elem in files]
        
        # Load all the thickness data
        mtxs=load_connectomes(files, Ndim, load_cnn)

        # Mean matrix across the z axis (subjects)
        mtxs_mean = np.mean(mtxs, axis=2)
        # Mean acros columns
        mtxs_colM = np.mean(mtxs_mean, axis=1)
        # map to labels
        mtxs_surf = map_to_labels(mtxs_colM, labels_c69, fill=np.nan, mask=mask_c69)
        
        # Plot the mean thickness 10mm on conte69 surface
        plot_hemispheres(inf_lh, inf_rh, array_name=mtxs_surf, cmap=cmap, nan_color=(0, 0, 0, 1),
                              zoom=1.3, size=(900, 750), embed_nb=True,
                              color_bar='right', layout_style='grid', color_range='sym',
                              label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                              screenshot=True, filename=f'{tmpDir}/micapipe_qc_{out_png}-roi.png')
        # Calculate the similarity matrix
        ## correlation matrix
        corr = np.corrcoef(np.mean(mtxs, axis=1).T)
        plot_connectome(corr, f'Subject similarity - {out_png}', xlab=None, ylab=None, col='Spectral_r', vmin=0.1, vmax=1,
                       yticklabels=bids_ids, xticklabels=bids_ids, save_path=f'{tmpDir}/micapipe_qc_{out_png}_matrix-roi.png')
        plt.close()
        # plot a barplot if the subjects mean similarity is below a threshold 0.2???
        # Exclude diagonal values
        np.fill_diagonal(corr, np.nan)
        
        # Calculate column means
        colmean = np.nanmean(corr, axis=0)
        
        # Get the position indices where colmean < 0.2
        indices = np.where(colmean < 0.25)[0]
        
        # Write html code
        _static_block = '<div style="page-break-after: always;"></div>'
        _static_block +=  report_module_header_template(module=f'{out_png} | ROI Schaeffer-400')
        _static_block += report_titleh2(f'{out_png} | ROI Schaeffer-400 Group Mean')
        _static_block += report_module_output_figure('', f'{tmpDir}/micapipe_qc_{out_png}-roi.png')
        _static_block += report_titleh2(f'{out_png} | between Subjects Similarity')
        _static_block += report_module_output_figure('', f'{tmpDir}/micapipe_qc_{out_png}_matrix-roi.png')
        
        # Print the position indices
        if len(indices) > 0:
            # Create a dictionary with the data for the table
            table_data = {
                "Subject": [bids_ids[i] for i in indices],
                "Mean similarity Value": colmean[indices]
            }
            
            # Create a DataFrame from the table data
            table_df = pd.DataFrame(table_data)
            
            # Apply the table styles
            table_style = [
                {'selector': 'th', 'props': [('text-align', 'left'), ('background-color', '#f2f2f2'), ('padding', '1px')]},
                {'selector': 'td', 'props': [('border', '1px solid #ddd'), ('padding', '1px'), ('white-space', 'nowrap')]}
            ]
            styled_table = table_df.style.set_table_styles(table_style)
            
            # Display the styled table
            _static_block += report_titleh2(f'{out_png} | subjects with low within group similarity')
            _static_block += styled_table.to_html()
    else:
        _static_block = ''
        
    return(_static_block)

def get_acqs(Str):
    acqs = sorted(glob.glob(out+'/'+dir_str+f'/{Str}/*'))
    acqs_uni = [entry.split('/')[-1] for entry in acqs]
    acqs_uni = list(set(acqs_uni))
    if Str == 'dwi' and not any('acq-' in entry for entry in acqs_uni):
        acqs_uni = ['-']
    
    return(acqs_uni)

def get_tracts(dwi_str):
    acqs = sorted(glob.glob(out+f'/{dwi_str}/connectomes/*schaefer-400*full-connectome*'))
    acqs_uni = [entry.split('iFOD2-')[-1] for entry in acqs]
    tracts = [entry.replace('-SIFT2_full-connectome.txt','').split('_')[0].split('-')[0] for entry in acqs_uni]
    
    tracts = list(set(tracts))
    
    return tracts

def report_micapipe():
    # Read the JSON data from file
    with open(f'{out}/dataset_description.json', 'r') as json_file:
        json_data = json.load(json_file)
    
    # Extract the relevant information from the JSON
    info = {
        "Info": ["Name", "Version", "Reference", "DOI", "URL", "GitHub", "Tag", "RunBy", "Workstation", "LastRun", "Processing"],
        "Description": [json_data["Name"],
                        json_data["GeneratedBy"][0]["Version"],
                        json_data["GeneratedBy"][0]["Reference"],
                        json_data["GeneratedBy"][0]["DOI"],
                        json_data["GeneratedBy"][0]["URL"],
                        json_data["GeneratedBy"][0]["GitHub"],
                        json_data["GeneratedBy"][0]["Container"]["Tag"],
                        json_data["GeneratedBy"][0]["RunBy"],
                        json_data["GeneratedBy"][0]["Workstation"],
                        json_data["GeneratedBy"][0]["LastRun"],
                        json_data["GeneratedBy"][0]["Processing"]]
    }
    
    # Create a DataFrame from the extracted information
    table_df = pd.DataFrame(info)
    
    # Remove the row index
    table_df.index = [""] * len(table_df)
           
    # Apply the table styles
    table_style = [
               {'selector': 'th', 'props': [('text-align', 'left'), ('background-color', '#f2f2f2'), ('padding', '1px')]},
               {'selector': 'td', 'props': [('border', '1px solid #ddd'), ('padding', '1px'), ('white-space', 'nowrap')]}
           ]
    styled_table = table_df.style.set_table_styles(table_style)
           
    # Display the styled table
    _static_block = report_titleh2('Processing details')
    _static_block += styled_table.to_html()
           
    return(_static_block)

def qc_group():

    # QC header
    _static_block = qc_header()
    # -------------------------------------------------------
    # Progress
    _static_block +=  report_module_header_template(module='Processing progress')
    
    print( 'Creating... plot tables')
    tables_png = sorted(glob.glob(tmpDir+'/micapipe_qc_module_status*.png'))
    for png in tables_png:
        _static_block += report_module_output_figure('', png)
    print( 'Creating... Module progress')
    _static_block += report_module_output_figure('', f'{tmpDir}/micapipe_qc_module_progress_plot.png')
    print( 'Creating... Processing times')
    _static_block += report_module_output_figure('', f'{tmpDir}/micapipe_qc_time.png')
    
    _static_block += styled_table.replace('mean_std', 'Time in minutes (mean | SD)')
    
    # -------------------------------------------------------
    print( 'Creating... Mean vertex-wise maps')
    # Vertex-wise Thickness
    _static_block += report_surface_similarity(out, f'{dir_str}/maps/*_hemi-L_surf-fsLR-5k_label-thickness.func.gii',
                              f'{dir_str}/maps/*_hemi-R_surf-fsLR-5k_label-thickness.func.gii', 'thickness', 'rocket', quantile=(0.075, 0.995))

    # Vertex-wise DWI derived map (dynamic)
    _static_block += report_surface_similarity(out, f'{dir_str}/maps/*_hemi-L_surf-fsLR-5k_label-white_ADC.func.gii',
                              f'{dir_str}/maps/*_hemi-R_surf-fsLR-5k_label-white_ADC.func.gii', 'ADC', 'cmo.ice')
    
    # Vertex-wise func (dynamic)
    for acq in get_acqs('func'):
        _static_block += report_surface_similarity(out, f'{dir_str}/func/{acq}/surf/*_hemi-L_surf-fsLR-5k.func.gii',
                                  f'{dir_str}/func/{acq}/surf/*_hemi-R_surf-fsLR-5k.func.gii', acq.replace('desc-',''), 'cmo.dense_r', quantile=(0.05, 0.9))
    
    # Vertex-wise MPC (dynamic)
    for acq in get_acqs('mpc'):
        acq_qmri=acq.replace('acq-','')
        _static_block += report_surface_similarity(out, f'{dir_str}/maps/*_hemi-L_surf-fsLR-5k_label-midthickness_{acq_qmri}.func.gii',
                                  f'{dir_str}/maps/*_hemi-R_surf-fsLR-5k_label-midthickness_{acq_qmri}.func.gii', 'MPC-'+acq_qmri, 'crest', quantile=(0.01, 0.99))
    
    # Vertex-wise flair
    _static_block += report_surface_similarity(out, f'{dir_str}/maps/*_hemi-L_surf-fsLR-5k_label-midthickness_flair.func.gii',
                              f'{dir_str}/maps/*_hemi-R_surf-fsLR-5k_label-midthickness_flair.func.gii', 'flair', 'afmhot', quantile=(0.15, 0.999))

    # -------------------------------------------------------
    # ROI based
    print( 'Creating... ROI based maps')
    # ROI GDs
    _static_block += report_roi_similarity(out, f'{dir_str}/dist/*_atlas-schaefer-400_GD.shape.gii', 'GD', 'vlag', load_gd)
    
    # ROI SC (dynamic)
    for acq in get_acqs('dwi'):
        acq_dir=f'{dir_str}/dwi/{acq}/'.replace('/-/','/')
        dwi_tracts=get_tracts(acq_dir)
        for tracts in dwi_tracts:
            dwi_name=f'SC_{acq}_{tracts}'.replace('_-','')
            print( f'   dwi id: {dwi_name}')
            cnn_files=f'{dir_str}/dwi/{acq}/connectomes/*atlas-schaefer-400_desc-iFOD2-{tracts}-SIFT2_full-connectome.shape.gii'.replace('/-/','/')
            _static_block += report_roi_similarity(out, cnn_files, dwi_name, 'flare_r', load_sc)
    
    # ROI func (dynamic)
    for acq in get_acqs('func'):
        func_name='FC'+acq.replace('desc-','')
        print( f'   func id: {func_name}')
        _static_block += report_roi_similarity(out, f'{dir_str}/func/'+acq+'/surf/*_atlas-schaefer-400_desc-FC.shape.gii', func_name, 'cmo.dense_r', load_fc)
    
    # ROI MPC (dynamic)
    for acq in get_acqs('mpc'):
        acq_qmri=acq.replace('acq-','')
        _static_block += report_roi_similarity(out, f'{dir_str}/mpc/'+acq+'/*_atlas-schaefer-400_desc-MPC.shape.gii', 'MPC-'+acq_qmri, 'crest', load_mpc)

    _static_block += '<div style="page-break-after: always;"></div>'
    _static_block += qc_header()
    _static_block += report_micapipe()
    
    print( 'Done')
    return _static_block

# --------------------------------------------------------------------
# Write pdf report
# --------------------------------------------------------------------
def convert_html_to_pdf(source_html, output_filename):
    # open output file for writing (truncated binary)
    result_file = open(output_filename, "w+b")

    # convert HTML to PDFt
    pisa_status = pisa.CreatePDF(
            source_html,                # the HTML to convert
            dest=result_file)           # file handle to recieve result

    # close output file
    result_file.close()                 # close output file
    print('Out pdf group level report: '+output_filename)
    # return True on success and False on errors
    return pisa_status.err

# Generate PDF report of Micapipe QC
convert_html_to_pdf(qc_group(), f'{out}/micapipe_group-QC.pdf')

# ------------------------------------------
# Delete the temporary directory and its contents
tmp_files=sorted(glob.glob(tmpDir+'/*'))
for x in tmp_files: os.remove(x)
os.rmdir(tmpDir)