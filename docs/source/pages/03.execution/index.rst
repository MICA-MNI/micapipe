.. _execution:

.. title:: How execute the micapipe

Execution of micapipe
======================================


Command-Line Arguments
----------------------
.. code-block:: text

        $ mica-pipe  -sub <subject_id> -out <outputDirectory> -bids <BIDS-directory>

Positional Arguments
~~~~~~~~~~~~~~~~~~~~~~
	↪ **<subject_id>** : the subject identification (without the "sub-")

	↪ **<outputDirectory>** : the output path for the outcomes files of the preprocessing

	↪ **<BIDS-directory>** : the input path, where the folder of your BIDS valid dataset is

Some essential flags 
~~~~~~~~~~~~~~~~~~~~~~
	↪ **-ses** : Number of session (Default is ses-01)

	↪ **-force** : To overwrite the subject directory (WARNING! This will suppress your subject directory)

	↪ **-quiet** : Do NOT print comments

	↪ **-nocleanup** : Do NOT delete temporal directory at script completion

Flags for processing
~~~~~~~~~~~~~~~~~~~~~~
🚩 Flags for first stages of structural processing: 

	↪ **-proc_structural** : Volumetric processing

	↪ **-proc_freesurfer** : Freesurfer recon-all processing

	↪ **-proc_dwi** : DWI-Diffusion weighted images processing with MRtrix3

	↪ **-proc_rsfmri** : Resting state Functional MRI processing 

.. admonition:: Important to know ☝🏼

     You can use 🚩 -proc to run all the first stages of micapipe 	
	 The first stages of structural processing correspond to all the -proc stages of micapipe. These steps are used to preprocess the images in order to make them usable for the -post stages. 	

🚩 Flags for second stages of structural processing:

	↪ **-post_structural** : Post structural volumetric processing

	↪ **-post_dwi** : Post tractography and connectome generation 

	↪ **-post_mpc** : Microstructural profiles and covariance analysis

.. admonition:: Important to know ☝🏼

     You can use 🚩 -post to run all the second stages of micapipe 	
	 The second stages of structural processing correspond to all the -post stages of micapipe. These steps generate connectomes, correlations and matrices. 	


