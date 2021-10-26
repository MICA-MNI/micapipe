.. _autotract:

.. title:: Automatic Bundle Segmentation

Automatic Bundle Segmentation
================================================

This section describes the integrated automatic virtual dissection of the main brain tracts, which is based in `auto_tracto <https://github.com/rcruces/auto_tracto/blob/master/Readme.md>`_, and was implemented using ANTs registration tools.
The script ``03_auto_tracts.sh`` main purpose is to split a tck tractogram into `35 white matter tracts <https://github.com/lconcha/auto_tracto/blob/master/tract_definitions.md>`_, using already stablished automatic dissection protocols, which were manually tuned for optimal performance.

Derived from a full brain tractogram, 35 bundles are virtually dissected using the `LANIREM protocols <https://github.com/rcruces/auto_tracto/tree/master/lanirem/protocols>`_. The quality of the full brain tractogram will determine the quality of bundle separation. It is highly recommended to provide a tractogram with more than one million streamlines, and one that has been checked for errors. Strategies such as anatomically-contstrained tractography (ACT) and spherical deconvolution informed filtering of tractograms (SIFT), available in MRTrix3, should aid in obtaining such high-quality tractograms.

.. image:: autotract.png
   :align: center

.. tabs::
       .. tab:: Processing steps

          - Non-linear (SyN) registration of the native FA map to the FA atlas (FMRIB58_FA_1mm)
          - Apply transformation to each bundle protocols to register them to the native FA space (DWI).
          - Verify stop criterion conditionals.
          - Filter each white matter bundle according to the dissection protocols.

       .. tab:: Usage

          There are two ways to run the *automatic bundle segmentation*, one is integrated within the main script ``mica-pipe`` under the ``-SC`` module, and second is using the stand alone script ``03_auto_tracts.sh``.

          .. code-block:: bash
            :linenos:
            :caption: Stand alone script requires that all inputs are in DWI space.

            03_auto_tracts.sh \
                 -tck sub-01_full_brain_tractogram.tck \
                 -outbase sub-01_tract \
                 -mask sub-01_binary_brain_mask.nii.gz \
                 -fa sub-01_FA_map.nii.gz \
                 -weights sub-01_full_brain_tractogram_weights.txt

          .. code-block:: bash
            :linenos:
            :caption: Integrated usage within ``mica-pipe -SC``
            :emphasize-lines: 6

            mica-pipe  \
                 -sub <subject> \
                 -ses <session> \
                 -out <outputDirectory> \
                 -bids <BIDS> \
                 -SC -autoTract \
                 -tracts 40M

          **Mandatory Arguments**
            -tck        ``path`` full tractogram file *tck* (ideally SIFTED).
            -outbase    ``str`` base name for all your outputs.
            -mask       ``path`` Binary mask in subject dwi space.
            -fa         ``path`` FA map in subject dwi space. Used for registration to template.

          **Optional Arguments**
            -keep_tmp                 Do not delete temporary directory
            -tmpDir                   Specify location of temporary directory.
            -minStreamlinesPerVoxel   ``num`` Streamlines are truncated if voxel contains less than this number of streamlines. Default is 1
            -robust                   This option to runs a more ROBUST SyN registration ( More computation time )
            -weights                  ``path`` Use this option if you calculated a weights file from SIFT2
            -threads                  ``num`` Number of threads (Default is 6)
            -version      	          Print software version

       .. tab:: Outputs

          Directories created by this script will be in the selected ``outbase``.
          The outputs generate by the micapipe integrated script are in *<outputDirectory>/micapipe/<sub>/dwi/auto_tract*

          .. parsed-literal::
              dwi/auto_tract/
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_AC.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_AF_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_AF_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_atlas2fa.nii.gz
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CC_MID.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CGFP_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CGFP_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CGH_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CGH_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CG_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CGR_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CGR_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CG_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CST_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_CST_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_FA_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_FA_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_FMA.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_FMI.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_FX_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_FX_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_IFOF_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_IFOF_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_ILF_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_ILF_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_MLF_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_MLF_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_OR_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_OR_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_SLF_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_SLF_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_summary.txt
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_TAPETUM.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_UF_L.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_UF_R.tck
              ├── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_VOF_L.tck
              └── <sub>_space-dwi_desc-iFOD2-40M-SIFT2_VOF_R.tck

          The description refers to the algorithm used to generate the tractogram (*iFOD2*), its number of streamlines (*40M*) and the filtering algorithm (*SIFT2*).

          +---------+--------------------------------------+
          | Acronym | Tract name                           |
          +=========+======================================+
          | AC      | Anterior Commissure                  |
          +---------+--------------------------------------+
          | AF      | Arcuate Fasciculus                   |
          +---------+--------------------------------------+
          | CC_MID  | Corpus Callosum middle portion       |
          +---------+--------------------------------------+
          | CGFP    | Cingulum, fronto-parietal portion    |
          +---------+--------------------------------------+
          | CGH     | Cingulum, parahippocampal portion    |
          +---------+--------------------------------------+
          | CG      | Cingulum, whole                      |
          +---------+--------------------------------------+
          | CGR     | Cingulum, rostral and subgenual      |
          +---------+--------------------------------------+
          | CST     | Corticospinal tract                  |
          +---------+--------------------------------------+
          | FA      | Frontal aslant                       |
          +---------+--------------------------------------+
          | FMA     | Forceps major of corpus callosum     |
          +---------+--------------------------------------+
          | FMI     |Forceps minor of corpus callosum      |
          +---------+--------------------------------------+
          | FX      | Fornix                               |
          +---------+--------------------------------------+
          | IFOF    | Inferior fronto-occipital fasciculus |
          +---------+--------------------------------------+
          | ILF     | Inferior longitudinal fasciculus     |
          +---------+--------------------------------------+
          | MLF     | Middle longitudinal fasciculus       |
          +---------+--------------------------------------+
          | OR      | Optic radiation                      |
          +---------+--------------------------------------+
          | SLF     | Superior longitudinal fasciculus     |
          +---------+--------------------------------------+
          | UF      | Uncinate fasciculus                  |
          +---------+--------------------------------------+
          | VOF     | Vertical occipital fasciculus        |
          +---------+--------------------------------------+
