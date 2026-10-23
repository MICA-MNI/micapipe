# V1 environment beta validation

The CPU engineering beta integrates FastSurfer 2.6.0, FreeSurfer 8.2.0,
MRtrix3 3.0.8, ANTs 2.6.5, FSL 6.0.7.8, pyonlm and NiiVue 0.69.0.
Each tool uses a separate pinned Python environment.

The existing FIX training model was converted with the official FSL
convert_fix_model command. A three-case MICs compatibility pilot compared both
classifiers using identical legacy features and also ran the full new extractor.
At threshold 20, label agreement ranged from 94.19% to 98.99%, and temporal
correlation after voxelwise demeaning ranged from 0.951763 to 0.998045.
None passed the combined engineering gate (label agreement >=95%, noise Jaccard
>=0.90, temporal correlation >=0.999, temporal NRMSE <=1%).
This is a compatibility test; held-out status and scientific accuracy were not established.

**Legacy FIX remains the default. R and Octave are retained.**
Conversion success does not establish scientific equivalence.

The beta.1 image completed native-version/import checks, upstream pyonlm tests,
offline CPU FastSurfer segmentation, structural processing, and offline NiiVue tests.
Further testing found that structural QC PDFs could omit anatomical panels.
Use the [beta.2 follow-up](v1-beta-2.md), which fixes capture and failure handling.

Subject-level records, anatomical previews and internal infrastructure details
are retained with the lab's server validation evidence and are excluded from
this public branch.
