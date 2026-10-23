# V1 beta.2 follow-up validation

Reviewed 2026-10-23. All builds and tests ran on Linux servers.
This is a CPU engineering beta for lab review.

## Corrected QC behavior

Beta.1 could return success when a non-executable image helper failed, producing
text-only structural QC reports. Beta.2 calls the helper through the selected
Python interpreter, checks capture/PDF failures, and preserves the failure
status through the QC wrapper and main pipeline. It also preserves an explicit
FreeSurfer license path and corrects the structural image-width CSS.

Two real structural reports were regenerated with the final SIF: each had
three pages, eight anatomical panels and nine embedded images.
A deliberately failed capture returned CLI exit 1.
All six pages were rendered on the server. Manual visual sign-off remains pending.

## Completed engineering checks

- Final-image native-version, dependency-isolation and real PNG-generation
  smoke checks passed on two Linux servers, including execution from shared storage.
- Both packaged FIX classifiers reproduced their respective reference labels.
- Fresh structural processing completed 09/09 stages inside the beta.1 SIF;
  its 5TT check passed. Beta.2 uses the same native binaries and regenerated its QC.
- A second structural pilot completed 09/09 with two tissue-fraction voxels
  flagged for manual review; no probabilities were modified.
- Real CPU LAMAReg/SynthSeg inference produced finite, nonempty integer labels.
- The packaged pyonlm CLI denoised a real NIfTI crop with finite output and
  preserved geometry; all 25 upstream tests passed.
- Offline FastSurfer CPU segmentation and offline NiiVue viewing passed.

## Diffusion pilot

A 41-volume acquisition completed MRtrix denoising/Gibbs correction, FSL
topup/eddy_cpu, ANTs bias correction, tensor fitting and eddy_quad.
Corrected data, FA and MD were finite, and a six-page eddy QC PDF was generated.
This standalone pilot is not the complete micapipe DWI workflow.

The physical-range check flagged 83 of 363,767 masked voxels with FA above 1
(maximum 1.22443), and 14 with negative MD. All FA outliers had a negative
tensor eigenvalue. A controlled refit of the same corrected data and mask with
MRtrix 3.0.1 from the older image reproduced FA exactly (RMSE 0), including all
outliers and negative MD values. The old tensor2metric also reproduced the new
tensor's FA exactly. This found no tensor-fit/metric regression in that case;
it does not validate preprocessing equivalence or remove the QC warning.
No values were clipped or changed.

## Limits and provenance

The [FIX compatibility gate](v1-beta.md) failed; legacy FIX, R and Octave remain.
Full surface reconstruction, GPU execution, the complete functional/DWI
workflows, full SnakeBIDS execution, manual anatomical QC and production/scientific
readiness are not established.

Python locks and requested native versions/archive hashes are pinned.
OS apt packages are not snapshot-pinned, so bit-for-bit rebuilds are not claimed.

The released beta.2 SIF was built from archived server implementation revision
2103d4630daa60fec18399bce618a525f8009aab.
SHA-256: d38144ef090da8b794d45eeab4f71f4e46f369d2618d9f938f2d5281a81ee973.
This public snapshot preserves the tested processing code, assets and locks.
Only the validation documents were abridged for publication. Detailed subject
records, previews and full development history remain in the lab's server archive.
