# Alternative FastSurfer surfaces with FastScrubber (`-alt_fs`)

**Experimental, intended mainly for 7T data.** `-alt_fs` changes the standard
FastSurfer workflow and has not been validated as thoroughly as the defaults:
FastSurfer (default) and FreeSurfer (`-freesurfer`) are unchanged.

```
micapipe -sub SUB -ses SES -bids BIDS -out OUT -fs_licence LICENSE \
    -proc_structural -proc_surf -alt_fs -post_structural
```

## What it does

1. FastSurfer runs as usual (or a FastSurfer directory given with `-surf_dir`
   is copied into the output; the original is never modified).
2. [FastScrubber](https://github.com/donnagift/fast_scrubber), an nnU-Net
   trained on 200 manually corrected 7T MP2RAGE brain masks, predicts a new
   brain mask from FastSurfer's `mri/orig.mgz`, excluding the meninges,
   vessels, tentorium and sinuses that FastSurfer's mask often keeps at 7T.
3. FastScrubber only decides the outer boundary of the brain. micapipe adds
   back every voxel FastSurfer's segmentation (`mri/aseg.auto.mgz`) labels as
   deep brain (all labels except cerebral and cerebellar cortex and CSF),
   keeps the largest connected piece and fills interior holes, then uses the
   result as `mri/mask.mgz`. On a 3T scan FastScrubber left a detached 75 mL
   non-brain piece, which makes `mri_normalize` fail outright.
   It then clears every label outside the new mask in `mri/aseg.auto.mgz` and
   `mri/aseg.presurf.mgz`, as FreeSurfer's aseg never leaves the brain mask.
   `mri_normalize` uses the aseg white matter as control points: a few hundred
   white matter voxels outside the mask give it a NaN white matter peak and
   the surfaces collapse (BIDS_CI sub-mri7T: mean thickness 2.15 -> 0.77 mm
   with 252 such voxels). FastSurfer 2.6.0 also labels a detached piece about
   110 mm below the brain as white matter on BIDS_CI sub-mri3T, which collapses
   the surfaces on the default path as well; reported upstream as
   [Deep-MI/FastSurfer#956](https://github.com/Deep-MI/FastSurfer/issues/956),
   where 2.5.4 behaves the same way and 2.4.2 does not. FastSurfer's
   segmentation stays unchanged in `mri/aparc.DKTatlas+aseg.deep.mgz`.
   FastSurfer's own mask is kept as `mri/mask.fastsurfer.mgz`, FastScrubber's
   raw prediction as `mri/mask.fastscrubber.mgz`, and the first surface logs as
   `scripts/recon-surf.fastsurfer.log` and `scripts/recon-all.fastsurfer.log`.
4. FastSurfer rebuilds the surfaces (`run_fastsurfer.sh --surf_only`), which
   recreates `norm.mgz` and `brainmask.mgz` from the new mask.
5. `-post_structural` and later modules use these surfaces like any FastSurfer
   output. `scripts/micapipe_alt_fs.json` records the options used.

The surfaces stay in the `fastsurfer/` derivatives directory.

## Inference settings

FastScrubber runs all 5 cross-validation folds **without** test-time
augmentation (TTA) on the GPU or the CPU. `-alt_fs_tta` adds TTA (8 mirrored
passes per fold). Measured on held-out data, dropping TTA changes very little
and saves most of the time:

| | 7T PNI (10 held-out sessions) | MICs 3T (6 sessions) |
| --- | --- | --- |
| Dice, 5 folds without vs with TTA | 0.9988-0.9992 | 0.992-0.999 |
| Time without / with TTA, GPU (RTX 2060 SUPER) | 2.4 / 13.5 min | |
| Time without / with TTA, 8 CPU threads | | 64-68 min / ~9.4 h |

Single folds are not offered: on 3T data, which FastScrubber was not trained
on, one fold left a 125 mL non-brain piece in one of six sessions.

## Containers

FastScrubber, its environment and its weights are included in both micapipe
v1 beta containers. The `cuda` variant runs FastSurfer and FastScrubber on an
NVIDIA GPU when one is available (`singularity run --nv`,
`docker run --gpus all`); the `nocuda` variant runs on the CPU only. Outside
the containers, set `MICAPIPE_FASTSCRUBBER` to the FastScrubber command, for
example `"python /path/to/fastscrubber.py"`.

## Known caveats

- FastScrubber masks can have interior holes. In 10 MICs 3T sessions processed
  with the raw masks, the ventricles shrank by about 2% and pericalcarine/lingual
  topology defects appeared in about 3 of 10 sessions, while cortical thickness
  error improved. micapipe now restores deep brain labels and fills interior
  holes, which has not yet been re-evaluated on that cohort.
- The model was trained on 7T MP2RAGE UNI images only. 3T results are not
  validated against manual masks.
- Review the surfaces visually before using them.
