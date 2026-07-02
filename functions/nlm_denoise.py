#!/usr/bin/env python3
"""
NIfTI Non-Local Means denoising — Python replacement for EZminc's `mincnlm`.

`mincnlm` (BIC-MNI/EZminc) is a C++ wrapper around Coupé et al.'s optimized
blockwise Non-Local Means (ONLM) filter for 3D MRI, reading/writing MINC.
This module reproduces the same algorithm on NIfTI inputs using nibabel for
I/O and dipy for the denoiser, which ships the same ONLM implementation
(Descoteaux/Coupé port) plus automatic noise estimation.

Parameter mapping (mincnlm flag -> here):
    -sigma   -> --sigma          0 = automatic estimation, as in mincnlm
    -v       -> --patch-radius   patch/neighborhood radius   (default 1)
    -d       -> --block-radius   search-volume radius        (default 5)
    -w       -> --noise-model    0/1 = gaussian, 2 = rician   (default gaussian)
    -block   -> --block / --pointwise  (blockwise ONLM is the mincnlm default)
    -mt      -> --threads

Note on `-beta`: mincnlm exposes beta as the smoothing-bandwidth scale
(h ~ beta * sigma). dipy's non_local_means fixes this internally, so beta is
accepted for CLI compatibility but only 1.0 is honored; a warning is emitted
otherwise.

Requires: nibabel (already in micapipe env) and dipy (must be added to
micapipe_environment.yml — it is not currently installed there).

dipy is imported lazily inside the functions that need it, so this module can
be imported (and its CLI/guards unit-tested) in environments without dipy.
"""

import argparse
import os
import sys
import warnings

import numpy as np
import nibabel as nib


def build_parser():
    """Build the argument parser (pure; no dipy needed)."""
    p = argparse.ArgumentParser(
        description="Non-Local Means denoising for NIfTI (mincnlm replacement).",
    )
    p.add_argument("input", help="input NIfTI image (.nii/.nii.gz)")
    p.add_argument("output", help="output denoised NIfTI image")
    p.add_argument(
        "--sigma", type=float, default=0.0,
        help="noise standard deviation; 0 = estimate automatically (mincnlm -sigma)",
    )
    p.add_argument(
        "--beta", type=float, default=1.0,
        help="smoothing bandwidth scale (mincnlm -beta); dipy honors 1.0 only",
    )
    p.add_argument(
        "--patch-radius", type=int, default=1,
        help="patch/neighborhood radius (mincnlm -v)",
    )
    p.add_argument(
        "--block-radius", type=int, default=5,
        help="search-volume radius (mincnlm -d)",
    )
    p.add_argument(
        "--noise-model", type=int, choices=(0, 1, 2), default=0,
        help="0/1 = gaussian (L2), 2 = rician (mincnlm -w)",
    )
    grp = p.add_mutually_exclusive_group()
    grp.add_argument(
        "--block", dest="blockwise", action="store_true", default=True,
        help="blockwise ONLM (mincnlm default)",
    )
    grp.add_argument(
        "--pointwise", dest="blockwise", action="store_false",
        help="voxelwise NLM instead of blockwise",
    )
    p.add_argument("--mask", help="optional NIfTI brain mask")
    p.add_argument(
        "--threads", type=int, default=None,
        help="number of threads (mincnlm -mt); default = all cores",
    )
    p.add_argument(
        "--clobber", action="store_true",
        help="overwrite the output if it already exists",
    )
    p.add_argument("--verbose", action="store_true")
    return p


def noise_model_is_rician(w):
    """mincnlm -w: 0/1 = gaussian (L2), 2 = rician."""
    return w == 2


def resolve_sigma(data, sigma_arg):
    """Return the noise sigma: the given value, or an auto-estimate if <= 0.

    Mirrors mincnlm's `-sigma 0` -> automatic estimation.
    """
    if sigma_arg > 0:
        return sigma_arg
    from dipy.denoise.noise_estimate import estimate_sigma
    return estimate_sigma(data, N=0)


def run_denoise(data, sigma, patch_radius, block_radius, rician,
                blockwise=True, mask=None, threads=None):
    """Apply (blockwise or voxelwise) NLM. dipy imported lazily.

    dipy's `nlmeans` implements Coupé et al.'s ONLM. Newer dipy (>=1.10) unifies
    it under a `method` kwarg (`'blockwise'` = mincnlm default, `'classic'` =
    voxelwise); older dipy splits it into `non_local_means` (blockwise) and
    `nlmeans` (voxelwise). Support both so the container's pinned version works.
    """
    import inspect
    from dipy.denoise.nlmeans import nlmeans

    common = dict(
        sigma=sigma, mask=mask,
        patch_radius=patch_radius, block_radius=block_radius,
        rician=rician, num_threads=threads,
    )
    if "method" in inspect.signature(nlmeans).parameters:
        return nlmeans(data, method="blockwise" if blockwise else "classic",
                       **common)
    if blockwise:
        try:
            from dipy.denoise.non_local_means import non_local_means
            return non_local_means(data, **common)
        except ImportError:
            pass  # fall through to voxelwise nlmeans
    return nlmeans(data, **common)


def main(argv=None):
    args = build_parser().parse_args(argv)

    if os.path.exists(args.output) and not args.clobber:
        sys.exit(f"error: {args.output} exists (use --clobber to overwrite)")

    if args.beta != 1.0:
        warnings.warn(
            "beta != 1.0 is not supported by dipy's non_local_means; "
            "proceeding with the built-in bandwidth (beta=1.0)."
        )

    img = nib.load(args.input)
    data = img.get_fdata(dtype=np.float32)
    if data.ndim != 3:
        sys.exit(f"error: expected a 3D image, got shape {data.shape}")

    mask = None
    if args.mask:
        mask = np.asarray(nib.load(args.mask).get_fdata() > 0, dtype=bool)

    rician = noise_model_is_rician(args.noise_model)
    sigma = resolve_sigma(data, args.sigma)
    if args.verbose:
        print(f"sigma: {np.mean(np.atleast_1d(sigma)):.4f}")
        print(
            f"NLM ({'blockwise' if args.blockwise else 'voxelwise'}), "
            f"patch_radius={args.patch_radius}, block_radius={args.block_radius}, "
            f"rician={rician}"
        )

    denoised = run_denoise(
        data, sigma=sigma,
        patch_radius=args.patch_radius, block_radius=args.block_radius,
        rician=rician, blockwise=args.blockwise, mask=mask, threads=args.threads,
    )

    out = nib.Nifti1Image(denoised.astype(np.float32), img.affine, img.header)
    nib.save(out, args.output)
    if args.verbose:
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
