#!/usr/bin/env python3
"""
NIfTI Non-Local Means denoising — Python replacement for EZminc's `mincnlm`.

`mincnlm` (BIC-MNI/EZminc) is a C++ wrapper around Coupé et al.'s optimized
blockwise Non-Local Means (ONLM) filter for 3D MRI, reading/writing MINC.
This module reproduces the same algorithm on NIfTI inputs using nibabel for
I/O and dipy for the denoiser, which ships the same ONLM implementation
(Descoteaux/Coupé port) plus automatic noise estimation.

Replacing MINC with NIfTI is the whole point: inputs and outputs are `.nii` /
`.nii.gz`, and the input affine + header geometry are carried through to the
output unchanged.

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

Requires: nibabel (already in micapipe env) and dipy (added to
micapipe_environment.yml alongside this tool).

dipy is imported lazily inside the functions that need it, so this module can
be imported (and its CLI/guards unit-tested) in environments without dipy.
"""

import argparse
import os
import sys
import warnings

import numpy as np
import nibabel as nib


class DenoiseError(ValueError):
    """Raised on invalid user input; main() turns it into a clean exit."""


def _positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {value}")
    return value


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
        "--patch-radius", type=_positive_int, default=1,
        help="patch/neighborhood radius (mincnlm -v)",
    )
    p.add_argument(
        "--block-radius", type=_positive_int, default=5,
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
    p.add_argument("--mask", help="optional NIfTI brain mask (same grid as input)")
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


def load_3d_volume(path):
    """Load a NIfTI file as a float32 3D array + its image object.

    Raises DenoiseError (not a raw traceback) on a missing/unreadable file or
    a non-3D image, so the CLI reports a clean message.
    """
    if not os.path.isfile(path):
        raise DenoiseError(f"input not found: {path}")
    try:
        img = nib.load(path)
    except Exception as exc:  # nibabel raises assorted types on bad files
        raise DenoiseError(f"could not read {path} as NIfTI: {exc}")
    data = np.asanyarray(img.dataobj, dtype=np.float32)
    if data.ndim != 3:
        raise DenoiseError(f"expected a 3D image, got shape {data.shape}")
    return img, data


def sanitize(data):
    """Replace non-finite voxels (NaN/Inf) with 0, warning if any were found."""
    finite = np.isfinite(data)
    if not finite.all():
        n = int((~finite).sum())
        warnings.warn(f"{n} non-finite voxel(s) set to 0 before denoising.")
        data = np.where(finite, data, 0.0).astype(np.float32)
    return data


def load_mask(path, shape):
    """Load a boolean mask and verify it lives on the same grid as the data."""
    if not os.path.isfile(path):
        raise DenoiseError(f"mask not found: {path}")
    mask = np.asanyarray(nib.load(path).dataobj)
    if mask.shape != shape:
        raise DenoiseError(
            f"mask shape {mask.shape} does not match image shape {shape}"
        )
    return mask > 0


def resolve_sigma(data, sigma_arg):
    """Return the noise sigma: the given value, or an auto-estimate if 0.

    Mirrors mincnlm's `-sigma 0` -> automatic estimation. A negative sigma is
    an error (only 0 means "estimate").
    """
    if sigma_arg < 0:
        raise DenoiseError(f"sigma must be >= 0 (0 = auto), got {sigma_arg}")
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


def build_output(denoised, template_img):
    """Wrap the denoised array in a NIfTI image that keeps the input geometry.

    The affine and header are copied from the input; the header's on-disk dtype
    is set to float32 so it matches the data we actually write.
    """
    header = template_img.header.copy()
    header.set_data_dtype(np.float32)
    return nib.Nifti1Image(denoised.astype(np.float32), template_img.affine, header)


def denoise_file(args):
    """End-to-end run from parsed args. Raises DenoiseError on bad input."""
    if os.path.exists(args.output) and not args.clobber:
        raise DenoiseError(f"{args.output} exists (use --clobber to overwrite)")

    if args.beta != 1.0:
        warnings.warn(
            "beta != 1.0 is not supported by dipy's non_local_means; "
            "proceeding with the built-in bandwidth (beta=1.0)."
        )

    img, data = load_3d_volume(args.input)
    data = sanitize(data)

    mask = load_mask(args.mask, data.shape) if args.mask else None
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

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    nib.save(build_output(denoised, img), args.output)
    if args.verbose:
        print(f"wrote {args.output}")


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        denoise_file(args)
    except DenoiseError as exc:
        sys.exit(f"error: {exc}")


if __name__ == "__main__":
    main()
