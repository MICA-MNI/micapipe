#!/usr/bin/env python3
"""
NIfTI Non-Local Means denoising via ANTs `DenoiseImage` — a C++ alternative to
the dipy-based `nlm_denoise.py`, for comparison / structural use.

ANTs `DenoiseImage` implements Manjón/Coupé adaptive Non-Local Means in C++
(the same algorithm family as EZminc's mincnlm), is NIfTI-native, fast, ships
in the micapipe base image, and is widely validated on structural MRI. This
wrapper mirrors `nlm_denoise.py`'s CLI so the two engines are drop-in
comparable.

Flag mapping (nlm_denoise flag -> DenoiseImage):
    --noise-model 2 -> -n Rician ;  0/1 -> -n Gaussian
    --patch-radius  -> -p   (patch radius)
    --block-radius  -> -r   (search radius; ANTs default is 2)
    --mask          -> -x   (mask image)
    --threads       -> ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
Notes:
    --sigma is accepted for CLI parity but ignored — DenoiseImage estimates the
    noise level internally (a warning is emitted if a value is passed).

Requires: ANTs `DenoiseImage` on PATH (present in the micapipe container).
"""

import argparse
import os
import shutil
import subprocess
import sys
import warnings


class DenoiseError(ValueError):
    """Raised on invalid user input; main() turns it into a clean exit."""


def _positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1, got {value}")
    return value


def build_parser():
    p = argparse.ArgumentParser(
        description="NLM denoising for NIfTI via ANTs DenoiseImage (C++).",
    )
    p.add_argument("input", help="input NIfTI image (.nii/.nii.gz)")
    p.add_argument("output", help="output denoised NIfTI image")
    p.add_argument(
        "--sigma", type=float, default=0.0,
        help="accepted for parity with nlm_denoise; ignored (ANTs auto-estimates)",
    )
    p.add_argument(
        "--patch-radius", type=_positive_int, default=1, help="ANTs -p (patch radius)",
    )
    p.add_argument(
        "--block-radius", type=_positive_int, default=2, help="ANTs -r (search radius)",
    )
    p.add_argument(
        "--noise-model", type=int, choices=(0, 1, 2), default=0,
        help="0/1 = Gaussian, 2 = Rician",
    )
    p.add_argument("--mask", help="optional NIfTI mask (ANTs -x)")
    p.add_argument("--threads", type=int, default=None, help="ITK thread count")
    p.add_argument("--clobber", action="store_true", help="overwrite existing output")
    p.add_argument("--verbose", action="store_true")
    return p


def build_command(args):
    """Assemble the DenoiseImage argv (pure; no ANTs needed for this)."""
    cmd = [
        "DenoiseImage", "-d", "3",
        "-i", args.input, "-o", args.output,
        "-n", "Rician" if args.noise_model == 2 else "Gaussian",
        "-p", str(args.patch_radius),
        "-r", str(args.block_radius),
    ]
    if args.mask:
        cmd += ["-x", args.mask]
    if args.verbose:
        cmd += ["-v", "1"]
    return cmd


def denoise_file(args, runner=subprocess.run):
    """Validate inputs and run DenoiseImage. `runner` is injectable for tests."""
    if not os.path.isfile(args.input):
        raise DenoiseError(f"input not found: {args.input}")
    if os.path.exists(args.output) and not args.clobber:
        raise DenoiseError(f"{args.output} exists (use --clobber to overwrite)")
    if args.mask and not os.path.isfile(args.mask):
        raise DenoiseError(f"mask not found: {args.mask}")
    if args.sigma:
        warnings.warn("--sigma is ignored by ANTs DenoiseImage (auto-estimated).")
    if shutil.which("DenoiseImage") is None:
        raise DenoiseError("DenoiseImage not on PATH (need ANTs / micapipe env).")

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)

    env = dict(os.environ)
    if args.threads is not None:
        env["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(args.threads)

    cmd = build_command(args)
    if args.verbose:
        print("running:", " ".join(cmd))
    result = runner(cmd, env=env)
    if getattr(result, "returncode", 0) != 0:
        raise DenoiseError(f"DenoiseImage failed (exit {result.returncode})")


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        denoise_file(args)
    except DenoiseError as exc:
        sys.exit(f"error: {exc}")


if __name__ == "__main__":
    main()
