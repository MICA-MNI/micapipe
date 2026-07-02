"""Tests for functions/nlm_denoise.py (mincnlm -> NIfTI/Python port).

Unit tests exercise the CLI and guards with no heavy deps. Integration tests
are gated on dipy (`importorskip`) and run the real Coupé ONLM filter on a
synthetic phantom, asserting the noise actually goes down.

Run:  pytest tests/test_nlm_denoise.py -v
"""

import importlib.util
import os
import sys
import warnings

import numpy as np
import nibabel as nib
import pytest

# Load functions/nlm_denoise.py by path (functions/ is not an importable pkg).
_HERE = os.path.dirname(os.path.abspath(__file__))
_MODPATH = os.path.join(_HERE, "..", "functions", "nlm_denoise.py")
_spec = importlib.util.spec_from_file_location("nlm_denoise", _MODPATH)
nlm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nlm)


# --------------------------------------------------------------------------- #
# Fixtures / helpers
# --------------------------------------------------------------------------- #
def _phantom(shape=(24, 24, 24), fg=100.0):
    """A simple 3D box phantom: bright cube on a dark background."""
    vol = np.zeros(shape, dtype=np.float32)
    s = tuple(slice(d // 4, 3 * d // 4) for d in shape)
    vol[s] = fg
    return vol


def _add_gaussian_noise(vol, sigma, seed=0):
    rng = np.random.default_rng(seed)
    return (vol + rng.normal(0, sigma, vol.shape)).astype(np.float32)


def _save_nii(path, data, affine=None):
    if affine is None:
        affine = np.diag([1.5, 1.5, 1.5, 1.0])
    nib.save(nib.Nifti1Image(data.astype(np.float32), affine), str(path))


# --------------------------------------------------------------------------- #
# Unit tests — no dipy required
# --------------------------------------------------------------------------- #
def test_parser_defaults_match_mincnlm():
    args = nlm.build_parser().parse_args(["in.nii", "out.nii"])
    assert args.sigma == 0.0          # 0 = auto estimate
    assert args.beta == 1.0
    assert args.patch_radius == 1     # -v
    assert args.block_radius == 5     # -d
    assert args.noise_model == 0      # gaussian
    assert args.blockwise is True     # -block is the default
    assert args.threads is None


def test_noise_model_mapping():
    assert nlm.noise_model_is_rician(0) is False
    assert nlm.noise_model_is_rician(1) is False
    assert nlm.noise_model_is_rician(2) is True


def test_block_and_pointwise_are_mutually_exclusive():
    assert nlm.build_parser().parse_args(["i", "o", "--pointwise"]).blockwise is False
    with pytest.raises(SystemExit):
        nlm.build_parser().parse_args(["i", "o", "--block", "--pointwise"])


def test_invalid_noise_model_rejected():
    with pytest.raises(SystemExit):
        nlm.build_parser().parse_args(["i", "o", "--noise-model", "7"])


def test_clobber_guard_blocks_existing_output(tmp_path):
    inp = tmp_path / "in.nii.gz"
    out = tmp_path / "out.nii.gz"
    _save_nii(inp, _phantom())
    out.write_bytes(b"existing")
    with pytest.raises(SystemExit) as e:
        nlm.main([str(inp), str(out)])
    assert "exists" in str(e.value)


def test_rejects_non_3d_input(tmp_path):
    inp = tmp_path / "4d.nii.gz"
    out = tmp_path / "out.nii.gz"
    _save_nii(inp, np.zeros((8, 8, 8, 4), dtype=np.float32))
    with pytest.raises(SystemExit) as e:
        nlm.main([str(inp), str(out)])
    assert "3D" in str(e.value)


def test_beta_not_one_warns_before_touching_dipy(tmp_path):
    # A 4-D input makes main() bail right after the beta check, so this test
    # verifies the warning without needing dipy installed.
    inp = tmp_path / "4d.nii.gz"
    out = tmp_path / "out.nii.gz"
    _save_nii(inp, np.zeros((8, 8, 8, 2), dtype=np.float32))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with pytest.raises(SystemExit):
            nlm.main([str(inp), str(out), "--beta", "0.5"])
    assert any("beta" in str(x.message) for x in w)


# --------------------------------------------------------------------------- #
# Integration tests — real denoising, gated on dipy
# --------------------------------------------------------------------------- #
def test_denoise_reduces_noise_and_preserves_geometry(tmp_path):
    pytest.importorskip("dipy")
    clean = _phantom()
    noisy = _add_gaussian_noise(clean, sigma=15.0)

    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    inp = tmp_path / "noisy.nii.gz"
    out = tmp_path / "denoised.nii.gz"
    _save_nii(inp, noisy, affine)

    nlm.main([str(inp), str(out), "--sigma", "15", "--verbose"])

    assert out.exists()
    res = nib.load(str(out))
    got = res.get_fdata(dtype=np.float32)

    # Geometry preserved.
    assert got.shape == clean.shape
    np.testing.assert_allclose(res.affine, affine, atol=1e-6)

    # The denoised image is closer to the clean phantom than the noisy input.
    rmse = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))
    assert rmse(got, clean) < rmse(noisy, clean)


def test_auto_sigma_estimation_runs(tmp_path):
    pytest.importorskip("dipy")
    noisy = _add_gaussian_noise(_phantom(), sigma=12.0)
    inp = tmp_path / "noisy.nii.gz"
    out = tmp_path / "out.nii.gz"
    _save_nii(inp, noisy)
    # sigma defaults to 0 -> automatic estimation path.
    nlm.main([str(inp), str(out)])
    assert out.exists()


def test_rician_and_pointwise_paths_run(tmp_path):
    pytest.importorskip("dipy")
    noisy = _add_gaussian_noise(_phantom(), sigma=10.0)
    inp = tmp_path / "noisy.nii.gz"
    out = tmp_path / "out.nii.gz"
    _save_nii(inp, noisy)
    nlm.main([str(inp), str(out), "--sigma", "10",
              "--noise-model", "2", "--pointwise"])
    assert out.exists()
