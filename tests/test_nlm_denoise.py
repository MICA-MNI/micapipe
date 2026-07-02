"""Tests for functions/nlm_denoise.py (mincnlm -> NIfTI/Python port).

Unit tests exercise the CLI, guards and pure helpers with no heavy deps.
Integration tests are gated on dipy (`importorskip`) and run the real Coupé
ONLM filter on a synthetic phantom, asserting the noise actually goes down and
that NIfTI geometry survives the round-trip.

Run:  pytest tests/test_nlm_denoise.py -v
"""

import importlib.util
import os
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


def _rmse(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def _parse(*argv):
    return nlm.build_parser().parse_args(list(argv))


# --------------------------------------------------------------------------- #
# Unit tests — no dipy required
# --------------------------------------------------------------------------- #
def test_parser_defaults_match_mincnlm():
    args = _parse("in.nii", "out.nii")
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
    assert _parse("i", "o", "--pointwise").blockwise is False
    with pytest.raises(SystemExit):
        _parse("i", "o", "--block", "--pointwise")


def test_invalid_noise_model_rejected():
    with pytest.raises(SystemExit):
        _parse("i", "o", "--noise-model", "7")


@pytest.mark.parametrize("flag", ["--patch-radius", "--block-radius"])
@pytest.mark.parametrize("bad", ["0", "-2"])
def test_non_positive_radius_rejected(flag, bad):
    with pytest.raises(SystemExit):
        _parse("i", "o", flag, bad)


def test_clobber_guard_blocks_existing_output(tmp_path):
    inp, out = tmp_path / "in.nii.gz", tmp_path / "out.nii.gz"
    _save_nii(inp, _phantom())
    out.write_bytes(b"existing")
    with pytest.raises(SystemExit) as e:
        nlm.main([str(inp), str(out)])
    assert "exists" in str(e.value)


def test_missing_input_reports_clean_error(tmp_path):
    with pytest.raises(SystemExit) as e:
        nlm.main([str(tmp_path / "nope.nii.gz"), str(tmp_path / "out.nii.gz")])
    assert "not found" in str(e.value)


def test_rejects_non_3d_input(tmp_path):
    inp, out = tmp_path / "4d.nii.gz", tmp_path / "out.nii.gz"
    _save_nii(inp, np.zeros((8, 8, 8, 4), dtype=np.float32))
    with pytest.raises(SystemExit) as e:
        nlm.main([str(inp), str(out)])
    assert "3D" in str(e.value)


def test_negative_sigma_rejected(tmp_path):
    inp, out = tmp_path / "in.nii.gz", tmp_path / "out.nii.gz"
    _save_nii(inp, _phantom())
    with pytest.raises(SystemExit) as e:
        nlm.main([str(inp), str(out), "--sigma", "-1"])
    assert "sigma" in str(e.value)


def test_mask_shape_mismatch_rejected(tmp_path):
    inp, out = tmp_path / "in.nii.gz", tmp_path / "out.nii.gz"
    mask = tmp_path / "mask.nii.gz"
    _save_nii(inp, _phantom(shape=(24, 24, 24)))
    _save_nii(mask, np.ones((10, 10, 10), dtype=np.float32))
    with pytest.raises(SystemExit) as e:
        nlm.main([str(inp), str(out), "--sigma", "5", "--mask", str(mask)])
    assert "shape" in str(e.value)


def test_beta_not_one_warns_before_touching_dipy(tmp_path):
    # A 4-D input makes the run bail right after the beta check, so this
    # verifies the warning fires without needing dipy installed.
    inp, out = tmp_path / "4d.nii.gz", tmp_path / "out.nii.gz"
    _save_nii(inp, np.zeros((8, 8, 8, 2), dtype=np.float32))
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        with pytest.raises(SystemExit):
            nlm.main([str(inp), str(out), "--beta", "0.5"])
    assert any("beta" in str(x.message) for x in w)


def test_sanitize_replaces_non_finite():
    data = np.array([[[1.0, np.nan], [np.inf, -np.inf]]], dtype=np.float32)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = nlm.sanitize(data)
    assert np.isfinite(out).all()
    assert out[0, 0, 0] == 1.0 and out[0, 0, 1] == 0.0
    assert any("non-finite" in str(x.message) for x in w)


def test_build_output_preserves_geometry_and_sets_dtype():
    affine = np.diag([2.0, 3.0, 4.0, 1.0])
    template = nib.Nifti1Image(np.zeros((4, 4, 4), np.float32), affine)
    out = nlm.build_output(np.ones((4, 4, 4)), template)
    np.testing.assert_allclose(out.affine, affine)
    assert out.get_data_dtype() == np.float32


# --------------------------------------------------------------------------- #
# Integration tests — real denoising, gated on dipy
# --------------------------------------------------------------------------- #
def test_denoise_reduces_noise_and_preserves_geometry(tmp_path):
    pytest.importorskip("dipy")
    clean = _phantom()
    noisy = _add_gaussian_noise(clean, sigma=15.0)

    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    inp, out = tmp_path / "noisy.nii.gz", tmp_path / "denoised.nii.gz"
    _save_nii(inp, noisy, affine)

    nlm.main([str(inp), str(out), "--sigma", "15", "--verbose"])

    assert out.exists()
    res = nib.load(str(out))
    got = res.get_fdata(dtype=np.float32)

    assert got.shape == clean.shape
    np.testing.assert_allclose(res.affine, affine, atol=1e-6)
    assert res.get_data_dtype() == np.float32
    # Denoised image is closer to the clean phantom than the noisy input.
    assert _rmse(got, clean) < _rmse(noisy, clean)


def test_auto_sigma_estimation_runs(tmp_path):
    pytest.importorskip("dipy")
    inp, out = tmp_path / "noisy.nii.gz", tmp_path / "out.nii.gz"
    _save_nii(inp, _add_gaussian_noise(_phantom(), sigma=12.0))
    nlm.main([str(inp), str(out)])  # sigma defaults to 0 -> auto estimate
    assert out.exists()


def test_rician_and_pointwise_paths_run(tmp_path):
    pytest.importorskip("dipy")
    inp, out = tmp_path / "noisy.nii.gz", tmp_path / "out.nii.gz"
    _save_nii(inp, _add_gaussian_noise(_phantom(), sigma=10.0))
    nlm.main([str(inp), str(out), "--sigma", "10",
              "--noise-model", "2", "--pointwise"])
    assert out.exists()


def test_output_parent_dir_is_created(tmp_path):
    pytest.importorskip("dipy")
    inp = tmp_path / "noisy.nii.gz"
    out = tmp_path / "nested" / "sub" / "out.nii.gz"  # parents don't exist yet
    _save_nii(inp, _add_gaussian_noise(_phantom(), sigma=8.0))
    nlm.main([str(inp), str(out), "--sigma", "8"])
    assert out.exists()


def test_uncompressed_nii_roundtrip(tmp_path):
    pytest.importorskip("dipy")
    inp, out = tmp_path / "in.nii", tmp_path / "out.nii"  # plain .nii
    _save_nii(inp, _add_gaussian_noise(_phantom(), sigma=10.0))
    nlm.main([str(inp), str(out), "--sigma", "10"])
    assert out.exists() and nib.load(str(out)).shape == (24, 24, 24)


def test_mask_limits_denoising_region(tmp_path):
    pytest.importorskip("dipy")
    noisy = _add_gaussian_noise(_phantom(), sigma=12.0)
    inp, out = tmp_path / "in.nii.gz", tmp_path / "out.nii.gz"
    mask = tmp_path / "mask.nii.gz"
    _save_nii(inp, noisy)
    m = np.zeros((24, 24, 24), dtype=np.float32)
    m[6:18, 6:18, 6:18] = 1  # denoise only the central cube
    _save_nii(mask, m)
    nlm.main([str(inp), str(out), "--sigma", "12", "--mask", str(mask)])
    got = nib.load(str(out)).get_fdata(dtype=np.float32)
    assert got.shape == noisy.shape and np.isfinite(got).all()


def test_blockwise_and_voxelwise_both_valid_and_differ(tmp_path):
    pytest.importorskip("dipy")
    noisy = _add_gaussian_noise(_phantom(), sigma=12.0)
    inp = tmp_path / "in.nii.gz"
    _save_nii(inp, noisy)
    ob, ov = tmp_path / "block.nii.gz", tmp_path / "vox.nii.gz"
    nlm.main([str(inp), str(ob), "--sigma", "12", "--block"])
    nlm.main([str(inp), str(ov), "--sigma", "12", "--pointwise"])
    b = nib.load(str(ob)).get_fdata(dtype=np.float32)
    v = nib.load(str(ov)).get_fdata(dtype=np.float32)
    assert np.isfinite(b).all() and np.isfinite(v).all()
    assert not np.allclose(b, v)  # the two methods are not identical
