"""Tests for functions/ants_denoise.py (ANTs DenoiseImage wrapper).

Unit tests cover CLI parsing, command mapping and guards with an injectable
runner (no ANTs needed). The integration test is gated on a real DenoiseImage
binary (skipped locally; runs in the micapipe container / CI).

Run:  pytest tests/test_ants_denoise.py -v
"""

import importlib.util
import os
import shutil
import types
import warnings

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_MODPATH = os.path.join(_HERE, "..", "functions", "ants_denoise.py")
_spec = importlib.util.spec_from_file_location("ants_denoise", _MODPATH)
ants = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ants)


def _parse(*argv):
    return ants.build_parser().parse_args(list(argv))


def _touch(p, content=b"x"):
    p.write_bytes(content)
    return str(p)


class _FakeRunner:
    """Records the command and returns a chosen exit code."""
    def __init__(self, returncode=0):
        self.returncode = returncode
        self.calls = []

    def __call__(self, cmd, env=None):
        self.calls.append((cmd, env))
        return types.SimpleNamespace(returncode=self.returncode)


# --------------------------------------------------------------------------- #
# Unit tests — no ANTs needed
# --------------------------------------------------------------------------- #
def test_parser_defaults():
    a = _parse("in.nii", "out.nii")
    assert a.patch_radius == 1 and a.block_radius == 2  # ANTs defaults
    assert a.noise_model == 0 and a.sigma == 0.0


def test_command_mapping_gaussian_vs_rician():
    g = ants.build_command(_parse("i.nii", "o.nii", "--noise-model", "0"))
    r = ants.build_command(_parse("i.nii", "o.nii", "--noise-model", "2"))
    assert g[g.index("-n") + 1] == "Gaussian"
    assert r[r.index("-n") + 1] == "Rician"
    assert g[:5] == ["DenoiseImage", "-d", "3", "-i", "i.nii"]


def test_command_maps_radii_and_mask_and_verbose():
    c = ants.build_command(
        _parse("i.nii", "o.nii", "--patch-radius", "2", "--block-radius", "4",
               "--mask", "m.nii", "--verbose"))
    assert c[c.index("-p") + 1] == "2"
    assert c[c.index("-r") + 1] == "4"
    assert c[c.index("-x") + 1] == "m.nii"
    assert "-v" in c


@pytest.mark.parametrize("bad", ["0", "-3"])
def test_non_positive_radius_rejected(bad):
    with pytest.raises(SystemExit):
        _parse("i", "o", "--patch-radius", bad)


def test_missing_input_errors(tmp_path):
    with pytest.raises(SystemExit) as e:
        ants.main([str(tmp_path / "nope.nii.gz"), str(tmp_path / "o.nii.gz")])
    assert "not found" in str(e.value)


def test_clobber_guard(tmp_path):
    inp = _touch(tmp_path / "in.nii.gz")
    out = _touch(tmp_path / "out.nii.gz")
    with pytest.raises(SystemExit) as e:
        ants.main([inp, out])
    assert "exists" in str(e.value)


def test_missing_mask_errors(tmp_path):
    inp = _touch(tmp_path / "in.nii.gz")
    with pytest.raises(SystemExit) as e:
        ants.main([inp, str(tmp_path / "o.nii.gz"), "--mask", str(tmp_path / "no.nii")])
    assert "mask not found" in str(e.value)


def test_missing_binary_errors_cleanly(tmp_path, monkeypatch):
    inp = _touch(tmp_path / "in.nii.gz")
    monkeypatch.setattr(ants.shutil, "which", lambda _: None)
    with pytest.raises(ants.DenoiseError) as e:
        ants.denoise_file(ants.build_parser().parse_args([inp, str(tmp_path / "o.nii.gz")]))
    assert "not on PATH" in str(e.value)


def test_successful_run_builds_correct_command(tmp_path, monkeypatch):
    inp = _touch(tmp_path / "in.nii.gz")
    out = tmp_path / "nested" / "out.nii.gz"  # parent doesn't exist yet
    monkeypatch.setattr(ants.shutil, "which", lambda _: "/usr/bin/DenoiseImage")
    fake = _FakeRunner(returncode=0)
    args = ants.build_parser().parse_args([inp, str(out), "--noise-model", "2",
                                           "--threads", "4"])
    ants.denoise_file(args, runner=fake)
    assert out.parent.is_dir()                      # parent auto-created
    cmd, env = fake.calls[0]
    assert "Rician" in cmd
    assert env["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] == "4"


def test_nonzero_exit_is_reported(tmp_path, monkeypatch):
    inp = _touch(tmp_path / "in.nii.gz")
    monkeypatch.setattr(ants.shutil, "which", lambda _: "/usr/bin/DenoiseImage")
    args = ants.build_parser().parse_args([inp, str(tmp_path / "o.nii.gz")])
    with pytest.raises(ants.DenoiseError):
        ants.denoise_file(args, runner=_FakeRunner(returncode=1))


def test_sigma_warns(tmp_path, monkeypatch):
    inp = _touch(tmp_path / "in.nii.gz")
    monkeypatch.setattr(ants.shutil, "which", lambda _: "/usr/bin/DenoiseImage")
    args = ants.build_parser().parse_args([inp, str(tmp_path / "o.nii.gz"), "--sigma", "5"])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ants.denoise_file(args, runner=_FakeRunner(0))
    assert any("sigma" in str(x.message) for x in w)


# --------------------------------------------------------------------------- #
# Integration — real DenoiseImage, gated
# --------------------------------------------------------------------------- #
def test_real_denoiseimage_roundtrip(tmp_path):
    if shutil.which("DenoiseImage") is None:
        pytest.skip("DenoiseImage not installed (runs in micapipe container/CI)")
    import numpy as np
    import nibabel as nib
    rng = np.random.default_rng(0)
    vol = np.zeros((24, 24, 24), np.float32); vol[6:18, 6:18, 6:18] = 100
    noisy = (vol + rng.normal(0, 12, vol.shape)).astype(np.float32)
    inp, out = tmp_path / "in.nii.gz", tmp_path / "out.nii.gz"
    nib.save(nib.Nifti1Image(noisy, np.eye(4)), str(inp))
    ants.main([str(inp), str(out), "--noise-model", "2"])
    assert out.exists() and nib.load(str(out)).shape == (24, 24, 24)
