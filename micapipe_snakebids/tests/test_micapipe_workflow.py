"""Tests for micapipe_snakebids/workflow/lib/micapipe_workflow.py.

Run with pytest, or directly: python micapipe_snakebids/tests/test_micapipe_workflow.py
"""

import json
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "workflow", "lib"))

import micapipe_workflow as mw  # noqa: E402

# Every option at a non-default value, the way the CLI would set them.
EVERYTHING = {
    "t1w_str": ["run-1_T1w"], "uni": True, "mf": 3.0,
    "freesurfer": True, "alt_fs": True, "alt_fs_tta": True,
    "surf_dir": "/data/out/fs/sub-01_ses-01", "surf_t1": "/data/bids/t1.nii.gz",
    "atlas": ["schaefer-100", "economo"],
    "flair_scan_str": "FLAIR",
    "dwi_main": ["acq-b700_dir-AP_dwi", "acq-b2000_dir-AP_dwi"], "dwi_rpe": ["dir-PA_dwi"],
    "dwi_processed": "/data/out/dwi.mif", "rpe_all": True, "b0thr": 61,
    "no_bvalue_scaling": True, "dwi_upsample": True, "dwi_acq": "b2000",
    "tracts": "10M", "keep_tck": True, "auto_tract": True,
    "weighted_sc": "/data/w.nii.gz", "tck": "/data/t.tck",
    "microstructural_img": "T1map", "mpc_acq": "T1map",
    "asl_scan_str": "asl", "m0_scan_str": "m0scan",
    "func_pe": "task-rest_acq-APse_bold", "func_rpe": "/data/bids/rpe.nii.gz",
    "main_scan_run": 1, "phase_reversal_run": 1,
    "topup_config": "/data/topup.cnf", "icafix_training": "/data/fix.RData",
    "ses_anat": "01", "nsr": True, "gsr": True, "no_fix": True, "drop_tr": True, "no_fc": True,
    "nocleanup": True,
}


def render(stage, config, scan=None, **kw):
    return mw.micapipe_arguments(
        stage, config, subject="01", session="01", bids_root="/bids", out_root="/out",
        threads=4, licence="/opt/licence.txt", tmp_root="/tmp", scan=scan, **kw,
    )


def test_modules_are_ordered_and_validated():
    assert mw.requested_modules(["sc", "structural"]) == ("structural", "sc")
    assert mw.requested_modules("dwi,func") == ("dwi", "func")
    try:
        mw.requested_modules(["bogus"])
    except mw.WorkflowConfigError as err:
        assert "bogus" in str(err)
    else:
        raise AssertionError("unknown module accepted")


def test_closure_pulls_in_dependencies_in_order():
    assert mw.module_closure(["sc"]) == ("structural", "surface", "post_structural", "dwi", "sc")
    assert mw.module_closure(["qc_subj"]) == ("structural", "qc_subj")
    assert mw.module_closure(["func"], skip=("structural", "surface", "post_structural")) == ("func",)
    chain = mw.module_closure(list(mw.MODULE_ORDER))
    assert chain == mw.MODULE_ORDER
    for stage in chain:
        for dep in mw.MODULE_DEPENDENCIES[stage]:
            assert chain.index(dep) < chain.index(stage)
    assert mw.predecessor("sc", chain) == "dwi"
    assert mw.predecessor("structural", chain) is None


def test_module_names_follow_micapipe():
    assert mw.module_name("surface", {}) == "proc_surf-fastsurfer"
    assert mw.module_name("surface", {"freesurfer": "TRUE"}) == "proc_surf-freesurfer"
    assert mw.module_name("surface", {"freesurfer": "FALSE"}) == "proc_surf-fastsurfer"
    assert mw.module_name("surface", {"alt_fs": True}) == "proc_surf-fastsurfer"
    assert mw.module_name("dwi", {}) == "proc_dwi"
    assert mw.module_name("dwi", {"dwi_acq": "acq-b2000"}) == "proc_dwi_acq-b2000"
    assert mw.module_name("sc", {}) == "SC-40M"
    assert mw.module_name("sc", {"tracts": "10M", "dwi_acq": "b2000"}) == "SC-10M_acq-b2000"
    assert mw.module_name("mpc", {}) == "MPC-qMRI"
    assert mw.module_name("mpc_swm", {"mpc_acq": "T1map"}) == "MPC-SWM-T1map"
    assert mw.module_name("func", {}, "task-rest_acq-AP_bold") == "proc_func-desc-se_task-rest_acq-AP_bold"
    assert (
        mw.module_name("func", {}, "task-x_echo-1_bold,task-x_echo-2_bold")
        == "proc_func-desc-me_task-x_bold"
    )
    assert mw.module_record_name("gd", {}) == "{stem}_module-GD.json"


def test_defaults_render_no_stray_words():
    """The shipped rules rendered "FALSE uni", "-freesurfer FALSE", ... and hung micapipe."""
    old_style = {
        "uni": "FALSE", "freesurfer": "FALSE", "surf_dir": "FALSE", "rpe_all": "FALSE",
        "keep_tck": "FALSE", "auto_tract": "FALSE", "weighted_sc": "FALSE", "tck": "FALSE",
        "dwi_main": "DEFAULT", "dwi_rpe": "DEFAULT", "dwi_processed": "FALSE",
        "t1w_str": "DEFAULT", "atlas": "DEFAULT", "mpc_acq": "FALSE",
    }
    for config in ({}, old_style):
        for stage in mw.MODULE_ORDER:
            args = render(stage, config, scan="task-rest_acq-AP_bold")
            mw.check_arguments(args)
            assert "FALSE" not in args and "DEFAULT" not in args, (stage, args)
    assert render("structural", {})[-1] == "-proc_structural"


def test_every_option_renders_as_micapipe_expects():
    seen = set()
    for stage in mw.MODULE_ORDER:
        args = render(stage, EVERYTHING, scan="task-rest_acq-AP_bold")
        mw.check_arguments(args)
        seen.update(a for a in args if a.startswith("-"))
    expected = {flag for _, flag, _ in mw.PARAMETERS.values()}
    expected |= {flag for _, flag, _ in mw.SHARED_PARAMETERS.values()}
    assert expected <= seen, sorted(expected - seen)

    structural = render("structural", EVERYTHING)
    assert structural[structural.index("-T1wStr") + 1] == "run-1_T1w"
    assert "-uni" in structural and structural[structural.index("-mf") + 1] == "3.0"
    surface = render("surface", EVERYTHING)
    assert {"-freesurfer", "-alt_fs", "-alt_fs_tta"} <= set(surface)
    by_name = render("surface", {"surf_t1": "T1w"})
    assert by_name[by_name.index("-T1") + 1] == "/bids/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii.gz"
    post = render("post_structural", EVERYTHING)
    assert post[post.index("-atlas") + 1] == "schaefer-100,economo"
    dwi = render("dwi", EVERYTHING)
    assert dwi[dwi.index("-dwi_main") + 1] == (
        "/bids/sub-01/ses-01/dwi/sub-01_ses-01_acq-b700_dir-AP_dwi.nii.gz,"
        "/bids/sub-01/ses-01/dwi/sub-01_ses-01_acq-b2000_dir-AP_dwi.nii.gz"
    )
    assert dwi[dwi.index("-dwi_acq") + 1] == "b2000"
    sc = render("sc", EVERYTHING)
    assert {"-keep_tck", "-autoTract"} <= set(sc) and "-dwi_acq" in sc
    mpc_swm = render("mpc_swm", EVERYTHING)
    assert mpc_swm[mpc_swm.index("-microstructural_img") + 1] == (
        "/bids/sub-01/ses-01/anat/sub-01_ses-01_T1map.nii.gz"
    )
    func = render("func", EVERYTHING, scan="task-rest_acq-AP_bold")
    assert func[func.index("-mainScanStr") + 1] == "task-rest_acq-AP_bold"
    assert func[func.index("-func_pe") + 1].endswith("/func/sub-01_ses-01_task-rest_acq-APse_bold.nii.gz")
    assert func[func.index("-func_rpe") + 1] == "/data/bids/rpe.nii.gz"
    # Options of other modules never leak into a call.
    assert "-tracts" not in render("dwi", EVERYTHING)
    assert "-atlas" not in render("structural", EVERYTHING)


def test_paths_are_translated_into_the_container():
    to_container = mw.PathMap([("/data/bids", "/bids"), ("/data/out", "/out")])
    surface = render("surface", EVERYTHING, to_container=to_container)
    assert surface[surface.index("-surf_dir") + 1] == "/out/fs/sub-01_ses-01"
    assert surface[surface.index("-T1") + 1] == "/bids/t1.nii.gz"
    assert to_container("/elsewhere/x") == "/elsewhere/x"
    assert to_container("/data/bids") == "/bids"


def test_check_arguments_catches_what_hangs_micapipe():
    for bad in (
        ["-sub", "01", "-proc_structural", "-mf", "3", "FALSE", "uni"],
        ["-sub", "01", "-proc_surf", "-freesurfer", "FALSE"],
        ["-sub", "01", "-SC", "-keep_tck", "FALSE"],
    ):
        try:
            mw.check_arguments(bad)
        except mw.WorkflowConfigError as err:
            assert "loop forever" in str(err)
        else:
            raise AssertionError(bad)


def _touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").close()


def test_session_discovery():
    with tempfile.TemporaryDirectory() as bids:
        _touch(f"{bids}/sub-HC001/ses-01/anat/sub-HC001_ses-01_T1w.nii.gz")  # no run entity
        _touch(f"{bids}/sub-HC001/ses-02/anat/sub-HC001_ses-02_run-1_T1w.nii.gz")
        _touch(f"{bids}/sub-PX001/ses-a1/anat/sub-PX001_ses-a1_UNIT1.nii.gz")
        _touch(f"{bids}/sub-PX002/ses-01/dwi/sub-PX002_ses-01_dwi.nii.gz")  # no anatomy
        assert mw.discover_sessions(bids) == [("HC001", "01"), ("HC001", "02"), ("PX001", "a1")]
        assert mw.discover_sessions(bids, participant_label=["sub-HC001"], session_label="02") == [
            ("HC001", "02")
        ]
        assert mw.discover_sessions(bids, exclude_participant_label=["HC001"]) == [("PX001", "a1")]
        try:
            mw.discover_sessions(bids, participant_label=["NOPE"])
        except mw.WorkflowConfigError as err:
            assert "NOPE" in str(err)
        else:
            raise AssertionError("empty selection accepted")
    with tempfile.TemporaryDirectory() as bids:
        _touch(f"{bids}/sub-01/anat/sub-01_T1w.nii.gz")
        try:
            mw.discover_sessions(bids)
        except mw.WorkflowConfigError as err:
            assert "ses-" in str(err)
        else:
            raise AssertionError("session-less dataset accepted")


def test_default_func_scan_matches_micapipe():
    with tempfile.TemporaryDirectory() as bids:
        func = f"{bids}/sub-01/ses-01/func"
        _touch(f"{func}/sub-01_ses-01_task-rest_acq-AP_bold.nii.gz")
        _touch(f"{func}/sub-01_ses-01_task-rest_acq-APse_bold.nii.gz")  # not *_acq-AP_*
        assert mw.default_func_scans(bids, "01", "01") == ["task-rest_acq-AP_bold"]
        _touch(f"{func}/sub-01_ses-01_task-rest_acq-AP_run-2_bold.nii.gz")
        try:
            mw.default_func_scans(bids, "01", "01")
        except mw.WorkflowConfigError as err:
            assert "--func-scans" in str(err)
        else:
            raise AssertionError("ambiguous rest scans accepted")
        assert mw.default_func_scans(bids, "01", "01", main_scan_run=2) == ["task-rest_acq-AP_run-2_bold"]


def test_module_record_check():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "r.json")
        record = {"Module": "proc_structural", "Subject": "sub-01", "Session": "01", "Status": "COMPLETED"}
        with open(path, "w") as stream:
            json.dump(record, stream)
        mw.check_module_record(path, "proc_structural", "01", "01")
        record["Status"] = "INCOMPLETE"
        with open(path, "w") as stream:
            json.dump(record, stream)
        for args in (("proc_structural", "01", "01"), ("proc_dwi", "01", "01")):
            try:
                mw.check_module_record(path, *args)
            except RuntimeError as err:
                assert "did not complete" in str(err)
            else:
                raise AssertionError(args)
        try:
            mw.check_module_record(os.path.join(tmp, "missing.json"), "GD", "01", "01")
        except RuntimeError as err:
            assert "no record" in str(err)
        else:
            raise AssertionError("missing record accepted")


if __name__ == "__main__":
    failures = 0
    for name, test in sorted(globals().items()):
        if name.startswith("test_") and callable(test):
            try:
                test()
                print(f"PASS {name}")
            except Exception as err:  # noqa: BLE001
                failures += 1
                print(f"FAIL {name}: {err!r}")
    sys.exit(1 if failures else 0)
