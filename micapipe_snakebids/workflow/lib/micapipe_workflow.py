"""Pure helpers for the micapipe Snakemake workflow.

Nothing here imports Snakemake, so it can be unit tested on its own
(``micapipe_snakebids/tests/test_micapipe_workflow.py``).

The workflow drives one ``micapipe`` call per module and session, and it takes
micapipe's own record as the proof that a module finished: the
``QC/<sub>_<ses>_module-<name>.json`` file every module writes, with
``"Status": "COMPLETED"``. micapipe runs each module as ``$COMMAND | tee``, so
its exit status cannot be trusted on its own.
"""

from __future__ import annotations

import glob
import json
import os
import re
import shlex
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

OUTPUT_DIRNAME = "micapipe_v0.2.0"

#: Every module, in the order micapipe needs them. A session's modules always
#: run one at a time in this order: micapipe rewrites session-level files
#: (QC records, the processed-subjects ledger) and two modules running at once
#: on one session can corrupt them.
MODULE_ORDER = (
    "structural",
    "surface",
    "post_structural",
    "flair",
    "dwi",
    "sc",
    "mpc",
    "mpc_swm",
    "gd",
    "asl",
    "func",
    "swm",
    "qc_subj",
)

#: The micapipe flag that runs each module.
MODULE_FLAGS = {
    "structural": "-proc_structural",
    "surface": "-proc_surf",
    "post_structural": "-post_structural",
    "flair": "-proc_flair",
    "dwi": "-proc_dwi",
    "sc": "-SC",
    "mpc": "-MPC",
    "mpc_swm": "-MPC_SWM",
    "gd": "-GD",
    "asl": "-proc_asl",
    "func": "-proc_func",
    "swm": "-SWM",
    "qc_subj": "-QC_subj",
}

#: What each module refuses to start without (micapipe_check_dependency calls
#: in functions/*.sh, plus the inputs each script reads). A requested module
#: pulls these in, so ``--modules sc`` runs the structural chain and DWI first,
#: and a module whose record already says COMPLETED is simply not re-run.
MODULE_DEPENDENCIES = {
    "structural": (),
    "surface": ("structural",),
    "post_structural": ("surface",),
    "flair": ("post_structural",),
    "dwi": ("post_structural",),
    "sc": ("dwi",),
    "mpc": ("post_structural",),
    "mpc_swm": ("post_structural",),
    "gd": ("post_structural",),
    "asl": ("post_structural",),
    "func": ("post_structural",),
    "swm": ("post_structural",),
    "qc_subj": ("structural",),
}

#: Modules that need a FreeSurfer licence.
LICENCE_MODULES = frozenset({"structural", "surface", "post_structural"})

# --------------------------------------------------------------------------
#  Parameters: config key -> (module, micapipe flag, kind)
#
#  kind is one of
#    flag  : true emits the bare flag, false/None omits it
#    value : a scalar, emitted as ``-flag value`` when set
#    list  : a list (or comma-separated string), emitted as ``-flag a,b``
#    path  : a path, translated into the container when one is used
#    paths : a list of paths, comma-joined after translation
#    scan  : a BIDS file of the session, given either as a path or as the part
#            of its name after ``sub-X_ses-Y_`` (``acq-b2000_dir-AP_dwi``);
#            resolved per session
#    scans : a list of ``scan`` values, comma-joined
#
#  Only keys whose module is part of a run are rendered, and a key left unset
#  is never rendered, so micapipe keeps its own default. Nothing is ever
#  rendered as the string "FALSE": micapipe's parser treats a bare word it
#  does not expect as neither an option nor a value and loops on it forever.
# --------------------------------------------------------------------------
PARAMETERS: Mapping[str, tuple[str, str, str]] = {
    # proc_structural
    "t1w_str": ("structural", "-T1wStr", "list"),
    "uni": ("structural", "-uni", "flag"),
    "mf": ("structural", "-mf", "value"),
    # proc_surf
    "freesurfer": ("surface", "-freesurfer", "flag"),
    "alt_fs": ("surface", "-alt_fs", "flag"),
    "alt_fs_tta": ("surface", "-alt_fs_tta", "flag"),
    "surf_dir": ("surface", "-surf_dir", "path"),
    "surf_t1": ("surface", "-T1", "scan"),
    # post_structural
    "atlas": ("post_structural", "-atlas", "list"),
    # proc_flair
    "flair_scan_str": ("flair", "-flairScanStr", "value"),
    # proc_dwi
    "dwi_main": ("dwi", "-dwi_main", "scans"),
    "dwi_rpe": ("dwi", "-dwi_rpe", "scans"),
    "dwi_processed": ("dwi", "-dwi_processed", "path"),
    "rpe_all": ("dwi", "-rpe_all", "flag"),
    "b0thr": ("dwi", "-b0thr", "value"),
    "no_bvalue_scaling": ("dwi", "-no_bvalue_scaling", "flag"),
    "dwi_upsample": ("dwi", "-dwi_upsample", "flag"),
    # SC
    "tracts": ("sc", "-tracts", "value"),
    "keep_tck": ("sc", "-keep_tck", "flag"),
    "auto_tract": ("sc", "-autoTract", "flag"),
    "weighted_sc": ("sc", "-weighted_SC", "path"),
    "tck": ("sc", "-tck", "path"),
    # proc_asl
    "asl_scan_str": ("asl", "-aslScanStr", "value"),
    "m0_scan_str": ("asl", "-m0ScanStr", "value"),
    # proc_func (-mainScanStr comes from func_scans, one job per scan)
    "func_pe": ("func", "-func_pe", "scan"),
    "func_rpe": ("func", "-func_rpe", "scan"),
    "main_scan_run": ("func", "-mainScanRun", "value"),
    "phase_reversal_run": ("func", "-phaseReversalRun", "value"),
    "topup_config": ("func", "-topupConfig", "path"),
    "icafix_training": ("func", "-icafixTraining", "path"),
    "ses_anat": ("func", "-sesAnat", "value"),
    "nsr": ("func", "-NSR", "flag"),
    "gsr": ("func", "-GSR", "flag"),
    "no_fix": ("func", "-noFIX", "flag"),
    "drop_tr": ("func", "-dropTR", "flag"),
    "no_fc": ("func", "-noFC", "flag"),
}

#: Parameters more than one module takes. dwi_acq names which DWI acquisition
#: proc_dwi processed and SC reads; mpc_acq names the MPC acquisition for both
#: MPC and MPC_SWM, and MPC_SWM takes the same microstructural image.
SHARED_PARAMETERS: Mapping[str, tuple[tuple[str, ...], str, str]] = {
    "dwi_acq": (("dwi", "sc"), "-dwi_acq", "value"),
    "mpc_acq": (("mpc", "mpc_swm"), "-mpc_acq", "value"),
    "microstructural_img": (("mpc", "mpc_swm"), "-microstructural_img", "scan"),
}

#: The BIDS datatype directory a ``scan`` value of each module lives in.
SCAN_DATATYPE = {"surface": "anat", "dwi": "dwi", "func": "func", "mpc": "anat", "mpc_swm": "anat"}


class WorkflowConfigError(ValueError):
    """A configuration the workflow refuses before scheduling anything."""


# --------------------------------------------------------------------------
#  Values from config and the CLI
# --------------------------------------------------------------------------

_TRUE = {"true", "yes", "on", "1"}
_FALSE = {"false", "no", "off", "0", "none", "default", ""}


def as_bool(value: Any, key: str = "value") -> bool:
    """Read a config boolean, accepting the old "TRUE"/"FALSE" strings."""
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return False
    raise WorkflowConfigError(f"{key} must be true or false, got {value!r}")


def is_unset(value: Any) -> bool:
    """True for values that mean "let micapipe use its default"."""
    if value is None or value is False:
        return True
    if isinstance(value, str) and value.strip().upper() in {"", "DEFAULT", "FALSE"}:
        return True
    if isinstance(value, (list, tuple)) and not value:
        return True
    return False


def as_list(value: Any) -> list[str]:
    """A list from a YAML list, a CLI list, or a comma-separated string."""
    if is_unset(value):
        return []
    if isinstance(value, (list, tuple)):
        items: list[str] = []
        for item in value:
            items.extend(as_list(item))
        return items
    return [part.strip() for part in str(value).split(",") if part.strip()]


def strip_prefix(value: Any, prefix: str) -> str:
    text = str(value)
    return text[len(prefix):] if text.startswith(prefix) else text


# --------------------------------------------------------------------------
#  Modules
# --------------------------------------------------------------------------


def requested_modules(value: Any) -> tuple[str, ...]:
    """The modules asked for, validated and in canonical order."""
    names = [name.strip().lower() for name in as_list(value)]
    if not names:
        raise WorkflowConfigError("no micapipe modules were requested (--modules)")
    unknown = sorted(set(names) - set(MODULE_ORDER))
    if unknown:
        raise WorkflowConfigError(
            f"unknown modules {unknown}; choose from {list(MODULE_ORDER)}"
        )
    return tuple(name for name in MODULE_ORDER if name in set(names))


def module_closure(modules: Iterable[str], skip: Iterable[str] = ()) -> tuple[str, ...]:
    """``modules`` plus everything they depend on, in canonical order.

    ``skip`` names modules this session does not run itself (the structural
    chain of a func-only session that borrows another session's anatomy).
    """
    skipped = set(skip)
    wanted: set[str] = set()
    pending = [name for name in modules if name not in skipped]
    while pending:
        name = pending.pop()
        if name in wanted:
            continue
        wanted.add(name)
        pending.extend(dep for dep in MODULE_DEPENDENCIES[name] if dep not in skipped)
    return tuple(name for name in MODULE_ORDER if name in wanted)


def predecessor(module: str, chain: Sequence[str]) -> str | None:
    """The module that runs just before ``module`` in a session's chain."""
    if module not in chain:
        return None
    index = list(chain).index(module)
    return chain[index - 1] if index else None


# --------------------------------------------------------------------------
#  Module records (QC/<stem>_module-<name>.json)
# --------------------------------------------------------------------------


def session_root(out_dir: str, subject: str, session: str) -> str:
    return f"{out_dir}/{OUTPUT_DIRNAME}/sub-{subject}/ses-{session}"


def surface_recon(config: Mapping[str, Any]) -> str:
    """FastSurfer unless -freesurfer; -alt_fs keeps FastSurfer's name."""
    return "freesurfer" if as_bool(config.get("freesurfer"), "freesurfer") else "fastsurfer"


def acq_suffix(value: Any) -> str:
    """micapipe's ``_acq-<label>`` suffix for -dwi_acq, or ''."""
    if is_unset(value):
        return ""
    return "_acq-" + strip_prefix(str(value).strip(), "acq-")


def module_name(module: str, config: Mapping[str, Any], scan: str | None = None) -> str:
    """The ``Module`` micapipe records for ``module`` under ``config``.

    Mirrors functions/*.sh: proc_surf-<recon>, proc_dwi[_acq-X],
    SC-<tracts>[_acq-X], MPC-<acq|qMRI>, MPC-SWM-<acq|qMRI>, and
    proc_func-desc-<se|me>_<scan>.
    """
    if module == "structural":
        return "proc_structural"
    if module == "surface":
        return f"proc_surf-{surface_recon(config)}"
    if module == "post_structural":
        return "post_structural"
    if module == "flair":
        return "proc_flair"
    if module == "dwi":
        return "proc_dwi" + acq_suffix(config.get("dwi_acq"))
    if module == "sc":
        tracts = "40M" if is_unset(config.get("tracts")) else str(config["tracts"])
        return f"SC-{tracts}" + acq_suffix(config.get("dwi_acq"))
    if module in ("mpc", "mpc_swm"):
        acq = "qMRI" if is_unset(config.get("mpc_acq")) else str(config["mpc_acq"])
        return ("MPC-" if module == "mpc" else "MPC-SWM-") + acq
    if module == "gd":
        return "GD"
    if module == "asl":
        return "proc_asl"
    if module == "func":
        if not scan:
            raise ValueError("func needs a scan")
        return "proc_func-" + func_tag(scan)
    if module == "swm":
        return "SWM"
    if module == "qc_subj":
        return "QC_subj"
    raise ValueError(f"unknown module {module!r}")


def module_record_name(module: str, config: Mapping[str, Any], scan: str | None = None) -> str:
    """The file name of a module's record (QC_subj's is its structural report)."""
    if module == "qc_subj":
        return "{stem}_module-proc_structural_qc-report.pdf"
    return "{stem}_module-" + module_name(module, config, scan) + ".json"


def func_tag(scan: str) -> str:
    """``desc-se_<scan>`` for one image, ``desc-me_<first echo>`` for several.

    02_proc-func.sh names a scan by its first image with ``echo-1_`` removed,
    and calls a run of several images multi-echo.
    """
    images = [name.strip() for name in scan.split(",") if name.strip()]
    first = images[0].replace("echo-1_", "")
    desc = "se" if len(images) == 1 else "me"
    return f"desc-{desc}_{first}"


def check_module_record(path: str, expected_module: str, subject: str, session: str) -> None:
    """Raise unless micapipe recorded ``expected_module`` as COMPLETED."""
    if not os.path.isfile(path):
        raise RuntimeError(f"micapipe wrote no record at {path}")
    with open(path) as stream:
        text = stream.read()
    try:
        record = json.loads(text)
    except json.JSONDecodeError as err:
        raise RuntimeError(f"{path} is not valid JSON: {err}") from err
    problems = []
    if record.get("Status") != "COMPLETED":
        problems.append(f"Status is {record.get('Status')!r} (progress {record.get('Progress')})")
    if record.get("Module") != expected_module:
        problems.append(f"Module is {record.get('Module')!r}, expected {expected_module!r}")
    if strip_prefix(record.get("Subject", ""), "sub-") != subject:
        problems.append(f"Subject is {record.get('Subject')!r}, expected {subject!r}")
    if strip_prefix(record.get("Session", ""), "ses-") != session:
        problems.append(f"Session is {record.get('Session')!r}, expected {session!r}")
    if problems:
        raise RuntimeError(f"micapipe did not complete {expected_module}: " + "; ".join(problems))


# --------------------------------------------------------------------------
#  Sessions
# --------------------------------------------------------------------------


def discover_sessions(
    bids_dir: str,
    participant_label: Any = None,
    exclude_participant_label: Any = None,
    session_label: Any = None,
) -> list[tuple[str, str]]:
    """(subject, session) pairs that hold a T1w-like anatomical image.

    Any image micapipe's proc_structural can start from counts (T1w, or an
    MP2RAGE UNIT1 for -uni), with or without run/acq entities: micapipe
    picks the image itself through -T1wStr.
    """
    wanted = {strip_prefix(s, "sub-") for s in as_list(participant_label)}
    excluded = {strip_prefix(s, "sub-") for s in as_list(exclude_participant_label)}
    sessions_wanted = {strip_prefix(s, "ses-") for s in as_list(session_label)}

    pairs: list[tuple[str, str]] = []
    sessionless: list[str] = []
    for subject_dir in sorted(glob.glob(os.path.join(bids_dir, "sub-*"))):
        if not os.path.isdir(subject_dir):
            continue
        subject = strip_prefix(os.path.basename(subject_dir), "sub-")
        if (wanted and subject not in wanted) or subject in excluded:
            continue
        session_dirs = sorted(glob.glob(os.path.join(subject_dir, "ses-*")))
        if not session_dirs and _has_anatomy(os.path.join(subject_dir, "anat")):
            sessionless.append(subject)
        for session_dir in session_dirs:
            session = strip_prefix(os.path.basename(session_dir), "ses-")
            if sessions_wanted and session not in sessions_wanted:
                continue
            if _has_anatomy(os.path.join(session_dir, "anat")):
                pairs.append((subject, session))

    if sessionless and not pairs:
        raise WorkflowConfigError(
            "these subjects have no ses-* directories, which this workflow does not "
            "support yet: " + ", ".join(f"sub-{s}" for s in sessionless[:10])
        )
    if not pairs:
        selection = []
        if wanted:
            selection.append(f"participant_label={sorted(wanted)}")
        if excluded:
            selection.append(f"exclude_participant_label={sorted(excluded)}")
        if sessions_wanted:
            selection.append(f"session_label={sorted(sessions_wanted)}")
        raise WorkflowConfigError(
            f"no subject/session with a T1w (or UNIT1) image under {bids_dir}"
            + (" for " + ", ".join(selection) if selection else "")
        )
    return pairs


def _has_anatomy(anat_dir: str) -> bool:
    return any(
        glob.glob(os.path.join(anat_dir, f"*_{suffix}.nii*")) for suffix in ("T1w", "UNIT1")
    )


def bids_session_dir(bids_dir: str, subject: str, session: str) -> str:
    return f"{bids_dir}/sub-{subject}/ses-{session}"


def default_func_scans(
    bids_dir: str, subject: str, session: str, main_scan_run: Any = None
) -> list[str]:
    """The scan micapipe picks without -mainScanStr, by its name.

    utilities.sh looks for ``<sub>_<ses>_task-rest_acq-AP_*bold.nii*``; with
    several matches it needs -mainScanRun (1-based, in ``ls`` order).
    """
    prefix = f"sub-{subject}_ses-{session}_"
    pattern = os.path.join(
        bids_session_dir(bids_dir, subject, session), "func", prefix + "task-rest_acq-AP_*bold.nii*"
    )
    found = sorted(glob.glob(pattern))
    scans = [_scan_of(path, prefix) for path in found]
    if len(scans) <= 1:
        return scans
    if is_unset(main_scan_run):
        raise WorkflowConfigError(
            f"sub-{subject} ses-{session} has {len(scans)} rest scans ({', '.join(scans)}); "
            "name the ones to process with --func-scans, or pick one with --main-scan-run"
        )
    index = int(main_scan_run)
    return [scans[index - 1] if 1 <= index <= len(scans) else scans[0]]


def _scan_of(path: str, prefix: str) -> str:
    name = os.path.basename(path)
    name = re.sub(r"\.nii(\.gz)?$", "", name)
    return name[len(prefix):] if name.startswith(prefix) else name


def resolve_scan(value: str, bids_root: str, subject: str, session: str, datatype: str) -> str:
    """A ``scan`` value as a path: kept if it is one, else built in the session."""
    if "/" in value:
        return value
    name = value if value.startswith(f"sub-{subject}_") else f"sub-{subject}_ses-{session}_{value}"
    if not re.search(r"\.nii(\.gz)?$", name):
        name += ".nii.gz"
    return f"{bids_root}/sub-{subject}/ses-{session}/{datatype}/{name}"


# --------------------------------------------------------------------------
#  Rendering one micapipe call
# --------------------------------------------------------------------------


class PathMap:
    """Translate host paths into the container's view of them."""

    def __init__(self, pairs: Sequence[tuple[str, str]] = ()):
        # Longest host prefix first, so /data/bids/sub wins over /data/bids.
        self.pairs = sorted(
            ((os.path.normpath(host), inner) for host, inner in pairs),
            key=lambda pair: len(pair[0]),
            reverse=True,
        )

    def __call__(self, path: str) -> str:
        norm = os.path.normpath(str(path))
        for host, inner in self.pairs:
            if norm == host:
                return inner
            if norm.startswith(host + os.sep):
                return inner + norm[len(host):]
        return str(path)


def render_parameters(
    module: str,
    config: Mapping[str, Any],
    *,
    subject: str,
    session: str,
    bids_root: str,
    to_container=lambda path: path,
) -> list[str]:
    """The micapipe arguments ``config`` asks for, for one module call."""
    args: list[str] = []
    entries = [(key, (mod,), flag, kind) for key, (mod, flag, kind) in PARAMETERS.items()]
    entries += [(key, mods, flag, kind) for key, (mods, flag, kind) in SHARED_PARAMETERS.items()]
    seen: set[tuple[str, str]] = set()
    for key, modules, flag, kind in entries:
        if module not in modules or (key, flag) in seen:
            continue
        seen.add((key, flag))
        value = config.get(key)
        if kind == "flag":
            if as_bool(value, key):
                args.append(flag)
            continue
        if is_unset(value):
            continue
        if kind == "value":
            args += [flag, str(value)]
        elif kind == "list":
            args += [flag, ",".join(as_list(value))]
        elif kind == "path":
            args += [flag, to_container(str(value))]
        elif kind == "paths":
            args += [flag, ",".join(to_container(v) for v in as_list(value))]
        elif kind in ("scan", "scans"):
            datatype = SCAN_DATATYPE[module]
            values = as_list(value) if kind == "scans" else [str(value)]
            paths = [
                to_container(resolve_scan(v, bids_root, subject, session, datatype))
                if "/" in v
                else resolve_scan(v, bids_root, subject, session, datatype)
                for v in values
            ]
            args += [flag, ",".join(paths)]
        else:  # pragma: no cover - table is static
            raise ValueError(kind)
    return args


def micapipe_arguments(
    module: str,
    config: Mapping[str, Any],
    *,
    subject: str,
    session: str,
    bids_root: str,
    out_root: str,
    threads: int,
    licence: str | None,
    tmp_root: str | None,
    to_container=lambda path: path,
    scan: str | None = None,
) -> list[str]:
    """The whole micapipe argument list for one module of one session."""
    args = [
        "-sub", subject,
        "-ses", session,
        "-bids", bids_root,
        "-out", out_root,
        "-threads", str(threads),
    ]
    if licence:
        args += ["-fs_licence", licence]
    if tmp_root:
        args += ["-tmpDir", tmp_root]
    if as_bool(config.get("nocleanup"), "nocleanup"):
        args.append("-nocleanup")
    args.append(MODULE_FLAGS[module])
    if module == "func":
        args += ["-mainScanStr", scan]
    args += render_parameters(
        module, config, subject=subject, session=session, bids_root=bids_root,
        to_container=to_container,
    )
    return args


def check_arguments(args: Sequence[str]) -> None:
    """Refuse any argument micapipe's parser would spin on.

    micapipe's ``while [ $# -gt 0 ]; case`` loop only shifts what it
    recognises: a word that is neither one of its options nor the value of the
    option before it is never consumed, and micapipe loops on it at 100% CPU.
    """
    value_options = {
        "-sub", "-ses", "-bids", "-out", "-threads", "-fs_licence", "-tmpDir",
        "-mainScanStr",
    }
    for key, (module, flag, kind) in PARAMETERS.items():
        if kind != "flag":
            value_options.add(flag)
    for key, (modules, flag, kind) in SHARED_PARAMETERS.items():
        if kind != "flag":
            value_options.add(flag)
    expecting_value = False
    for arg in args:
        if expecting_value:
            if arg == "" or arg is None:
                raise WorkflowConfigError("an option was rendered with an empty value")
            expecting_value = False
            continue
        if not str(arg).startswith("-"):
            raise WorkflowConfigError(
                f"micapipe would loop forever on the stray argument {arg!r} in: "
                + " ".join(shlex.quote(str(a)) for a in args)
            )
        expecting_value = arg in value_options
    if expecting_value:
        raise WorkflowConfigError(f"the last option {args[-1]!r} has no value")


def shell_join(args: Iterable[str]) -> str:
    return " ".join(shlex.quote(str(arg)) for arg in args)
