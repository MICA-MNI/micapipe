class BIDSName:
    """
    Builds a BIDS-compliant filename from entity key-value pairs and a suffix.

    Based on BIDS specification v1.11.1
    https://bids-specification.readthedocs.io/en/stable/appendices/entity-table.html
    https://bids-specification.readthedocs.io/en/stable/appendices/entities.html
    """

    # -------------------------------------------------------------------------
    # All valid suffixes, organised by modality (from the BIDS entity table).
    # A given workflow may wish to restrict this list to its own modality.
    # -------------------------------------------------------------------------
    ALLOWED_SUFFIXES = [
        # --- MRI : anatomical (anat) ---
        "T1w", "T2w", "PDw", "T2starw", "FLAIR", "inplaneT1", "inplaneT2",
        "PDT2", "angio", "T2star", "FLASH", "PD",
        # qMRI maps
        "T1map", "T2map", "T2starmap", "R1map", "R2map", "R2starmap",
        "PDmap", "MTRmap", "MTsat", "UNIT1", "T1rho", "MWFmap", "MTVmap",
        "Chimap", "S0map", "M0map",
        # qMRI source images / protocols
        "IRT1", "MEGRE", "MESE", "MP2RAGE", "MPM", "MTS", "MTR", "VFA",
        # other anat
        "defacemask",

        # --- MRI : diffusion (dwi) ---
        "dwi", "sbref", "ADC", "FA", "colFA", "expADC", "trace",
        # S0map already listed above

        # --- MRI : functional (func) ---
        "bold", "cbv", "phase", "noRF",
        # sbref already listed above

        # --- MRI : field maps (fmap) ---
        "phasediff", "phase1", "phase2",
        "magnitude", "magnitude1", "magnitude2",
        "fieldmap", "epi", "m0scan",
        # B1 field maps
        "TB1DAM", "TB1EPI", "TB1AFI", "TB1TFL", "TB1RFM", "TB1SRGE",
        "TB1map", "RB1COR", "RB1map",

        # --- MRI : perfusion (perf) ---
        "asl", "aslcontext", "asllabeling",
        # m0scan and noRF already listed above

        # --- Physiological / stimulus recordings (any modality) ---
        "physio", "physioevents", "stim",

        # --- Common cross-modality ---
        "events", "channels", "coordsystem", "electrodes", "photo",

        # --- EEG ---
        "eeg",

        # --- iEEG ---
        "ieeg",

        # --- MEG ---
        "meg", "markers", "headshape",

        # --- PET ---
        "pet", "blood",

        # --- NIRS ---
        "nirs", "optodes",

        # --- MRS ---
        "svs", "mrsi", "unloc", "mrsref",

        # --- Behavioural (beh) ---
        "beh",

        # --- Microscopy (micr) ---
        "TEM", "SEM", "uCT", "BF", "DF", "PC", "DIC", "FLUO", "CONF",
        "PLI", "CARS", "2PE", "MPE", "SR", "NLO", "OCT", "SPIM", "XPCT",

        # --- Motion ---
        "motion",
    ]

    # -------------------------------------------------------------------------
    # Canonical entity order as defined by the BIDS entity table (v1.11.1).
    # Entities that only appear in a subset of modalities are included so that
    # this single class can handle any datatype.
    #
    # Unified order (MRI as backbone; other-modality entities inserted at their
    # spec-defined positions relative to the shared entities):
    #
    #   sub → ses → sample → task → tracksys → acq → ce → trc → stain
    #   → nuc → voi → rec → dir → run → mod → echo → flip → inv → mt
    #   → part → proc → space → split → recording → chunk
    #
    # Entity key  | Full name                  | Modalities
    # ------------|----------------------------|---------------------------
    # sub         | Subject                    | all
    # ses         | Session                    | all
    # sample      | Sample                     | micr
    # task        | Task                       | func, eeg, ieeg, meg, beh, motion, perf, nirs, mrs
    # tracksys    | Tracking system            | motion
    # acq         | Acquisition                | all
    # ce          | Contrast enhancing agent   | anat, fmap
    # trc         | Tracer                     | pet
    # stain       | Stain                      | micr
    # nuc         | Nucleus                    | mrs
    # voi         | Volume of interest         | mrs
    # rec         | Reconstruction             | anat, dwi, fmap, perf, pet, mrs
    # dir         | Phase-encoding direction   | dwi, func, fmap, perf
    # run         | Run index                  | all
    # mod         | Corresponding modality     | anat (defacemask), func (noRF), perf (noRF)
    # echo        | Echo index                 | anat, func, fmap, mrs
    # flip        | Flip angle index           | anat, fmap
    # inv         | Inversion time index       | anat, fmap, mrs
    # mt          | Magnetization transfer     | anat, fmap
    # part        | Part (mag/phase/real/imag) | anat, dwi, func, fmap
    # proc        | Processed (on device)      | meg
    # space       | Space                      | eeg, ieeg, meg (derivatives)
    # split       | Split index                | meg
    # recording   | Recording label            | func, perf, dwi, anat, pet, beh, eeg, ieeg, meg, nirs, motion
    # chunk       | Chunk index                | anat, dwi, micr
    # -------------------------------------------------------------------------

    # Keys listed in canonical BIDS order; "suffix" is handled separately.
    ENTITY_ORDER = [
        "sub", "ses", "sample", "task", "tracksys",
        "acq", "ce", "trc", "stain", "nuc", "voi",
        "rec", "dir", "run", "mod",
        "echo", "flip", "inv", "mt", "part",
        "proc", "space", "split", "recording", "chunk",
    ]

    def __init__(self, **kwargs):
        """
        Pass entities and suffix as keyword arguments, e.g.:

            BIDSName(sub="01", ses="01", task="rest", suffix="bold").build()
            # → 'sub-01_ses-01_task-rest_bold'

        All entity keys must come from ENTITY_ORDER.
        The 'suffix' key is required and must be in ALLOWED_SUFFIXES.
        """
        unknown = set(kwargs) - set(self.ENTITY_ORDER) - {"suffix"}
        if unknown:
            raise ValueError(
                f"Unknown entity/key(s): {unknown}. "
                f"Allowed entities: {self.ENTITY_ORDER}"
            )

        if "suffix" not in kwargs:
            raise ValueError("'suffix' is required.")

        suffix = kwargs["suffix"]
        if suffix not in self.ALLOWED_SUFFIXES:
            raise ValueError(
                f"Invalid suffix '{suffix}'. "
                f"Must be one of:\n{self.ALLOWED_SUFFIXES}"
            )

        self.values = kwargs

    def build(self) -> str:
        """Return the underscore-joined BIDS filename stem (no extension)."""
        parts = []
        for key in self.ENTITY_ORDER:
            if key in self.values:
                value = self.values[key]
                parts.append(f"{key}-{value}")
        # Suffix appended last, without a key prefix
        parts.append(self.values["suffix"])
        return "_".join(parts)

    def __str__(self) -> str:
        return self.build()

    def __repr__(self) -> str:
        kv = ", ".join(f"{k}={v!r}" for k, v in self.values.items())
        return f"BIDSName({kv})"


# ---------------------------------------------------------------------------
# Quick smoke-test (remove or move to a test file for production use)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    examples = [
        BIDSName(sub="01", ses="01", task="rest", run="1", suffix="bold"),
        BIDSName(sub="01", acq="highres", suffix="T1w"),
        BIDSName(sub="01", ses="02", task="nback", acq="AP", dir="AP", run="2", echo="1", part="mag", suffix="bold"),
        BIDSName(sub="01", ses="01", acq="MTon", mt="on", suffix="MTR"),
        BIDSName(sub="01", ses="01", sample="A", acq="stitched", chunk="1", suffix="SPIM"),
        BIDSName(sub="01", task="rest", trc="FDG", run="1", suffix="pet"),
    ]
    for b in examples:
        print(b.build())