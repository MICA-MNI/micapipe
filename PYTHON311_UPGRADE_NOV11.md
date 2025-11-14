# Python 3.11 Upgrade & Additional Tools - November 11, 2025

## Executive Summary

This document tracks the major upgrade from Python 3.9 to Python 3.11 in the micapipe Docker container, along with the addition of critical neuroimaging tools for handling missing reverse phase encoding data.

## Motivation

### Critical Issues Addressed
1. **Python 3.9 End-of-Life**: Python 3.9 reached EOL in October 2025
   - No more security updates
   - Package maintainers dropping support
   - Security vulnerabilities will not be patched

2. **LAMAReg Compatibility**: LAMAReg requires Python >=3.10
   - `pyproject.toml` specifies: `requires-python = ">=3.10"`
   - Essential for cross-modality registration

3. **Future-Proofing**: Python 3.11 supported until October 2027
   - 2+ years of guaranteed security updates
   - Better package ecosystem support
   - Improved performance (10-60% faster than 3.9)

## Changes Implemented

### 1. Python Version Upgrade

#### Miniconda Environment
- **Before**: `python=3.9` with Miniconda3-py39_22.11.1-1
- **After**: `python=3.11` with Miniconda3-py311_25.9.1-3

#### Designer Environment
- **Before**: `python=3.8`
- **After**: `python=3.11`
- **Rationale**: Consistency across environments, modern package support

### 2. Miniconda Installer Update

**Old Installer:**
```dockerfile
https://repo.anaconda.com/miniconda/Miniconda3-py39_22.11.1-1-Linux-x86_64.sh
```

**New Installer:**
```dockerfile
https://repo.anaconda.com/miniconda/Miniconda3-py311_23.10.0-1-Linux-x86_64.sh
```

**Version Jump**: py39_22.11.1-1 (Nov 2022) → py311_23.10.0-1 (Oct 2023)
- Python 3.11 base
- Compatible with Ubuntu 18.04 (GLIBC 2.27)
- Note: Newer versions (24.x, 25.x) require GLIBC ≥2.28

### 3. New Tools Added

#### synb0-DISCO
- **Purpose**: Synthesize b0 images for DWI when reverse phase encoding is not present
- **Repository**: https://github.com/MASILab/Synb0-DISCO
- **Installation Path**: `/opt/Synb0-DISCO`
- **Environment Variable**: `SYNB0_DISCO_HOME="/opt/Synb0-DISCO"`

#### SynBOLD-DisCo
- **Purpose**: Synthesize reverse phase encoding for fMRI distortion correction
- **Repository**: https://github.com/MASILab/SynBOLD-DisCo
- **Installation Path**: `/opt/SynBOLD-DisCo`
- **Environment Variable**: `SYNBOLD_DISCO_HOME="/opt/SynBOLD-DisCo"`

**Use Case**: Both tools address common neuroimaging scenarios where reverse phase encoding data is unavailable, enabling distortion correction that would otherwise be impossible.

## Package Compatibility Verification

### Core Scientific Packages (Python 3.11 Compatible)

| Package | Current Version | Python 3.11 Support | Notes |
|---------|----------------|---------------------|-------|
| numpy | 1.21.5 → 2.3.4 | ✅ Full support | May need version bump |
| scipy | Latest | ✅ Full support | Works with Python 3.11 |
| pandas | 1.4.4 | ✅ Full support | Stable |
| matplotlib | 3.4.3 | ✅ Full support | Stable |
| vtk | 9.2.2 → 9.5.2 | ✅ Full support | 3 year gap identified |
| fury | 0.8.0 → 0.12.0 | ✅ Full support | 2+ year gap identified |
| pyvista | Latest (0.46.4) | ✅ Full support | Unpinned strategy good |
| nibabel | 4.0.2 | ✅ Full support | Stable |
| dipy | 1.4.1 | ✅ Full support | Stable |

### LAMAReg Dependencies (Python >=3.10 Required)
- antspyx ✅
- tensorflow ✅
- keras ✅
- statsmodels ✅
- All LAMAReg requirements verified compatible

### MRtrix3 & Neuroimaging Tools
- mrtrix3 3.0.7 ✅
- FreeSurfer 7.4.1 ✅
- FastSurfer 2.4.2 ✅
- AFNI ✅
- ANTs 2.3.4 ✅
- FSL 6.0.2 ✅

## Testing Checklist

### Critical Tests Required Before Production

- [ ] **Build Test**: Full Docker build completes without errors
- [ ] **Import Tests**: All Python packages import successfully
  ```python
  import numpy, scipy, pandas, matplotlib
  import nibabel, dipy, vtk, fury, pyvista
  import antspyx
  from lamareg import cli
  ```
- [ ] **Tool Availability**: All neuroimaging tools accessible
  ```bash
  which mri_convert  # FreeSurfer
  which mrconvert    # MRtrix3
  which 3dinfo       # AFNI
  which antsRegistration  # ANTs
  which fslinfo      # FSL
  ```
- [ ] **LAMAReg Test**: Run basic LAMAReg registration
  ```bash
  lamareg --version
  lamareg synthseg --help
  ```
- [ ] **synb0-DISCO Test**: Verify tool is accessible
- [ ] **SynBOLD-DisCo Test**: Verify tool is accessible
- [ ] **DESIGNER Test**: Activate designer environment and test
- [ ] **Pipeline Test**: Run minimal micapipe workflow

### Performance Tests
- [ ] Build time comparison (expect similar or faster)
- [ ] Runtime performance (Python 3.11 should be 10-60% faster)
- [ ] Memory usage (should be comparable or better)

## Known Version Gaps (Not Blocking but Recommended)

### VTK Version Gap
- **Current**: vtk==9.2.2 (from 2022)
- **Latest**: vtk==9.5.2 (Sep 17, 2025)
- **Gap**: 3 years
- **Recommendation**: Test upgrade to 9.5.2 in separate branch
- **Risk**: Moderate - may have API changes

### Fury Version Gap
- **Current**: fury==0.8.0 (from 2022)
- **Latest**: fury==0.12.0 (Dec 11, 2024)
- **Gap**: 2+ years
- **Recommendation**: Test upgrade to 0.12.0 after VTK upgrade
- **Risk**: Moderate - depends on VTK version

### NumPy Version Consideration
- **Current**: numpy=1.21.5 (pinned in Dockerfile)
- **Latest**: numpy 2.3.4
- **Note**: NumPy 2.0+ has breaking changes
- **Recommendation**: Stay on 1.x series until full testing
- **LAMAReg requirement**: numpy<2 (already handled)

## Neurodocker Startup Script Evaluation

**Decision**: KEEP the neurodocker/startup.sh script

**Rationale**:
- Provides standardized Docker entry point
- Minimal overhead (simple bash script)
- Enables flexible container execution
- Compatible with existing workflows
- No security or maintenance concerns

**Script Purpose**:
```bash
#!/bin/bash
set +e
if [ -n "$1" ]; then
    exec "$@"
else
    /bin/bash
fi
```
- Allows custom command execution
- Falls back to interactive bash
- Industry standard pattern

## File Changes Summary

### Modified Files
1. `Dockerfile.mamba-base`
   - Line 47-51: Updated Miniconda installer copy check
   - Line 158-164: Updated Miniconda installation
   - Line 176: Changed to `python=3.11` for micapipe env
   - Line 233: Changed to `python=3.11` for designer env
   - Added synb0-DISCO installation (lines ~375-380)
   - Added SynBOLD-DisCo installation (lines ~381-386)

### Pre-Download Requirements
Update your pre-download script to fetch:
```bash
# New Miniconda installer (compatible with Ubuntu 18.04 GLIBC 2.27)
wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.10.0-1-Linux-x86_64.sh
```

**Important:** Newer Miniconda versions (24.x, 25.x) require GLIBC ≥2.28 and are incompatible with Ubuntu 18.04 (bionic).

## Migration Notes

### For Server Builds
If you have pre-downloaded dependencies:
1. Download new Miniconda installer: `Miniconda3-py311_25.9.1-3-Linux-x86_64.sh`
2. Update build scripts to reference new installer
3. Old installer will be ignored; build will download if not found

### For CI/CD Builds
- First build after upgrade will take longer (rebuilding Python environments)
- Subsequent builds should be faster due to Python 3.11 performance
- Expect ~5-10 minutes additional time for initial mamba environment creation

## Rollback Plan

If Python 3.11 causes issues:

1. **Revert Python versions**:
   ```dockerfile
   python=3.11 → python=3.10  # Minimum for LAMAReg
   ```

2. **Revert Miniconda**:
   ```dockerfile
   Miniconda3-py311_23.10.0-1 → Miniconda3-py310_23.10.0-1
   ```

3. **Keep new tools**: synb0-DISCO and SynBOLD-DisCo are Python-version independent

**Note**: Python 3.9 is NOT recommended for rollback due to EOL status

## Timeline

| Date | Event |
|------|-------|
| Oct 2025 | Python 3.9 reaches EOL |
| Nov 5, 2025 | antspy conda-forge removal fixed |
| Nov 11, 2025 | Python 3.11 upgrade + new tools added |
| Oct 2027 | Python 3.11 EOL (future) |

## Related Documentation

- [DOCKERFILE_FIXES_NOV5.md](./DOCKERFILE_FIXES_NOV5.md) - antspy fix
- [BUILD_LOCATION_UPDATE.md](./BUILD_LOCATION_UPDATE.md) - Path migration
- [DOCKERFILE_BUG_REPORT_NOV11.md](./DOCKERFILE_BUG_REPORT_NOV11.md) - Systematic bug analysis

## References

1. LAMAReg Requirements: https://github.com/MICA-MNI/LAMAReg/blob/main/pyproject.toml
2. Python 3.9 EOL: https://devguide.python.org/versions/
3. Python 3.11 Release: https://docs.python.org/3.11/whatsnew/3.11.html
4. synb0-DISCO: https://github.com/MASILab/Synb0-DISCO
5. SynBOLD-DisCo: https://github.com/MASILab/SynBOLD-DisCo
6. VTK Releases: https://pypi.org/project/vtk/
7. Fury Releases: https://pypi.org/project/fury/

## Contact

For issues or questions:
- GitHub Issues: https://github.com/MICA-MNI/micapipe/issues
- LAMAReg Issues: https://github.com/MICA-MNI/LAMAReg/issues

---
**Status**: ✅ READY FOR TESTING
**Next Step**: Build Docker image and run test suite
**Priority**: HIGH (Python 3.9 EOL)
