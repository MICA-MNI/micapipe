# Changes from Working October Version (79464ef)

## Critical Differences That Could Cause Issues

### 1. **Python Version Upgrade** (REQUIRED for LAMAReg)
- **Old**: Python 3.9 (Oct version)
- **New**: Python 3.11 (current)
- **Why**: LAMAReg requires Python >=3.10
- **Status**: ✅ Correctly implemented

### 2. **Miniconda Version**
- **Old**: `Miniconda3-py39_22.11.1-1` (Nov 2022)
- **New**: `Miniconda3-py311_23.10.0-1` (Oct 2023)
- **Why**: Python 3.11 support, GLIBC 2.27 compatibility
- **Status**: ✅ Compatible with Ubuntu 18.04

### 3. **System Package Installation Strategy**
- **Old**: One large RUN command with all apt packages
- **New**: Split into 4 separate RUN commands (Essential, Dev tools, Graphics, Scientific)
- **Why**: Prevent OOM (Out of Memory) during build
- **Potential Issue**: More Docker layers, slightly larger image
- **Status**: ⚠️ Should help with memory but adds layers

### 4. **Package Installation Strategy**
- **Old**: Large batches with `--threads 16 --retry-delay 3 --retry-attempts 2`
- **New**: Smaller chunks, removed threading options
- **Why**: Prevent OOM, better error isolation
- **Potential Issue**: Slightly slower but more stable
- **Status**: ✅ Better for reliability

### 5. **Removed Package Version Pins**
- **Old**: `dipy==1.4.1`, `fury==0.8.0`, specific numpy/scipy versions
- **New**: Let conda-forge resolve versions for Python 3.11
- **Why**: Old pins incompatible with Python 3.11
- **Potential Issue**: May get different versions than expected
- **Status**: ⚠️ Necessary but needs testing

### 6. **Removed Options**
- **Old**: `--threads 16 --retry-delay 3 --retry-attempts 2` on mamba install
- **New**: No threading options
- **Why**: Caused instability in some environments
- **Status**: ✅ More stable

### 7. **Removed libmamba Solver Configuration**
- **Old**: `conda config --set solver libmamba || echo "⚠️ libmamba solver not available"`
- **New**: Removed this line
- **Why**: Potentially unstable, mamba already uses libmamba
- **Status**: ✅ Simplification

### 8. **New Tools Added** (per requirements)
- **synb0-DISCO**: For DWI without reverse phase encoding
- **SynBOLD-DisCo**: For fMRI without reverse phase encoding
- **Status**: ✅ New requirements met

### 9. **New pip Packages**
- **antspyx**: Added explicitly (was in old version via different method)
- **brainspace==0.1.10**: Neuroimaging analysis
- **tedana==0.0.12**: fMRI denoising
- **duecredit**: Citation management
- **vtk==9.2.2**: Fixed version (was unpinned)
- **Status**: ✅ Enhanced functionality

### 10. **GUI/X11 Libraries**
- **New**: Dedicated RUN layer for X11/visualization libraries
- **Why**: Support for pyvista, fury, vtk visualization
- **Status**: ✅ Better organization

### 11. **Server Path Change**
- **Old**: `/host/cassio/export03/data/enning`
- **New**: `/export03/data/enning`
- **Why**: Different server or mount point
- **Status**: ⚠️ Verify this is correct for bb-compxg-01

### 12. **Designer Environment**
- **Old**: Python 3.8
- **New**: Python 3.11
- **Why**: Consistency with main environment
- **Status**: ⚠️ Verify DESIGNER works with Python 3.11

## Non-Issues (Working as Expected)

1. ✅ `COPY . /tmp/build_context/` - Still present, correct
2. ✅ Pre-download file checks - Same logic, different filenames
3. ✅ FSL, FreeSurfer, AFNI installation - Unchanged
4. ✅ ANTs, FastSurfer, SWM - Unchanged
5. ✅ R installation - Unchanged

## Potential Bugs to Watch For

### High Priority
1. **Docker Storage Full** - New issue since October
   - **Fix**: Run `./check_docker_space.sh` and clean up
   
2. **Server Path** - Changed from `/host/cassio/export03/data/enning` to `/export03/data/enning`
   - **Verify**: Correct path for bb-compxg-01 server
   
3. **Package Versions** - No longer pinned
   - **Risk**: May get incompatible versions
   - **Fix**: Test thoroughly after build

### Medium Priority
4. **DESIGNER with Python 3.11** - Upgraded from 3.8
   - **Risk**: DESIGNER scripts may expect Python 3.8
   - **Fix**: Test DESIGNER after build

5. **More Docker Layers** - Split RUN commands
   - **Risk**: Slightly larger image, more layers
   - **Impact**: Minimal, worth it for stability

### Low Priority
6. **No Threading Options** - Removed `--threads 16`
   - **Risk**: Slower builds
   - **Impact**: More stable, acceptable tradeoff

## Recommendations

### Before Building
1. ✅ Run `./check_docker_space.sh` on server
2. ✅ Clean up Docker: `docker system prune -a`
3. ✅ Verify correct Miniconda file: `Miniconda3-py311_23.10.0-1-Linux-x86_64.sh`
4. ⚠️ Verify server path is correct for your system

### After Building
1. Test LAMAReg works with Python 3.11
2. Test DESIGNER with Python 3.11
3. Test synb0-DISCO and SynBOLD-DisCo
4. Verify package versions are compatible
5. Run full micapipe test suite

## Build Script Changes

### Added in build_comprehensive_base_server.sh
- `DOCKER_TMPDIR=/export03/data/enning/docker_tmp` - Use server storage for Docker temp

## Summary

The current version has:
- ✅ **Python 3.11** - Required for LAMAReg
- ✅ **Better memory management** - Split RUN commands
- ✅ **New tools** - synb0-DISCO, SynBOLD-DisCo
- ✅ **Enhanced packages** - antspyx, brainspace, tedana
- ⚠️ **Docker storage issue** - NEW problem since October
- ⚠️ **Path change** - Verify `/export03/data/enning` is correct
- ⚠️ **Unpinned versions** - Necessary for Python 3.11 but needs testing

**Main Issue**: Docker storage full (not related to code changes)
**Solution**: Clean Docker storage on server before building
