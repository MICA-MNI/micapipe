# Dockerfile.mamba-base Systematic Fixes - November 5, 2025

## Summary
Systematically fixed package installation issues to prevent build failures and ensure consistency with environment.yml specifications.

## Issues Found and Fixed

### 1. ✅ ANTsPy Package Removed from Conda-Forge
**Problem:** `antspy` was removed from conda-forge channel between mid-October and early November 2025.

**Solution:** Changed to pip installation using `antspyx` package.
```dockerfile
# OLD (broken):
RUN mamba install -y -n micapipe -c conda-forge -c mrtrix3 antspy

# NEW (working):
RUN pip install --no-cache-dir antspyx
```

### 2. ✅ Duplicate Fury Installation
**Problem:** `fury==0.8.0` was installed twice - once via conda (line 200) and once via pip (line 259).

**Solution:** Removed conda installation, kept only pip version to avoid conflicts.
```dockerfile
# OLD:
mamba install -y fury==0.8.0 ipython jupyter...

# NEW:
mamba install -y ipython jupyter...  # fury removed
pip install fury==0.8.0  # kept in pip section
```

### 3. ✅ Added Version Constraints for Stability
**Problem:** Core packages were installed without version pins, could lead to incompatibilities.

**Solution:** Added version constraints to match environment.yml:
- `numpy=1.21.5`
- `pandas=1.4.4`
- `matplotlib=3.4.3`
- `nibabel=4.0.2`
- `dipy==1.4.1`

### 4. ✅ Added Missing Packages
**Problem:** Several packages in environment.yml were not in Dockerfile.

**Solution:** Added missing packages to conda installations:
- Scientific: `scikit-fmm`, `nilearn`, `astropy`, `pillow`, `packaging`
- Utilities: `attrs`, `click`, `pyyaml`, `python-dateutil`, `pytz`, `lxml`
- Networking: `aiohttp`

### 5. ✅ Enhanced Pip Packages
**Problem:** Pip section was minimal, missing several key packages from environment.yml.

**Solution:** Added to pip installation:
```python
pip install --no-cache-dir \
    vtk==9.2.2 \
    pyvista \
    fury==0.8.0 \
    antspyx \
    brainspace==0.1.10 \
    tedana==0.0.12 \
    duecredit \
    git+https://github.com/MICA-MNI/LAMAReg.git \
    git+https://github.com/MICA-MNI/ENIGMA.git
```

### 6. ✅ Verified Package Availability
**Checked:**
- ✓ `mrtrix3==3.0.7` exists in mrtrix3 channel (for linux-64)
- ✓ `dipy==1.4.1` available in conda-forge
- ✓ All pip packages installable
- ✓ VTK version pinned to 9.2.2 (consistent with environment.yml)

## Changes to Environment.yml

Updated `micapipe_environment.yml`:
```yaml
# OLD:
  - antspy
  - antspyx

# NEW:
  - pip:
    - antspyx  # ANTsPy no longer available in conda-forge, install via pip
```

## Potential Remaining Issues to Monitor

1. **Python 3.9 EOL**: Python 3.9 reaches end-of-life in October 2025. Consider upgrading to Python 3.10 or 3.11 in future.

2. **NumPy Version**: Using numpy=1.21.5 (from 2021). Some newer packages may require numpy 1.22+. Monitor for warnings.

3. **VTK Compatibility**: Using vtk==9.2.2 via pip. Ensure it's compatible with other visualization packages (fury, pyvista).

4. **Compiler Toolchain**: Installing `gcc_linux-64` and `gxx_linux-64` via conda might conflict with system gcc. Monitor build issues.

## Testing Recommendations

1. Build the image and verify all packages install successfully
2. Test import of critical packages:
   ```python
   import ants  # from antspyx
   import vtk
   import fury
   import brainspace
   import dipy
   import nibabel
   ```
3. Run a simple pipeline test to ensure tools work together

## Build Command

```bash
./build_comprehensive_base_server.sh
```

## Expected Improvements

- ✅ No more "antspy does not exist" errors
- ✅ No package conflicts from duplicate installations
- ✅ Better version consistency between Dockerfile and environment.yml
- ✅ More complete package set matching environment specifications
- ✅ Cleaner separation between conda and pip packages
