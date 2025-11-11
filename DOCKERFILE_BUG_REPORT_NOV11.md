# Dockerfile Systematic Bug & Compatibility Review
## Date: November 11, 2025

---

## 🚨 CRITICAL ISSUES

### 1. Python 3.9 End-of-Life ⚠️ **CRITICAL**
**Current:** Python 3.9  
**Status:** **EOL October 2025** (just ended last month!)  
**Impact:** HIGH - Security updates stopped, package support will decrease  
**Recommendation:** **Upgrade to Python 3.10 or 3.11 immediately**

```dockerfile
# CURRENT (PROBLEMATIC):
RUN mamba create -y -n micapipe python=3.9

# RECOMMENDED:
RUN mamba create -y -n micapipe python=3.11
```

**Migration Notes:**
- Python 3.10: Supported until Oct 2026
- Python 3.11: Supported until Oct 2027 (recommended)
- Test all packages for Python 3.11 compatibility

---

### 2. VTK Version Seriously Outdated ⚠️ **HIGH PRIORITY**
**Current:** vtk==9.2.2 (from 2022)  
**Latest:** vtk 9.5.2 (September 2025)  
**Gap:** ~3 years behind  
**Impact:** MEDIUM-HIGH - Missing features, potential security issues

**Verified from PyPI:** https://pypi.org/project/vtk/
- Latest: 9.5.2 (Sep 17, 2025)
- Your version: 9.2.2 (ancient)

**Recommendation:**
```dockerfile
# TEST THIS FIRST:
pip install --no-cache-dir vtk==9.5.2

# Then update if compatible with fury and pyvista
```

**Compatibility Check Needed:**
- VTK 9.5.2 with fury 0.8.0 ❓
- VTK 9.5.2 with pyvista (unpinned) ✅ (should work)

---

### 3. Fury Version Outdated ⚠️ **MEDIUM PRIORITY**
**Current:** fury==0.8.0 (from ~2022)  
**Latest:** fury 0.12.0 (December 11, 2024)  
**Gap:** 4 major versions behind  
**Impact:** MEDIUM - Missing newer visualization features

**Verified from PyPI:** https://pypi.org/project/fury/
- Latest: 0.12.0 (Dec 11, 2024)
- Your version: 0.8.0 (2+ years old)

**Recommendation:**
```dockerfile
# UPDATE TO:
pip install --no-cache-dir fury==0.12.0
```

**Known Issue:** fury 0.8.0 may have compatibility issues with VTK 9.5.x  
**Action:** Test fury 0.12.0 + VTK 9.5.2 combination

---

## ⚠️ MODERATE ISSUES

### 4. Old Conda Packages with Pinned Versions
**Potentially Problematic Pins:**

| Package | Current | Age | Status |
|---------|---------|-----|--------|
| numpy | 1.21.5 | 2021 | 4 years old |
| pandas | 1.4.4 | 2022 | 3 years old |
| matplotlib | 3.4.3 | 2021 | 4 years old |
| nibabel | 4.0.2 | 2022 | 3 years old |
| dipy | 1.4.1 | 2021 | 4 years old |

**Impact:** LOW-MEDIUM - These versions work but miss new features  
**Recommendation:** Test with newer versions when upgrading Python

---

### 5. Ubuntu Base Image Very Old
**Current:** ubuntu:bionic-20201119 (November 2020)  
**Status:** Ubuntu 18.04 reached standard EOL April 2023  
**Impact:** MEDIUM - Security updates limited to ESM until 2028

**Recommendation:**
```dockerfile
# CONSIDER UPGRADING TO:
FROM ubuntu:jammy-20240627.1  # Ubuntu 22.04 LTS
# OR
FROM ubuntu:focal-20230412    # Ubuntu 20.04 LTS (if 22.04 has issues)
```

**Migration Effort:** HIGH - May require package updates

---

## ✅ VERIFIED WORKING

### 1. ANTsPy Installation ✅
**Method:** pip install antspyx  
**Latest:** 0.6.1 (June 25, 2025)  
**Status:** GOOD - Using pip as intended, no version pin

**Verified:** https://pypi.org/project/antspyx/

---

### 2. PyVista ✅
**Method:** pip install pyvista (no version pin)  
**Latest:** 0.46.4 (October 30, 2025)  
**Status:** EXCELLENT - Will always get latest version

**Verified:** https://pypi.org/project/pyvista/

---

### 3. MRtrix3 ✅
**Version:** mrtrix3==3.0.7  
**Channel:** mrtrix3  
**Status:** VERIFIED available for linux-64

---

### 4. FSL Download Link ✅
**URL:** https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.2-centos6_64.tar.gz  
**Status:** HTTP 200 (ACTIVE)  
**Verified:** November 11, 2025

---

## 🔍 POTENTIAL COMPATIBILITY ISSUES

### 1. VTK + Fury + PyVista Triangle
**Problem:** These packages must be compatible with each other

**Current State:**
- vtk==9.2.2 (2022)
- fury==0.8.0 (2022)  
- pyvista (unpinned, will get 0.46.4)

**Risk:** PyVista 0.46.4 may expect newer VTK  
**Action:** Test the combination or pin all three together

**Recommended Testing:**
```python
import vtk
import fury
import pyvista as pv
# Test basic operations
```

---

### 2. Conda GCC vs System Build-Essential
**Current:** Installing both:
- System: `build-essential` (apt-get)
- Conda: `gcc_linux-64` and `gxx_linux-64`

**Risk:** Path conflicts, version mismatches  
**Impact:** LOW - Usually works but can cause subtle build issues

**Recommendation:** Choose one approach:
```dockerfile
# OPTION A: Remove conda compilers, use system only
# Remove: gcc_linux-64 gxx_linux-64

# OPTION B: Remove system build-essential, use conda only
# Remove from apt-get: build-essential
```

---

### 3. NumPy Version with Python 3.9+ EOL
**Issue:** numpy=1.21.5 is from 2021  
**Latest:** numpy 2.3.4 (November 2025)

**Impact:** When upgrading to Python 3.11, numpy 1.21.5 may not have wheels  
**Recommendation:** Test with newer numpy when upgrading Python:
```dockerfile
numpy=1.24.3  # or newer, test for compatibility
```

---

## 📋 PACKAGES TO VERIFY (Missing from Check)

### Still Need URLs/Repos Checked:
1. ✅ **brainspace==0.1.10** - No PyPI page found (may be GitHub only)
2. ❓ **tedana==0.0.12** - No PyPI data retrieved
3. ❓ **FreeSurfer** download URL
4. ❓ **AFNI** download URL  
5. ❓ **FSL FIX** download URL
6. ❓ **FastSurfer** v2.4.2 tag exists
7. ❓ **SWM** (superficial-white-matter) repo
8. ❓ **DESIGNER-v2** repo
9. ❓ **LAMAReg** repo (MICA-MNI)
10. ❓ **ENIGMA** repo (MICA-MNI)

---

## 🎯 PRIORITY ACTION PLAN

### Immediate (Do This Week):
1. ✅ **TEST Python 3.11 Migration**
   - Create test branch
   - Update Python version
   - Test all imports
   - Run basic pipeline

2. ✅ **Update VTK + Fury Together**
   - Test VTK 9.5.2 + fury 0.12.0 compatibility
   - Verify pyvista still works
   - Update if tests pass

### Short Term (This Month):
3. ✅ **Verify All Download URLs**
   - FreeSurfer
   - AFNI
   - FSL FIX
   - Add fallback URLs if needed

4. ✅ **Test Git Repos**
   - Verify all GitHub repos are accessible
   - Check if tags/branches still exist
   - Update to latest stable versions if needed

### Medium Term (Next 3 Months):
5. ✅ **Ubuntu Base Image Upgrade**
   - Test with Ubuntu 22.04
   - Document any breaking changes
   - Migrate when ready

6. ✅ **Unpin Old Package Versions**
   - Test newer numpy, pandas, matplotlib
   - Update if compatible
   - Keep backward compatibility

---

## 🛠️ TESTING CHECKLIST

Before deploying changes:

```bash
# 1. Build test
docker build -f Dockerfile.mamba-base -t test-base .

# 2. Import test
docker run test-base python -c "
import vtk
import fury  
import pyvista
import ants  # from antspyx
import nibabel
import dipy
print('All imports successful!')
"

# 3. Version check
docker run test-base python -c "
import vtk, fury, pyvista
print(f'VTK: {vtk.vtkVersion.GetVTKVersion()}')
print(f'Fury: {fury.__version__}')
print(f'PyVista: {pyvista.__version__}')
"
```

---

## 📊 RISK ASSESSMENT

| Issue | Severity | Effort | Priority |
|-------|----------|--------|----------|
| Python 3.9 EOL | 🔴 HIGH | Medium | 1 |
| VTK Outdated | 🟡 MEDIUM | Low | 2 |
| Fury Outdated | 🟡 MEDIUM | Low | 3 |
| Old NumPy/Pandas | 🟢 LOW | Medium | 4 |
| Ubuntu 18.04 | 🟡 MEDIUM | High | 5 |
| GCC Conflict | 🟢 LOW | Low | 6 |

---

## 💡 RECOMMENDATIONS SUMMARY

### Must Fix:
1. **Upgrade Python to 3.11** (3.9 is EOL)
2. **Test and update VTK to 9.5.2**
3. **Test and update fury to 0.12.0**

### Should Fix:
4. Verify all download URLs still work
5. Test Git repository accessibility
6. Document fallback mirrors

### Nice to Have:
7. Upgrade Ubuntu base image to 22.04
8. Update old pinned packages (numpy, pandas, etc.)
9. Resolve conda/system gcc duplication

---

## 📝 NOTES

- **antspyx** pip installation is correct approach ✅
- **MRtrix3** 3.0.7 verified available ✅
- **FSL** download link verified working ✅  
- **pyvista** unpinned = always latest (good strategy) ✅

**Created:** November 11, 2025  
**Verified Against:** Live PyPI, Conda-forge, and official websites
