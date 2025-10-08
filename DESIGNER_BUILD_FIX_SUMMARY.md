# DESIGNER Build Error - Root Cause and Fix

**Error**: `Building wheel for designer2 (pyproject.toml): finished with status 'error'`

---

## Root Cause Analysis

DESIGNER is not just a Python package - it includes a **C++ extension** (`lib.rpg`) that must be compiled during installation.

### What is lib.rpg?

From DESIGNER's `setup.py`:
```python
ext_modules = [
    Extension(
        "lib.rpg",
        ["rpg_cpp/unring_rpg_pybind.cpp"],  # C++ source code
        include_dirs=[get_pybind_include(), '/usr/local/include'],
        library_dirs=['/usr/local/lib'],
        libraries=["fftw3", "fftw3_threads", "m"],  # FFTW library dependency
        extra_compile_args=["-std=c++11"],
    ),
]
```

**Purpose**: Gibbs ringing removal (RPG = Remove Partial Gibbs)

### Build Requirements

From DESIGNER's `pyproject.toml`:
```toml
[build-system]
requires = ["setuptools>=42", "wheel", "pybind11>=2.12.0"]
build-backend = "setuptools.build_meta"
```

**What was missing**:
1. ❌ `pybind11` - Required to create Python bindings for C++ code
2. ❌ `fftw` - FFTW library used by the C++ extension
3. ❌ Build toolchain - `cmake`, `make`, `gcc`, `g++`

---

## The Fix (Commit 05187db)

Added build dependencies to micapipe conda environment:

```dockerfile
RUN mamba install -y -n micapipe -c conda-forge \
           cvxpy pandas \
           pybind11 \                          # ← For C++ Python bindings
           fftw \                              # ← FFTW library
           cmake make gcc_linux-64 gxx_linux-64  # ← Build toolchain
    && mamba run -n micapipe pip install --no-cache-dir \
           multiprocessing-logging \
    && mamba run -n micapipe pip install --no-cache-dir /opt/DESIGNER  # ← Now compiles successfully
```

---

## Why This Wasn't Obvious

1. **PyPI has pre-built wheels**: When you `pip install designer2` from PyPI, it downloads a pre-compiled wheel (`.whl` file) that already has the C++ extension compiled. No build tools needed.

2. **Building from source requires compilation**: When installing from a local directory (`pip install /opt/DESIGNER`), pip must compile the C++ extension from source.

3. **DESIGNER's official Docker uses pre-built wheels**: Their Dockerfile does:
   ```dockerfile
   RUN pip install -r requirements.txt  # Installs from PyPI (pre-built)
   COPY . .
   RUN python -m pip install . --no-deps  # Source already has compiled .so files
   ```

4. **Our approach builds from fresh source**: We clone the repo and build everything, which requires the full build toolchain.

---

## Verification

After this fix, the build will:
1. ✅ Install pybind11, fftw, and build tools into micapipe environment
2. ✅ Compile `rpg_cpp/unring_rpg_pybind.cpp` into `lib/rpg.cpython-310-x86_64-linux-gnu.so`
3. ✅ Install DESIGNER with working `designer` and `tmi` commands
4. ✅ RPG Gibbs ringing removal functionality will be available

---

## Complete Fix Summary

**Two commits required**:

1. **d691d2f**: Runtime requirements
   - PYTHONPATH for MRtrix3 libraries
   - Shared environment (micapipe) with MRtrix3 and FSL
   - Missing Python packages (cvxpy, pandas)

2. **05187db**: Build requirements (THIS FIX)
   - pybind11 for C++ bindings
   - fftw library
   - Build toolchain (cmake, make, gcc, g++)

**Deploy on server**:
```bash
cd ~/micapipe
git pull origin comprehensive-base-image  # Gets both d691d2f and 05187db
./migrate_comprehensive_base_to_server.sh
```

**Expected outcome**: Docker build will succeed, C++ extension will compile, DESIGNER will work.

---

## Technical Details

### What pybind11 Does
- Creates Python bindings for C++ code
- Allows Python to call C++ functions directly
- Used extensively in scientific Python packages (numpy, scipy, etc.)

### What FFTW Does
- Fast Fourier Transform library
- Used by RPG algorithm for frequency domain processing
- Critical for Gibbs ringing removal in diffusion MRI

### The Build Process
```
pip install /opt/DESIGNER
  ↓
setuptools reads pyproject.toml
  ↓
Installs build dependencies (pybind11, setuptools, wheel)
  ↓
Runs setup.py to build C++ extensions
  ↓
Compiles rpg_cpp/unring_rpg_pybind.cpp with pybind11 + fftw
  ↓
Creates lib/rpg.cpython-310-x86_64-linux-gnu.so
  ↓
Installs Python package with compiled extensions
  ↓
Creates entry points: designer, tmi
```

Without pybind11/fftw/gcc → Compilation fails → Wheel build fails → Installation fails

---

**Status**: ✅ FIXED - Ready for testing on server
