# Finding and Preparing Test Data for LAMAReg Tests

## 🎯 Quick Answer

**You need**: Outputs from a **completed micapipe run** (processed subject data), NOT raw BIDS input data.

## 📍 What the Tests Actually Need

The tests are designed to validate **registration steps** within the micapipe pipeline. They expect intermediate processed files that are generated DURING a micapipe run, specifically the files that get **registered together**.

### Understanding the Data Flow

```
Raw BIDS Data → micapipe Processing → Intermediate Files → Registration → Final Output
                                           ↑
                                    Tests check HERE
```

## 🗂️ Where to Find Test Data on Your Server

### Option 1: Use Existing Micapipe Derivatives (RECOMMENDED)

If you have subjects that have already been processed through micapipe, you can find the test data in the derivatives folder:

```bash
# Typical micapipe output structure
/path/to/derivatives/micapipe/
└── sub-001/
    ├── anat/
    │   ├── sub-001_space-nativepro_T1w.nii.gz          # ← For most registrations
    │   └── sub-001_space-nativepro_desc-brain_T1w.nii.gz
    ├── dwi/
    │   ├── sub-001_space-dwi_model-CSD_map-FOD.nii.gz   # ← For DWI test
    │   └── sub-001_space-dwi_desc-T1w_nativepro.nii.gz  # ← For DWI test
    ├── func/
    │   └── sub-001_task-rest_space-func_desc-brain_bold.nii.gz  # ← For func test
    ├── flair/
    │   └── sub-001_FLAIR_preproc.nii.gz                 # ← For FLAIR test
    └── xfm/
        └── [transformation files generated during registration]
```

### Option 2: Run Micapipe First to Generate Test Data

If you don't have processed data yet:

```bash
# Run micapipe on one test subject
micapipe -sub 001 \
         -out /path/to/derivatives \
         -bids /path/to/raw/BIDS \
         -proc_structural \
         -proc_dwi \
         -proc_func

# This will generate all the intermediate files you need
```

## 📋 Specific Files Needed for Each Test

### 1. DWI Registration Test (`test_dwi_registration.sh`)

**What it tests**: T1w → DWI space registration

**Files needed**:
```bash
test_data_dwi/
├── T1w_in_dwi_brain.nii.gz    # T1w image already transformed to DWI space (brain-extracted)
└── dwi_fod.nii.gz              # DWI fiber orientation distribution (FOD) map
```

**Where to find on server**:
```bash
# T1w in DWI space (brain-extracted)
derivatives/micapipe/sub-XXX/dwi/sub-XXX_space-dwi_desc-T1w_nativepro-brain.nii.gz

# DWI FOD (use b0 volume of FOD)
derivatives/micapipe/sub-XXX/dwi/sub-XXX_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz
```

### 2. Functional MRI Test (`test_func_registration.sh`)

**What it tests**: Functional MRI → T1nativepro registration

**Files needed**:
```bash
test_data_func/
├── func_brain.nii.gz          # Functional MRI (brain-extracted, mean volume)
└── T1_nativepro.nii.gz        # T1w nativepro
```

**Where to find on server**:
```bash
# Functional MRI brain
derivatives/micapipe/sub-XXX/func/sub-XXX_task-*_space-func_desc-brain_bold.nii.gz
# or the mean volume
derivatives/micapipe/sub-XXX/func/volumetric/sub-XXX_task-*_space-func_desc-mean_bold.nii.gz

# T1 nativepro
derivatives/micapipe/sub-XXX/anat/sub-XXX_space-nativepro_T1w.nii.gz
```

### 3. FLAIR Test (`test_flair_registration.sh`)

**What it tests**: FLAIR → T1nativepro registration

**Files needed**:
```bash
test_data_flair/
├── FLAIR_preproc.nii.gz       # Preprocessed FLAIR
└── T1_nativepro.nii.gz        # T1w nativepro
```

**Where to find on server**:
```bash
# FLAIR preprocessed
derivatives/micapipe/sub-XXX/anat/sub-XXX_space-flair_FLAIR.nii.gz

# T1 nativepro  
derivatives/micapipe/sub-XXX/anat/sub-XXX_space-nativepro_T1w.nii.gz
```

### 4. MPC Test (`test_mpc_registration.sh`)

**What it tests**: Microstructural map → FreeSurfer space registration

**Files needed**:
```bash
test_data_mpc/
├── microstructural_map.nii.gz # Any microstructural map (FA, MD, etc.)
└── T1_fsnative.nii.gz         # FreeSurfer T1 (orig.mgz converted)
```

**Where to find on server**:
```bash
# Microstructural map (e.g., FA)
derivatives/micapipe/sub-XXX/dwi/sub-XXX_space-dwi_model-DTI_map-FA.nii.gz

# FreeSurfer T1
derivatives/micapipe/sub-XXX/anat/surfaces/freesurfer/sub-XXX/mri/orig.mgz
# Convert with: mri_convert orig.mgz T1_fsnative.nii.gz
```

### 5. MPC-SWM Test (`test_mpc_swm_registration.sh`)

**What it tests**: Microstructural map → T1nativepro registration

**Files needed**:
```bash
test_data_mpc_swm/
├── microstructural_map.nii.gz # Any microstructural map
└── T1_nativepro.nii.gz        # T1w nativepro
```

**Where to find on server**:
```bash
# Same as MPC test but different reference
# Microstructural map
derivatives/micapipe/sub-XXX/dwi/sub-XXX_space-dwi_model-DTI_map-FA.nii.gz

# T1 nativepro
derivatives/micapipe/sub-XXX/anat/sub-XXX_space-nativepro_T1w.nii.gz
```

## 🛠️ Step-by-Step: Setting Up Test Data

### Step 1: Find a Processed Subject

```bash
# On your server, find a subject with complete processing
cd /path/to/derivatives/micapipe
ls -d sub-*

# Pick one subject (e.g., sub-001)
SUB="sub-001"
```

### Step 2: Create Test Data Directory

```bash
# Create a directory for test data
mkdir -p ~/lamareg_test_data
cd ~/lamareg_test_data
```

### Step 3: Copy Required Files

```bash
# Set up paths
DERIV="/path/to/derivatives/micapipe"
SUB="sub-001"
TEST_DIR="~/lamareg_test_data"

# For DWI test
mkdir -p $TEST_DIR/dwi
cp $DERIV/$SUB/dwi/${SUB}_space-dwi_desc-T1w_nativepro-brain.nii.gz \
   $TEST_DIR/dwi/T1w_in_dwi_brain.nii.gz
cp $DERIV/$SUB/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz \
   $TEST_DIR/dwi/dwi_fod.nii.gz

# For functional test
mkdir -p $TEST_DIR/func
cp $DERIV/$SUB/func/${SUB}_task-*_space-func_desc-brain_bold.nii.gz \
   $TEST_DIR/func/func_brain.nii.gz
cp $DERIV/$SUB/anat/${SUB}_space-nativepro_T1w.nii.gz \
   $TEST_DIR/func/T1_nativepro.nii.gz

# For FLAIR test (if available)
mkdir -p $TEST_DIR/flair
cp $DERIV/$SUB/anat/${SUB}_space-flair_FLAIR.nii.gz \
   $TEST_DIR/flair/FLAIR_preproc.nii.gz 2>/dev/null || echo "FLAIR not available"
cp $DERIV/$SUB/anat/${SUB}_space-nativepro_T1w.nii.gz \
   $TEST_DIR/flair/T1_nativepro.nii.gz

# For MPC tests
mkdir -p $TEST_DIR/mpc
cp $DERIV/$SUB/dwi/${SUB}_space-dwi_model-DTI_map-FA.nii.gz \
   $TEST_DIR/mpc/microstructural_map.nii.gz

# T1 nativepro for MPC-SWM
cp $DERIV/$SUB/anat/${SUB}_space-nativepro_T1w.nii.gz \
   $TEST_DIR/mpc/T1_nativepro.nii.gz

# FreeSurfer T1 for MPC (convert from mgz)
mri_convert $DERIV/$SUB/anat/surfaces/freesurfer/$SUB/mri/orig.mgz \
   $TEST_DIR/mpc/T1_fsnative.nii.gz
```

### Step 4: Verify Test Data

```bash
# Check what you have
find ~/lamareg_test_data -name "*.nii.gz" -exec ls -lh {} \;

# Quick validation
cd ~/lamareg_test_data
for dir in */; do
    echo "=== $dir ==="
    ls -lh $dir/*.nii.gz
    echo ""
done
```

## 🧪 Running the Tests

### Method 1: Test Individual Modules

```bash
cd /path/to/micapipe/tests/lamareg_tests

# Test DWI registration
./test_dwi_registration.sh ~/lamareg_test_data/dwi

# Test functional registration
./test_func_registration.sh ~/lamareg_test_data/func

# Test FLAIR registration
./test_flair_registration.sh ~/lamareg_test_data/flair

# Test MPC
./test_mpc_registration.sh ~/lamareg_test_data/mpc

# Test MPC-SWM
./test_mpc_swm_registration.sh ~/lamareg_test_data/mpc
```

### Method 2: Run All Tests at Once

```bash
cd /path/to/micapipe/tests/lamareg_tests

# Note: Each test expects its specific subdirectory
# So we need to point to parent directory
./run_all_tests.sh ~/lamareg_test_data
```

### Method 3: Syntax Validation Only (No Data Needed)

```bash
# Just check command syntax without running registration
./test_dwi_registration.sh /tmp/dummy

# This validates:
# ✓ LAMAReg is installed
# ✓ Command syntax is correct
# ✓ Arguments are properly formatted
# ✗ Skips actual registration (no input files)
```

## 📊 What to Expect

### With Real Data
```
========================================
DWI Registration Test
========================================

Starting tests...

✓ PASS: LAMAReg installation
✓ PASS: ANTs installation
✓ PASS: Input: Moving image
✓ PASS: Input: Fixed image
✓ PASS: LAMAReg command syntax
✓ PASS: Transformation chain

========================================
Test Summary
========================================
Total tests: 6
Passed: 6
Failed: 0

✓ All tests passed!
```

### Without Real Data (Syntax Check Only)
```
✓ PASS: LAMAReg installation
✓ PASS: LAMAReg command syntax
✗ FAIL: Input: Moving image - File not found

Total tests: 8
Passed: 6
Failed: 2

Note: Using synthetic test data. Replace with real data for accurate testing.
```

## 🎯 Quick Setup Script

Save this as `setup_test_data.sh`:

```bash
#!/bin/bash
# Quick setup script for LAMAReg test data

DERIV="$1"  # Path to micapipe derivatives
SUB="$2"    # Subject ID (e.g., sub-001)
OUT="${3:-./lamareg_test_data}"  # Output directory

if [ $# -lt 2 ]; then
    echo "Usage: $0 <derivatives_dir> <subject_id> [output_dir]"
    echo "Example: $0 /data/derivatives/micapipe sub-001 ~/test_data"
    exit 1
fi

mkdir -p $OUT/{dwi,func,flair,mpc}

echo "Setting up test data for $SUB..."

# DWI
echo "  Copying DWI files..."
cp $DERIV/$SUB/dwi/${SUB}_space-dwi_desc-T1w_nativepro-brain.nii.gz \
   $OUT/dwi/T1w_in_dwi_brain.nii.gz 2>/dev/null
cp $DERIV/$SUB/dwi/${SUB}_space-dwi_model-CSD_map-FOD_desc-wmNorm.nii.gz \
   $OUT/dwi/dwi_fod.nii.gz 2>/dev/null

# Functional
echo "  Copying functional files..."
cp $DERIV/$SUB/func/${SUB}_task-*_space-func_desc-brain_bold.nii.gz \
   $OUT/func/func_brain.nii.gz 2>/dev/null
cp $DERIV/$SUB/anat/${SUB}_space-nativepro_T1w.nii.gz \
   $OUT/func/T1_nativepro.nii.gz 2>/dev/null

# FLAIR
echo "  Copying FLAIR files..."
cp $DERIV/$SUB/anat/${SUB}_space-flair_FLAIR.nii.gz \
   $OUT/flair/FLAIR_preproc.nii.gz 2>/dev/null
cp $DERIV/$SUB/anat/${SUB}_space-nativepro_T1w.nii.gz \
   $OUT/flair/T1_nativepro.nii.gz 2>/dev/null

# MPC
echo "  Copying MPC files..."
cp $DERIV/$SUB/dwi/${SUB}_space-dwi_model-DTI_map-FA.nii.gz \
   $OUT/mpc/microstructural_map.nii.gz 2>/dev/null
cp $DERIV/$SUB/anat/${SUB}_space-nativepro_T1w.nii.gz \
   $OUT/mpc/T1_nativepro.nii.gz 2>/dev/null

# FreeSurfer T1
if [ -f "$DERIV/$SUB/anat/surfaces/freesurfer/$SUB/mri/orig.mgz" ]; then
    echo "  Converting FreeSurfer T1..."
    mri_convert $DERIV/$SUB/anat/surfaces/freesurfer/$SUB/mri/orig.mgz \
       $OUT/mpc/T1_fsnative.nii.gz 2>/dev/null
fi

echo ""
echo "Done! Test data ready in: $OUT"
echo ""
echo "Summary:"
find $OUT -name "*.nii.gz" | wc -l | xargs echo "  Files created:"
echo ""
echo "Run tests with:"
echo "  cd tests/lamareg_tests"
echo "  ./test_dwi_registration.sh $OUT/dwi"
echo "  ./test_func_registration.sh $OUT/func"
echo "  # etc..."
```

## 💡 Pro Tips

1. **Use a subject with ALL modalities** for comprehensive testing
2. **Check file sizes** - files should be reasonably large (>1MB typically)
3. **Verify headers** with `fslinfo` before running tests
4. **Start with one test** (e.g., DWI) before running full suite
5. **Keep test data small** - one subject is enough for validation

## ❓ Common Questions

**Q: Do I need raw BIDS data?**  
A: No! You need **processed** micapipe outputs.

**Q: Can I test without any data?**  
A: Yes! Tests will validate command syntax even without input files.

**Q: Which subject should I use?**  
A: Any subject that completed micapipe processing successfully.

**Q: Do I need all modalities?**  
A: No, test each modality you have. Skip tests for missing modalities.

**Q: How long do tests take?**  
A: Syntax validation: ~1 second. Full registration: ~15-25 min per test.
