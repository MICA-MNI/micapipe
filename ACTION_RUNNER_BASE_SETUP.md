# Action Runner Base Image Setup

## 🎯 **Problem**
The new CI uses `base + main` Docker structure, but your action runner needs the base image to exist first.

## 🛠️ **Solution**

### **Step 1: Prepare your action runner (one-time setup)**

```bash
# On your server, run this once:
git pull origin comprehensive-base-image
./prepare_action_runner.sh
```

**What this does:**
- Checks if your action runner has the base image
- Builds the base image inside the runner if missing (~45-60 minutes one-time)
- Tests that the base image works

### **Step 2: CI becomes fast**

After base image is prepared:
- ✅ **CI checks**: "Base image exists? Yes!"
- ⚡ **CI builds**: Only main image (~5-10 minutes)
- 🚀 **Total CI time**: ~15-20 minutes instead of 60+ minutes

## 📋 **Timeline**

### **First time (with base image build):**
```
CI Run 1: 60+ minutes (builds base + main + tests)
```

### **Subsequent runs (base cached):**
```
CI Run 2+: 15-20 minutes (builds main + tests)
```

## 🔍 **How it Works**

### **Action Runner Setup:**
1. Action runner contains: Ubuntu + Singularity + GitHub runner
2. **NEW**: Base image cached inside runner
3. CI jobs reuse the base image

### **CI Workflow:**
```yaml
- name: Prepare base image
  run: |
    if docker image inspect "ghcr.io/mica-mni/micapipe-base:latest"; then
      echo "✅ Base image found - CI will be fast!"
    else  
      echo "❌ Building base image (one-time, ~45-60 min)"
      docker build -f Dockerfile.base -t ghcr.io/mica-mni/micapipe-base:latest .
    fi

- name: Build main image (always fast ~5-10 min)
  run: docker build -f Dockerfile.main ...
```

## 🚀 **Commands Summary**

```bash
# 1. One-time setup (run on your server)
./prepare_action_runner.sh

# 2. CI automatically becomes fast
# (no more commands needed)
```

## ✅ **Benefits**

- **One-time setup**: Base image built once, reused forever
- **Fast CI**: 15-20 minutes instead of 60+ minutes  
- **No registry dependency**: Everything local in action runner
- **Backward compatible**: Existing runner setup unchanged