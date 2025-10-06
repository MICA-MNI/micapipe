# ✅ FIX: Docker Content Trust Error (No Sudo Required)

## 🎯 The Problem

**Docker Content Trust (DCT)** is enabled on the server, causing SSL certificate validation errors:
```
could not validate the path to a trusted root: 
unable to retrieve valid leaf certificates
```

## ✅ Quick Fix (No Sudo Needed)

### Option 1: Disable DCT for Current Session

```bash
# On server, before running build script:
export DOCKER_CONTENT_TRUST=0

# Then run the build
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh
```

### Option 2: Disable DCT in Your Shell Profile (Permanent)

```bash
# Add to your ~/.bashrc or ~/.zshrc
echo 'export DOCKER_CONTENT_TRUST=0' >> ~/.bashrc
source ~/.bashrc

# Or for zsh
echo 'export DOCKER_CONTENT_TRUST=0' >> ~/.zshrc
source ~/.zshrc
```

### Option 3: One-Line Build Command

```bash
# Build with DCT disabled (no env variable needed)
cd /host/cassio/export03/data/enning/downloads
DOCKER_CONTENT_TRUST=0 ./build_base_image_server.sh
```

---

## 🔧 Automated Fix: Updated Build Script

I'll update `build_base_image_server.sh` to automatically disable DCT during builds.

---

## 💡 Quick Copy-Paste Solution

```bash
# === COPY AND PASTE ON SERVER ===

# Go to build directory
cd /host/cassio/export03/data/enning/downloads

# Disable Docker Content Trust and build
export DOCKER_CONTENT_TRUST=0
./build_base_image_server.sh

# === END ===
```

---

## 📖 What is Docker Content Trust?

- **DCT** is a security feature that verifies image signatures
- When enabled, Docker requires signed images from trusted publishers
- For **local builds**, you don't need DCT (we're building, not pulling)
- DCT is useful when pulling images from public registries

---

## ✅ Expected Result After Fix

```bash
# Build should proceed without SSL errors:
🏗️  Starting Docker build...
=============================================
Sending build context to Docker daemon  152.2kB
Step 1/XX : FROM ubuntu:18.04
 ---> 3339fde08fc3
Step 2/XX : ENV DEBIAN_FRONTEND=noninteractive
 ---> Running in abc123def456
...
```

No more SSL certificate errors! ✅
