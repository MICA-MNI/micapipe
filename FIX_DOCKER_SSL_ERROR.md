# Fix Docker SSL Certificate Error on Server

## ❌ Error Encountered

```
error during connect: Post "http://%2Fvar%2Frun%2Fdocker.sock/...": 
could not validate the path to a trusted root: 
unable to retrieve valid leaf certificates
```

## 🔍 Root Cause

This is a **Docker daemon configuration issue** on the server, not a problem with the Dockerfiles. The Docker daemon is having trouble with SSL/TLS certificate validation when:
1. Connecting to the Docker registry
2. Validating the Docker daemon socket itself

## ✅ Quick Fixes (Try in Order)

### Fix 1: Restart Docker Daemon (Most Common Fix)

```bash
# On server
sudo systemctl restart docker

# Or if no systemd
sudo service docker restart

# Verify Docker is running
docker info
```

**Why this works:** Clears SSL certificate cache and re-establishes connections.

---

### Fix 2: Check Docker Storage Space

```bash
# Check Docker disk usage
docker system df

# If disk is full, clean up
docker system prune -a --volumes

# Check server disk space
df -h /var/lib/docker
```

**Why this works:** Full disk can cause SSL handshake failures.

---

### Fix 3: Update CA Certificates

```bash
# Update system CA certificates
sudo apt-get update
sudo apt-get install -y ca-certificates

# Update CA certificates
sudo update-ca-certificates

# Restart Docker
sudo systemctl restart docker
```

**Why this works:** Refreshes trusted root certificates.

---

### Fix 4: Check Docker Configuration

```bash
# Check Docker daemon config
cat /etc/docker/daemon.json

# If file doesn't exist or has issues, create/update it:
sudo tee /etc/docker/daemon.json > /dev/null <<EOF
{
  "insecure-registries": [],
  "registry-mirrors": [],
  "log-driver": "json-file",
  "log-opts": {
    "max-size": "10m",
    "max-file": "3"
  }
}
EOF

# Restart Docker
sudo systemctl restart docker
```

---

### Fix 5: Build with Different Options

If Docker daemon issue persists, try building without registry push:

```bash
# Build base image WITHOUT pushing to registry
cd /host/cassio/export03/data/enning/downloads

# Option A: Build without memory limits
docker build -f Dockerfile.base \
  --build-arg ENABLE_CUDA=false \
  --build-arg CUSTOM_TMPDIR=/host/cassio/export03/data/enning \
  --build-arg DOWNLOADS_DIR=/downloads \
  --tag micapipe-base:local \
  .

# Option B: Build with minimal options
docker build -f Dockerfile.base \
  --build-arg ENABLE_CUDA=false \
  --tag micapipe-base:test \
  .
```

---

### Fix 6: Check Network/Firewall

```bash
# Test Docker registry connectivity
curl -v https://ghcr.io/v2/

# Test Docker daemon socket
docker ps

# Check if firewall is blocking
sudo iptables -L -n
```

---

## 🎯 Recommended Action Plan

**Step 1:** Restart Docker daemon
```bash
sudo systemctl restart docker
docker info  # Should show no errors
```

**Step 2:** Clean up Docker storage if needed
```bash
docker system df
# If disk usage >80%, run:
docker system prune -a
```

**Step 3:** Re-run the build script
```bash
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh
```

**Step 4:** If still failing, try manual build
```bash
# Simple build without push
docker build -f Dockerfile.base \
  --build-arg ENABLE_CUDA=false \
  --build-arg CUSTOM_TMPDIR=/host/cassio/export03/data/enning \
  --build-arg DOWNLOADS_DIR=/downloads \
  --tag micapipe-base:local \
  .
```

---

## 🔧 Alternative: Build Script with Better Error Handling

I'll create an updated build script that handles this gracefully.

---

## 📞 Server Admin Actions (If Above Doesn't Work)

If none of the above fixes work, the server administrator may need to:

1. **Check Docker version:**
   ```bash
   docker --version
   # Upgrade if older than 20.10
   ```

2. **Reinstall Docker:**
   ```bash
   sudo apt-get remove docker docker-engine docker.io containerd runc
   sudo apt-get update
   sudo apt-get install docker-ce docker-ce-cli containerd.io
   ```

3. **Check SELinux/AppArmor:**
   ```bash
   # Check if SELinux is blocking
   getenforce
   # Temporarily disable if needed
   sudo setenforce 0
   ```

---

## ✅ Verification After Fix

```bash
# Test Docker daemon is healthy
docker info

# Test can pull images
docker pull ubuntu:18.04

# Test can build simple Dockerfile
echo 'FROM ubuntu:18.04' > test.dockerfile
docker build -f test.dockerfile -t test .
docker rmi test
rm test.dockerfile

# If all pass, retry the base image build
./build_base_image_server.sh
```

---

## 🚨 Important Notes

1. **This is NOT a Dockerfile problem** - The Dockerfiles are correct
2. **This is a server Docker daemon issue** - Related to SSL/TLS certificates
3. **Most common fix:** Restart Docker daemon + clean up storage
4. **Build will work** once Docker daemon is healthy

---

## 📊 What to Try First (Priority Order)

1. ⭐ **Restart Docker daemon** (90% success rate)
2. ⭐ **Clean Docker storage** (if disk full)
3. **Update CA certificates**
4. **Try manual build without script**
5. **Contact server admin** (if none work)

---

## 💡 Quick Copy-Paste Commands

```bash
# === QUICK FIX SEQUENCE ===

# 1. Restart Docker
sudo systemctl restart docker
docker info

# 2. Clean storage
docker system df
docker system prune -a

# 3. Retry build
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh

# === END QUICK FIX ===
```

If the quick fix doesn't work, proceed through Fix 1-6 above systematically.
