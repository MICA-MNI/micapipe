# Server Build Constraints & Requirements

**CRITICAL INFORMATION - READ THIS FIRST**

## Server Environment

**Server:** venice.bic.mni.mcgill.ca  
**User:** eyang  
**Build Location:** `/host/cassio/export03/data/enning/downloads`

---

## ⚠️ CONSTRAINTS (MUST REMEMBER!)

### 1. **NO SUDO ACCESS**
- User `eyang` does **NOT** have sudo/root privileges
- Cannot run any command requiring `sudo`
- Cannot install system packages on host
- Cannot modify `/etc/` files on host
- Cannot clean `/var/cache/apt/archives/` on host

**Implications:**
- ❌ Cannot use: `sudo apt-get`, `sudo rm`, `sudo chmod` on host
- ✅ CAN use: Docker commands, user-owned files/directories
- ✅ CAN clean: Docker resources, temp files in user space

### 2. **HOME DIRECTORY DISK SPACE**
- `~/` (home directory) has **LIMITED disk space** (~50-100GB)
- **NEVER** build Docker images in `~/micapipe`
- **ALWAYS** build in `/host/cassio/export03/data/enning/downloads` (unlimited space)

**Implications:**
- ❌ Don't run: `docker build` from `~/micapipe`
- ✅ DO run: `docker build` from `/host/cassio/export03/data/enning/downloads`
- Build script enforces this with error check

### 3. **DOCKER PERMISSIONS**
- User `eyang` IS in docker group (can run docker commands)
- Can build, run, stop, remove Docker containers/images
- Can clean Docker build cache, networks, volumes

**Implications:**
- ✅ CAN run: `docker build`, `docker system prune`, `docker rmi`
- ✅ CAN clean: Docker-managed resources

---

## 🔧 Space Management (Without Sudo)

### What You CAN Clean:

```bash
# Docker resources (most effective)
docker system prune -a -f           # Remove all unused Docker data
docker builder prune -f             # Clear build cache
docker volume prune -f              # Remove unused volumes (careful!)
docker image prune -a -f            # Remove unused images

# User-owned temp files
rm -rf /tmp/docker-*
rm -rf /tmp/build-*
rm -rf ~/.cache/*

# Build directory temp files
cd /host/cassio/export03/data/enning/downloads
rm -rf tmp/ temp/ *.log build_*.log
```

### What You CANNOT Clean (Requires Sudo):

```bash
# These will FAIL without sudo:
sudo rm -rf /var/cache/apt/archives/*     # ❌ No sudo
sudo apt-get clean                        # ❌ No sudo
sudo rm -rf /var/lib/apt/lists/*          # ❌ No sudo
```

### Pre-Build Cleanup Routine:

```bash
# Always run before building:
cd /host/cassio/export03/data/enning/downloads
./cleanup_docker_space.sh   # Cleans Docker resources (no sudo needed)
```

---

## 📁 Directory Structure

```
/home/eyang/                                    # HOME - LIMITED SPACE ⚠️
  └── micapipe/                                 # Git repo - DO NOT BUILD HERE!

/host/cassio/export03/data/enning/              # UNLIMITED SPACE ✅
  └── downloads/                                # BUILD HERE!
      ├── micapipe/                             # Migrated code
      ├── Dockerfile.base
      ├── build_base_image_server.sh
      ├── cleanup_docker_space.sh
      ├── fsl-6.0.2-centos6_64.tar.gz          # Pre-downloaded (2.6GB)
      ├── freesurfer-*.tar.gz                   # Pre-downloaded (2.7GB)
      ├── afni-linux_openmp_64.tgz             # Pre-downloaded (800MB)
      └── Miniconda3-*.sh                       # Pre-downloaded (400MB)
```

---

## 🚫 Common Mistakes to Avoid

### Mistake 1: Building in Home Directory
```bash
# ❌ WRONG - will run out of space
cd ~/micapipe
docker build -f Dockerfile.base .

# ✅ CORRECT - unlimited space
cd /host/cassio/export03/data/enning/downloads
docker build -f Dockerfile.base .
```

### Mistake 2: Using Sudo Commands
```bash
# ❌ WRONG - no sudo access
sudo apt-get install something
sudo docker system prune

# ✅ CORRECT - no sudo needed
docker system prune -f
```

### Mistake 3: Manual Confirmation Prompts
```bash
# ❌ WRONG - requires manual input
read -p "Continue? (y/N): "

# ✅ CORRECT - fully automatic
echo "⚠️ Warning: continuing automatically (Ctrl+C to cancel)..."
```

---

## 🎯 Build Script Requirements

All build scripts MUST:

1. **Validate build location** (not home directory)
2. **No sudo commands** anywhere
3. **No manual prompts** (fully automatic)
4. **Clean Docker resources** before building
5. **Use pre-downloaded files** from build directory

Example validation:
```bash
# In build script
if [[ "$PWD" == *"$HOME"* ]]; then
    echo "❌ ERROR: Building from HOME directory!"
    exit 1
fi
```

---

## 📋 Pre-Build Checklist

Before ANY Docker build:

- [ ] Building from `/host/cassio/export03/data/enning/downloads`
- [ ] NOT building from `~/micapipe` (home)
- [ ] Run `./cleanup_docker_space.sh` first
- [ ] Pre-downloaded files present in current directory
- [ ] No sudo commands in scripts
- [ ] No manual confirmation prompts
- [ ] Git branch up to date (`git pull`)

---

## 🆘 Troubleshooting

### "No space left on device"
1. Run `./cleanup_docker_space.sh`
2. Check you're in `/host/cassio/.../downloads` not `~/`
3. Run `docker system df` to see Docker space usage
4. Run `df -h` to see overall disk usage

### "Permission denied" 
- You're trying to do something requiring sudo
- Check if operation is in "CANNOT Clean" list above
- Use Docker cleanup instead of system cleanup

### "Cannot write to ..."
- Path is outside your permissions
- Don't try to modify system directories
- Stay in user space or Docker-managed space

---

## 💡 Remember

**GOLDEN RULE:** If you don't have sudo, you can't modify system files. 
But you CAN:
- Build Docker images (they're your containers)
- Clean Docker resources (they're yours)
- Manage files in your directories
- Run anything in containers (containers are isolated)

**Docker handles system operations INSIDE containers** - that's why we can install packages in Dockerfile even without sudo on host!

---

**Last Updated:** October 7, 2025  
**Maintainer:** Remember these constraints ALWAYS!
