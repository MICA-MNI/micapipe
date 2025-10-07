# QUICK START - Fixed Docker Build

## Problem FIXED ✅
"no space left on device" error → Build happens ON SERVER with more space

## What Changed
- ⚠️  **CRITICAL: Builds happen ON SERVER, NOT in home directory!**
- Build location: `/host/cassio/export03/data/enning/downloads/`
- Pre-downloaded files + build files in SAME directory (simple!)

## Run These Commands

### 1. Pull Latest Fix (on local machine)
```bash
cd ~/micapipe
git checkout comprehensive-base-image
git pull origin comprehensive-base-image
```

### 2. Free Up Server Space (if needed)
```bash
cd ~/micapipe
# Copy cleanup script to server
cp cleanup_server_space.sh /host/cassio/export03/data/enning/downloads/
cd /host/cassio/export03/data/enning/downloads
./cleanup_server_space.sh
```

### 3. Migrate to Server (copies files TO server)
```bash
cd ~/micapipe
./migrate_comprehensive_base_to_server.sh
```

### 4. Build Base Image (⚠️  FROM SERVER LOCATION!)
```bash
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh
```

⏱️ **Expected time:** 45-90 minutes

## ⚠️ CRITICAL - Build Location

**✅ CORRECT (on server with space):**
```bash
cd /host/cassio/export03/data/enning/downloads
./build_base_image_server.sh
```

**❌ WRONG (in home directory - will fail!):**
```bash
cd ~/micapipe  # DON'T DO THIS!
docker build ...  # Will cause "no space left" error
```

## ⚠️ Important Notes

1. **⚠️  MUST build from server location** (`/host/.../downloads/` NOT `~/micapipe`)
2. **Uses pre-downloaded files** (FSL, FreeSurfer, AFNI already on server)
3. **Script will ERROR if you try to build from home directory**
4. **Run cleanup script if low on space**

## Commits Applied
- `0b34b56` - Separate build directory fix
- `cc0347a` - Documentation  
- `80010f7` - GPG key import with resilient fallbacks (fixes keyserver errors)
- `1efe3bb` - Revert to simple approach + cleanup script
- `LATEST` - Enforce build location checks (prevents home directory builds)

## What to Expect

✅ **Build runs from server location** (has space + pre-downloaded files)  
✅ **Uses pre-downloaded files** (FreeSurfer 2.7GB, FSL 2.6GB already there)  
✅ **No "no space left on device" errors**  
✅ **Build completes successfully** after 45-90 minutes  
✅ **Script blocks if you try to build from home directory**  

## If Issues Persist

1. Check you're in SERVER directory (NOT home):
   ```bash
   pwd  # Should be: /host/cassio/export03/data/enning/downloads
   ```

2. Run cleanup script:
   ```bash
   cd /host/cassio/export03/data/enning/downloads
   ./cleanup_server_space.sh
   ```

3. Check disk space on server:
   ```bash
   df -h /host/cassio/export03/data/enning
   ```

4. Clean Docker:
   ```bash
   docker system prune -a
   ```

## Success!
When build completes, you'll have:
- `micapipe-base:latest` (Stage 1 - comprehensive base with all tools)
- Ready to build Stage 2 with `./build_main_image_server.sh` (3-5 min)

---
📖 Full details: [DISK_SPACE_FIX_SUMMARY.md](./DISK_SPACE_FIX_SUMMARY.md)
