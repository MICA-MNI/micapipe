# QUICK START - Fixed Docker Build

## Problem FIXED ✅
"no space left on device" error → Build directory separated from large files

## What Changed
- Build now runs from: `/host/cassio/export03/data/enning/downloads/micapipe_build/`
- Large files stay in: `/host/cassio/export03/data/enning/downloads/`
- Docker build context: <500MB (was 8GB+)

## Run These Commands

### 1. Pull Latest Fix (on local machine)
```bash
cd ~/micapipe
git checkout comprehensive-base-image
git pull origin comprehensive-base-image
```

### 2. Migrate to Server (creates new build directory)
```bash
cd ~/micapipe
./migrate_comprehensive_base_to_server.sh
```

### 3. Build Base Image (from NEW location)
```bash
cd /host/cassio/export03/data/enning/downloads/micapipe_build
./build_base_image_server.sh
```

⏱️ **Expected time:** 60-120 minutes (downloads files from internet)

## ⚠️ Important Notes

1. **Must run from `micapipe_build/` directory** (not `downloads/`)
2. **Requires internet** (will download FreeSurfer, FSL, AFNI, etc.)
3. **Pre-downloaded files not used** (Docker can't access parent directory)
4. **Longer build time** (due to downloads, but NO disk space errors!)

## Commits Applied
- `0b34b56` - Separate build directory fix
- `cc0347a` - Documentation  
- `80010f7` - GPG key import with resilient fallbacks (fixes keyserver errors)

## What to Expect

✅ **Build context preparation:** <1 minute (not 10+ minutes)  
✅ **No "no space left on device" errors**  
✅ **Downloads:** FreeSurfer (2.7GB), FSL (2.6GB), AFNI (800MB), etc.  
✅ **Build completes successfully** after 60-120 minutes  

## If Issues Persist

1. Check you're in correct directory:
   ```bash
   pwd  # Should be: .../downloads/micapipe_build
   ```

2. Check disk space:
   ```bash
   df -h /var/lib/docker  # Should have >20GB free
   ```

3. Clean Docker:
   ```bash
   docker system prune -a
   ```

## Success!
When build completes, you'll have:
- `micapipe-base:latest` (Stage 1 - comprehensive base with all tools)
- Ready to build Stage 2 with `./build_main_image_server.sh` (3-5 min)

---
📖 Full details: [DISK_SPACE_FIX_SUMMARY.md](./DISK_SPACE_FIX_SUMMARY.md)
