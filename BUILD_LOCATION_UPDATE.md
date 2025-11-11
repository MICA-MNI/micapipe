# Build Location Change - November 11, 2025

## Summary
Updated all Docker build scripts and Dockerfiles to use the new server location: `/export03/data/enning/downloads/`

## Previous Location
```
/export02/data/enning/downloads/
```

## New Location
```
/export03/data/enning/downloads/
```

## Files Modified

### 1. **Dockerfile.mamba-base**
- Updated `CUSTOM_TMPDIR` build argument
- Updated compatibility comment

### 2. **build_comprehensive_base_server.sh**
- Updated expected server location checks
- Updated build argument `CUSTOM_TMPDIR`
- Updated all path references in messages

### 3. **build_main_image_server.sh**
- Updated expected server location checks
- Updated all path references in comments and messages

### 4. **prepare_build_context.sh**
- Updated `SOURCE_DIR` variable

### 5. **migrate_comprehensive_base_to_server.sh**
- Updated `SERVER_BASE_DIR` variable
- Updated all path references in documentation

## Verification

✅ All references to `/export02/data/enning` have been removed from:
- Shell scripts (*.sh)
- Dockerfiles

✅ All references now point to `/export03/data/enning`

## Usage

The scripts should now be run from the new location:

```bash
cd /export03/data/enning/downloads
./build_comprehensive_base_server.sh
```

or

```bash
cd /export03/data/enning/downloads
./build_main_image_server.sh
```

## Notes

- All pre-downloaded files should be at: `/export03/data/enning/`
- Build context will be at: `/export03/data/enning/downloads/`
- Scripts will verify the correct location before building
