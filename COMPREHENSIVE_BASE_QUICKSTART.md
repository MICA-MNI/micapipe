# MICApipe Comprehensive Base Image Strategy - Quick Start

## 🎯 **95% Faster CI Builds**

This strategy separates neuroimaging tools (stable) from MICApipe code (frequently changing) for ultra-fast development cycles.

## 🚀 **Quick Start**

### **Option 1: Interactive Workflow (Recommended)**
```bash
cd ~/micapipe
git checkout comprehensive-base-image
./workflow.sh
```
Choose from the interactive menu:
- **5** for complete setup (first time)
- **6** for development workflow (regular use)

### **Option 2: Manual Steps**
```bash
# 1. Migrate source to server (from ~/micapipe)
./migrate_comprehensive_base_to_server.sh

# 2. Build comprehensive base image (45-90 min, one-time)
cd /host/cassio/export03/data/enning/downloads
./build_comprehensive_base_server.sh

# 3. Fast CI builds (3-5 min, regular use)
./build_fast_ci_server.sh
```

## 📁 **Architecture**

```
~/micapipe (development)  →  /host/cassio/export03/data/enning/downloads (build)
     ↓                                           ↓
[Source Code]              →              [Source + Pre-downloads]
     ↓                                           ↓
Comprehensive Base Image (45-90 min, one-time)
     ↓
Fast Main Image Builds (3-5 min, always)
```

## 🔧 **Key Benefits**

✅ **95% Build Time Reduction** for code changes  
✅ **Uses Pre-Downloaded Files** from server  
✅ **Unlimited Space** on server for builds  
✅ **Home Directory** for development  
✅ **Automatic Sync** between environments  

## 📚 **Files Overview**

- `workflow.sh` - Interactive workflow manager
- `migrate_comprehensive_base_to_server.sh` - Source migration script
- `Dockerfile.mamba-base` - Comprehensive neuroimaging base image
- `Dockerfile.minimal` - Fast CI main image
- `build_comprehensive_base_server.sh` - Base image builder
- `build_fast_ci_server.sh` - Fast CI builder
- `SERVER_DEPLOYMENT_GUIDE.md` - Complete documentation

## 🎯 **Performance**

| Build Type | Time | Frequency |
|------------|------|-----------|
| Traditional | 45-90 min | Every change |
| Base Image | 45-90 min | One-time only |
| **Fast CI** | **3-5 min** | **Every change** |

🚀 **Start with `./workflow.sh` for the easiest experience!**