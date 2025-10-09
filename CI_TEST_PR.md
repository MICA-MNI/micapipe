# CI Test - Pull Request Trigger

This file was created to test the CI workflow with pull request triggers.

## Test Details
- **Branch**: test-ci-pr-trigger
- **Purpose**: Verify CI runs on PR instead of push
- **Expected**: CI should build using existing base image quickly (~10-15 minutes)

## What this tests:
1. ✅ CI triggers on pull request only
2. ✅ Uses existing ghcr.io/mica-mni/micapipe-base:latest image
3. ✅ Builds main image quickly (3-5 minutes)
4. ✅ Creates Singularity image
5. ✅ Runs micapipe tests

## After CI completes:
- Check GitHub Actions tab for build time
- Should complete much faster than v1 (10-15 min vs 60+ min)
- Base image should be reused, not rebuilt