# Fix for Rocky Linux X11 Dependency Error

## The Problem

When running `uv run pdb download 1tgr` on Rocky Linux (Google Cloud), you encountered:
```
ImportError: libXrender.so.1: cannot open shared object file: No such file or directory
```

This happened because RDKit's drawing module requires X11 libraries, even though the PDB download command doesn't need visualization features.

## The Solution

### Option 1: Pull the Latest Code (Recommended)

The issue has been fixed in the codebase. Simply pull the latest changes:

```bash
cd ATOMICA
git pull
```

The fix makes the RDKit drawing import "lazy" - it only loads when actually needed for visualization, not for simple PDB operations.

### Option 2: Install System Libraries (If You Need Visualization)

If you need molecule visualization features, install the required X11 libraries on Rocky Linux:

```bash
sudo yum install libXrender libXext libSM libXt
```

## After the Fix

The following commands will work **without** installing X11 libraries:
- `uv run pdb download <pdb_id>` - Download PDB structures
- `uv run pdb convert <file>` - Convert PDB/CIF files
- Most data processing and embedding operations

Visualization features (like `.to_SVG()` in molecule processing) will only require X11 libraries when actually used.

## More Information

For detailed information about system dependencies, see: [docs/SYSTEM_DEPENDENCIES.md](docs/SYSTEM_DEPENDENCIES.md)

## Verification

Test that it works:
```bash
cd ATOMICA
uv run pdb download 1tgr
```

This should now work without requiring X11 libraries!

