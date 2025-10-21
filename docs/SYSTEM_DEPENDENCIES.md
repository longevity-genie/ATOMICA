# System Dependencies

This document describes system-level dependencies that may be required for certain ATOMICA features.

## X11 Libraries (for RDKit visualization)

Some RDKit visualization features require X11 libraries. These are **only needed if you use molecule visualization functions** (e.g., `molecule.to_SVG()`).

### When X11 libraries are NOT needed

The following operations work without X11 libraries:
- PDB file downloading (`uv run pdb download`)
- PDB file conversion (`uv run pdb convert`)
- Protein embedding generation
- Most data processing tasks

### When X11 libraries ARE needed

Only needed for:
- Molecule visualization (`.to_SVG()` method in `data/tokenizer/molecule.py`)
- Creating molecular structure images

### Installing X11 Libraries

If you encounter an error like:
```
ImportError: libXrender.so.1: cannot open shared object file: No such file or directory
```

Install the required libraries based on your system:

#### Rocky Linux / CentOS / RHEL
```bash
sudo yum install libXrender libXext libSM libXt
```

#### Ubuntu / Debian
```bash
sudo apt-get install libxrender1 libxext6 libsm6 libxt6
```

#### Fedora
```bash
sudo dnf install libXrender libXext libSM libXt
```

#### Alpine Linux
```bash
apk add libxrender libxext libsm libxt
```

### For Headless Servers

If you're running on a headless server and need visualization, you can use a virtual X server:

```bash
# Rocky Linux / CentOS / RHEL
sudo yum install xorg-x11-server-Xvfb

# Ubuntu / Debian
sudo apt-get install xvfb

# Then run your command with xvfb-run:
xvfb-run uv run python your_script.py
```

### Docker Considerations

If you're using Docker, add the following to your Dockerfile:

```dockerfile
# For Rocky Linux / CentOS base images
RUN yum install -y libXrender libXext libSM libXt

# For Ubuntu / Debian base images
RUN apt-get update && apt-get install -y libxrender1 libxext6 libsm6 libxt6
```

## Implementation Notes

The RDKit visualization import (`rdkit.Chem.Draw.rdMolDraw2D`) is implemented as a lazy import in `data/tokenizer/molecule.py`. This means it's only loaded when actually needed, allowing other ATOMICA features to work without X11 libraries installed.

If you encounter issues with other dependencies, please open an issue on the project repository.

