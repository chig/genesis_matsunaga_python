# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GENESIS (GENeralized-Ensemble SImulation System) is a molecular dynamics simulation software for biomolecular systems. This repository is focused on developing the **Python interface** for GENESIS analysis tools, located in `src/analysis/interface/python_interface/`.

## Build Commands

### Quick Start (Python Interface Development)

```bash
# Set up Python environment
uv venv --python=python3.11
source .venv/bin/activate
uv pip install torch torchvision torchaudio nglview numpy mdtraj MDAnalysis plotly jupyterlab py3Dmol scikit-learn gdown

# Build GENESIS (first time or after configure.ac changes)
autoreconf -fi
./configure --disable-mpi CC=gcc-14 FC=gfortran LAPACK_LIBS="-L/usr/local/lib -llapack -lblas"
make -j$(sysctl -n hw.ncpu)

# Install Python package (editable mode)
pip install -e .
```

### Rebuilding After Changes

| Change Type | Required Action |
|-------------|-----------------|
| Python files (.py) | No rebuild needed (instant reflection) |
| Fortran files (.fpp) | `make` |
| New Fortran files added | `make clean && make` |
| configure.ac changes | `autoreconf -fi && ./configure ... && make` |

### Demo Notebook

```bash
uv run jupyter lab demo/demo.ipynb
```

## PyPI Package (genepy)

The package is distributed via PyPI as `genepy`. Users can install with:

```bash
pip install genepy
```

### Package Contents

| Component | Description |
|-----------|-------------|
| `genepy` Python package | Python interface with `libpython_interface.so` |
| `atdyn` | MD engine (CLI command) |
| 43 analysis tools | CLI commands (rmsd_analysis, trj_analysis, etc.) |

**Note**: `spdyn` requires MPI and is NOT included in PyPI. Users needing spdyn should build from source.

### Key Package Files

| File | Purpose |
|------|---------|
| `pyproject.toml` | Package metadata, dependencies, CLI entry points |
| `setup.py` | Forces platform-specific wheel (not pure Python) |
| `src/genepy/cli.py` | CLI entry point functions for all binaries |
| `MANIFEST.in` | Specifies files to include in source distribution |

### CLI Commands Available After Installation

```bash
# MD Engine
atdyn INP

# Analysis Tools (examples)
rmsd_analysis INP
trj_analysis INP
mbar_analysis INP
wham_analysis INP
crd_convert INP
kmeans_clustering INP
# ... and 37 more
```

## GitHub Actions (CI/CD)

### Workflow: `.github/workflows/build-wheels.yml`

Multi-stage build process for PyPI wheel distribution:

| Stage | Purpose | Platforms |
|-------|---------|-----------|
| `build-genesis` | Build Fortran binaries (atdyn, analysis tools, libpython_interface) | Linux x86_64, macOS arm64, macOS x86_64 |
| `build-wheels` | Create Python wheels with bundled binaries | Same as above |
| `build-sdist` | Create source distribution | Linux |
| `publish` | Upload to PyPI (on tag push `v*`) | Linux |
| `publish-testpypi` | Upload to TestPyPI (on workflow_dispatch) | Linux |

### Triggers

- **Push tags `v*`**: Build and publish to PyPI
- **Pull requests to main**: Build and test (no publish)
- **workflow_dispatch**: Build and publish to TestPyPI

### Publishing to PyPI

```bash
# Create and push a version tag
git tag v0.1.0
git push origin v0.1.0
# GitHub Actions will automatically build and publish
```

### Testing on TestPyPI

Use "Run workflow" button on GitHub Actions page to trigger `publish-testpypi`.

## Python Interface Architecture

### Directory: `src/analysis/interface/python_interface/`

The Python interface uses ctypes to call Fortran functions compiled into `libpython_interface.so`.

### Core Components

| File | Purpose |
|------|---------|
| `genesis_exe.py` | Main API - Python functions calling GENESIS analysis tools |
| `libgenesis.py` | C type definitions for Fortran subroutine bindings |
| `libloader.py` | Loads `libpython_interface.so` shared library |
| `s_molecule.py` | `SMolecule` class - Python representation of molecular structure |
| `s_molecule_c.py` | `SMoleculeC` class - ctypes Structure for C interface |
| `s_trajectories.py` | `STrajectories` class - trajectory data handling |
| `s_trajectories_c.py` | `STrajectoriesC` class - ctypes Structure for C interface |
| `c2py_util.py` | C data → numpy array conversion utilities |
| `py2c_util.py` | numpy array → C data conversion utilities |
| `ctrl_files.py` | Generates temporary control files for GENESIS functions |

### Fortran Interface Files (.fpp)

| Pattern | Purpose |
|---------|---------|
| `*_c_mod.fpp` | C-callable wrappers for GENESIS Fortran routines |
| `*_analysis.fpp` | Analysis algorithm implementations |
| `conv_f_c_util.fpp` | Fortran ↔ C data conversion utilities |
| `s_molecule_c_mod.fpp` | s_molecule C structure allocation/deallocation |
| `s_trajectories_c_mod.fpp` | s_trajectories C structure handling |

### Available Analysis Functions (in `genesis_exe.py`)

- `crd_convert()` - Coordinate/trajectory conversion
- `trj_analysis()` - Distance, angle, dihedral analysis
- `rmsd_analysis()` - RMSD calculation
- `drms_analysis()` - Distance RMSD calculation
- `rg_analysis()` - Radius of gyration
- `msd_analysis()` - Mean squared displacement
- `diffusion_analysis()` - Diffusion coefficient calculation
- `hb_analysis()` - Hydrogen bond analysis
- `avecrd_analysis()` - Average coordinate calculation
- `wham_analysis()` - WHAM free energy analysis
- `mbar_analysis()` - MBAR free energy analysis
- `kmeans_clustering()` - K-means trajectory clustering

### Integration with MDTraj and MDAnalysis

`SMolecule` has conversion methods:
- `to_mdtraj_topology()` / `from_mdtraj_topology()`
- `to_mdanalysis_universe()` / `from_mdanalysis_universe()`

## Testing

### Run All Python Interface Tests
```bash
cd src/analysis/interface/python_interface
./all_run.sh
```

### Run Individual Tests
```bash
python rmsd_analysis.py
python trj_analysis.py
python mbar_analysis_umbrella_1d.py
# etc.
```

## Adding a New Analysis Function

1. Create Fortran wrapper in `*_c_mod.fpp` with `bind(C)` interface
2. Create analysis implementation in `*_analysis.fpp` (or reuse existing GENESIS code)
3. Add function signature to `libgenesis.py`
4. Create Python wrapper function in `genesis_exe.py`
5. Write regression test as `<analysis_name>.py`
6. Add test to `all_run.sh`
7. Update `Makefile.am` to include new `.fpp` files

## Key Patterns

### Calling Fortran from Python
```python
# In genesis_exe.py
mol_c = molecule.to_SMoleculeC()  # Convert Python → C structure
try:
    with suppress_stdout_simple():  # Suppress Fortran output
        LibGenesis().lib.some_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(result_c),
            # ... other args
        )
    result = c2py_util.conv_double_ndarray(result_c, size)  # C → numpy
finally:
    LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
```

### Control File Generation
Analysis functions use temporary control files (GENESIS INP format):
```python
with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
    ctrl_files.write_ctrl_output(ctrl, rmsfile="dummy.rms")
    ctrl_files.write_ctrl_selection(ctrl, selection_group, selection_mole_name)
    # ... call Fortran function with ctrl.name
```

## Code Structure (Full Repository)

- `src/spdyn/` - Domain decomposition MD engine (MPI parallel)
- `src/atdyn/` - Atom decomposition MD engine
- `src/lib/` - Shared Fortran library
- `src/analysis/` - Analysis tools
  - `libana/` - Analysis library
  - `interface/python_interface/` - **Main development focus**
  - `trj_analysis/`, `free_energy/`, `mode_analysis/`, etc.
- `src/genepy/` - Python package for PyPI distribution
