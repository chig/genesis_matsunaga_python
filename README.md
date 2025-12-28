# GENESIS

**GENeralized-Ensemble SImulation System** - A molecular dynamics simulation software for biomolecular systems.

## Python Interface (genepie)

The `genepie` package provides a Python interface to GENESIS analysis tools and the ATDYN MD engine.

### Installation

```bash
pip install genepie
```

### Quick Start

```python
from genepie import genesis_exe, SMolecule

# Load molecular structure
mol = SMolecule.from_file(pdbfile="protein.pdb", psffile="protein.psf")
print(f"Loaded {mol.num_atoms} atoms")

# Load trajectory and calculate RMSD
traj = genesis_exe.crd_convert(
    psffile="protein.psf",
    pdbfile="protein.pdb",
    dcdfile="trajectory.dcd",
    selection_group="an:CA",
)
rmsd = genesis_exe.rmsd_analysis(molecule=mol, trajectories=traj)
print(f"RMSD: {rmsd.mean():.2f} Å")

# Run MD simulation
energies, coords = genesis_exe.run_atdyn_md(
    prmtopfile="protein.prmtop",
    ambcrdfile="protein.inpcrd",
    nsteps=1000,
    ensemble="NVT",
    temperature=300.0,
)
```

### Available Analysis Functions

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

### MD Engine Functions

- `run_atdyn_md()` - Run MD simulation
- `run_atdyn_min()` - Run energy minimization
- `run_atdyn_md_isolated()` - Run MD in subprocess (crash-safe)
- `run_atdyn_min_isolated()` - Run minimization in subprocess

### Supported File Formats

| Format | Topology | Coordinates | Parameters |
|--------|----------|-------------|------------|
| AMBER | `prmtopfile` | `ambcrdfile` | (in prmtop) |
| GROMACS | `grotopfile` | `grocrdfile` | (in grotop) |
| CHARMM | `psffile` | `pdbfile`/`crdfile` | `parfile`, `strfile` |

## Building from Source

```bash
# Set up Python environment
python -m venv .venv
source .venv/bin/activate
pip install numpy

# Build GENESIS
autoreconf -fi
./configure --disable-mpi CC=gcc FC=gfortran
make -j$(nproc)

# Install Python package
pip install -e .
```

## Documentation

- [GENESIS Website](https://www.r-ccs.riken.jp/labs/cbrt/)
- [CLAUDE.md](CLAUDE.md) - Developer guide for Claude Code

## License

See LICENSE file for details.
