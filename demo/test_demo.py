#!/usr/bin/env python
"""
Comprehensive test script for GENESIS Python interface (genepie).

This script tests the main functionality demonstrated in demo.ipynb,
including molecular loading, trajectory processing, and various analyses.

Usage:
    python test_demo.py

Requirements:
    - chignolin.pdb, chignolin.psf, chignolin.dcd in the same directory
    - genepie package installed
"""

import os
import sys
import pathlib
import numpy as np
import warnings

# Change to script directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Test counter
tests_passed = 0
tests_failed = 0


def test(name, condition, details=""):
    """Helper function to run a test and report results."""
    global tests_passed, tests_failed
    if condition:
        print(f"  [PASS] {name}")
        tests_passed += 1
        return True
    else:
        print(f"  [FAIL] {name}")
        if details:
            print(f"         {details}")
        tests_failed += 1
        return False


def section(title):
    """Print section header."""
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


# =============================================================================
# Setup and Imports
# =============================================================================
section("1. Imports")

try:
    from genepie import (
        SMolecule, STrajectories, STrajectoriesArray,
        genesis_exe, LibGenesis, ctrl_files
    )
    from genepie import (
        GenesisError, GenesisFortranError, GenesisValidationError
    )
    test("genepie imports", True)
except ImportError as e:
    test("genepie imports", False, str(e))
    print("\nFATAL: Cannot import genepie. Exiting.")
    sys.exit(1)

# Define file paths
PDB_PATH = pathlib.Path("chignolin.pdb")
PSF_PATH = pathlib.Path("chignolin.psf")
DCD_PATH = pathlib.Path("chignolin.dcd")

# Check input files exist
test("chignolin.pdb exists", PDB_PATH.exists())
test("chignolin.psf exists", PSF_PATH.exists())
test("chignolin.dcd exists", DCD_PATH.exists())

if not all([PDB_PATH.exists(), PSF_PATH.exists(), DCD_PATH.exists()]):
    print("\nFATAL: Required input files not found. Exiting.")
    sys.exit(1)


# =============================================================================
# SMolecule Tests
# =============================================================================
section("2. SMolecule - Loading Molecular Information")

mol_allatom = SMolecule.from_file(pdb=PDB_PATH, psf=PSF_PATH)
test("SMolecule.from_file()", mol_allatom is not None)
test("num_atoms > 0", mol_allatom.num_atoms > 0,
     f"num_atoms = {mol_allatom.num_atoms}")
test("atom_name array", len(mol_allatom.atom_name) == mol_allatom.num_atoms)
test("atom_coord array", mol_allatom.atom_coord.shape == (mol_allatom.num_atoms, 3))

print(f"\n  Molecule info: {mol_allatom.num_atoms} atoms")

# Test subset_atoms
mol_subset = mol_allatom.subset_atoms(np.array([0, 1, 2]))
test("subset_atoms()", mol_subset.num_atoms == 3)


# =============================================================================
# crd_convert Tests
# =============================================================================
section("3. crd_convert - Trajectory Loading")

traj_params_list = [ctrl_files.TrajectoryParameters(trjfile=str(DCD_PATH))]

# Test 1: CA atom selection
trajs_array_ca, mol_ca = genesis_exe.crd_convert(
    molecule=mol_allatom,
    traj_params=traj_params_list,
    trj_format="DCD",
    trj_type="COOR+BOX",
    selection_group=["an:CA"],
    fitting_method="TR+ROT",
    fitting_atom=1,
    pbc_correct="NO"
)
traj_ca = trajs_array_ca[0]
test("CA selection - trajectory loaded", traj_ca.nframe > 0)
test("CA selection - correct atom count", traj_ca.natom == mol_ca.num_atoms)
ca_names = [''.join(n).strip() for n in mol_ca.atom_name]
test("CA selection - all atoms are CA", all(n == 'CA' for n in ca_names))

print(f"\n  CA trajectory: {traj_ca.nframe} frames, {traj_ca.natom} atoms")

# Test 2: PROA segment selection
trajs_array, mol = genesis_exe.crd_convert(
    molecule=mol_allatom,
    traj_params=traj_params_list,
    trj_format="DCD",
    trj_type="COOR+BOX",
    selection_group=["segid:PROA"],
    fitting_method="TR+ROT",
    fitting_atom=1,
    pbc_correct="NO"
)
traj = trajs_array[0]
test("PROA selection - trajectory loaded", traj.nframe > 0)
test("PROA selection - natom matches mol", traj.natom == mol.num_atoms)

print(f"  PROA trajectory: {traj.nframe} frames, {traj.natom} atoms")

# Test 3: Heavy atoms selection
trajs_heavy, mol_heavy = genesis_exe.crd_convert(
    molecule=mol_allatom,
    traj_params=traj_params_list,
    trj_format="DCD",
    trj_type="COOR+BOX",
    selection_group=["heavy"],
    fitting_method="TR+ROT",
    fitting_atom=1,
    pbc_correct="NO"
)
heavy_names = [''.join(n).strip() for n in mol_heavy.atom_name]
h_count = sum(1 for n in heavy_names if n.startswith('H'))
test("Heavy selection - no hydrogen atoms", h_count == 0)

print(f"  Heavy atoms: {mol_heavy.num_atoms} atoms")


# =============================================================================
# trj_analysis Tests
# =============================================================================
section("4. trj_analysis - Distance/Angle/Dihedral")

distances_def = ["PROA:1:GLY:CA  PROA:2:TYR:CA"]
angles_def = ["PROA:1:GLY:CA  PROA:2:TYR:CA  PROA:3:ASP:CA"]
torsions_def = ["PROA:1:GLY:CA  PROA:2:TYR:CA  PROA:3:ASP:CA  PROA:4:PRO:CA"]

analysis_results = genesis_exe.trj_analysis(
    molecule=mol, trajs=traj,
    distance=distances_def,
    angle=angles_def,
    torsion=torsions_def
)

test("distance result shape", analysis_results.distance.shape[0] == traj.nframe)
test("angle result shape", analysis_results.angle.shape[0] == traj.nframe)
test("torsion result shape", analysis_results.torsion.shape[0] == traj.nframe)
test("distance values reasonable", 2.0 < analysis_results.distance.mean() < 6.0,
     f"mean = {analysis_results.distance.mean():.3f}")

print(f"\n  Distance mean: {analysis_results.distance.mean():.3f} A")
print(f"  Angle mean: {analysis_results.angle.mean():.3f} deg")
print(f"  Torsion mean: {analysis_results.torsion.mean():.3f} deg")


# =============================================================================
# rg_analysis Tests
# =============================================================================
section("5. rg_analysis - Radius of Gyration")

rg_result = genesis_exe.rg_analysis(
    molecule=mol, trajs=traj,
    selection_group=["an:CA"],
    mass_weighted=True
)

test("rg result shape", rg_result.rg.shape[0] == traj.nframe)
test("rg values reasonable", 4.0 < rg_result.rg.mean() < 15.0,
     f"mean = {rg_result.rg.mean():.3f}")

print(f"\n  Rg mean: {rg_result.rg.mean():.3f} A")


# =============================================================================
# rmsd_analysis Tests
# =============================================================================
section("6. rmsd_analysis - RMSD")

rmsd_result = genesis_exe.rmsd_analysis(
    molecule=mol, trajs=traj,
    selection_group=["an:CA"],
    fitting_method="TR+ROT",
    fitting_atom=1,
    analysis_atom=1
)

test("rmsd result shape", rmsd_result.rmsd.shape[0] == traj.nframe)
test("rmsd first frame is small", rmsd_result.rmsd[0] < 3.0,
     f"first frame = {rmsd_result.rmsd[0]:.6f}")
test("rmsd values non-negative", np.all(rmsd_result.rmsd >= 0))

print(f"\n  RMSD mean: {rmsd_result.rmsd.mean():.3f} A")
print(f"  RMSD max: {rmsd_result.rmsd.max():.3f} A")


# =============================================================================
# avecrd_analysis Tests
# =============================================================================
section("7. avecrd_analysis - Average Structure")

avecrd_result = genesis_exe.avecrd_analysis(
    mol, traj,
    selection_group=["segid:PROA and heavy"],
    fitting_method="TR+ROT",
    fitting_atom=1,
    check_only=False,
    num_iterations=5,
    analysis_atom=1,
)

test("average PDB generated", avecrd_result.pdb is not None)
test("average PDB not empty", len(avecrd_result.pdb) > 100)
test("average PDB contains ATOM", "ATOM" in avecrd_result.pdb)

print(f"\n  Average PDB size: {len(avecrd_result.pdb)} chars")


# =============================================================================
# msd_analysis Tests
# =============================================================================
section("8. msd_analysis - Mean Square Displacement")

msd_result = genesis_exe.msd_analysis(
    molecule=mol, trajs=traj,
    selection_group=["an:CA"],
    oversample=True,
    delta=min(traj.nframe - 1, 1000)  # Limit for faster testing
)

test("msd result not None", msd_result.msd is not None)
test("msd values non-negative", np.all(msd_result.msd >= 0))

print(f"\n  MSD shape: {msd_result.msd.shape}")


# =============================================================================
# diffusion_analysis Tests
# =============================================================================
section("9. diffusion_analysis - Diffusion Coefficient")

diffusion_coeffs = genesis_exe.diffusion_analysis(
    msd_data=msd_result.msd,
    time_step=1.0,
    start="20%"
)

test("diffusion coeffs shape", diffusion_coeffs.shape[1] > 0)

print(f"\n  Diffusion coeffs shape: {diffusion_coeffs.shape}")


# =============================================================================
# wham_analysis Tests
# =============================================================================
section("10. wham_analysis - Free Energy")

# Create test CV data
os.makedirs("cv_wham_data", exist_ok=True)
num_cv_frames = 200
rows = np.arange(1, num_cv_frames + 1)
np.random.seed(42)  # For reproducibility
for i, loc in enumerate([1.0, 2.0, 3.0], 1):
    values = np.random.normal(loc=loc, scale=0.3, size=num_cv_frames)
    data = np.column_stack((rows, values))
    np.savetxt(f"cv_wham_data/window_{i}.dat", data, fmt="%d %.5f")

cv_file_template = "cv_wham_data/window_{:d}.dat"
pmf_wham = genesis_exe.wham_analysis(
    cvfile=cv_file_template,
    dimension=1,
    nblocks=1,
    temperature=300.0,
    tolerance=1.0E-7,
    rest_function=[1],
    grids=[(0.0, 4.0, 81)],
    constant=[(30.0, 30.0, 30.0)],
    reference=[(1.0, 2.0, 3.0)],
    is_periodic=[False]
)

test("wham result not None", pmf_wham is not None)
test("wham result has data", pmf_wham.shape[0] > 0 if pmf_wham is not None else False)

if pmf_wham is not None:
    print(f"\n  WHAM PMF shape: {pmf_wham.shape}")


# =============================================================================
# MDTraj Interface Tests
# =============================================================================
section("11. MDTraj Interface")

try:
    import mdtraj as md

    mdtraj_top = mol.to_mdtraj_topology()
    test("to_mdtraj_topology()", mdtraj_top.n_atoms == mol.num_atoms)

    mdtraj_traj = traj.to_mdtraj_trajectory(mol)
    test("to_mdtraj_trajectory()", mdtraj_traj.n_frames == traj.nframe)

    print(f"\n  MDTraj topology: {mdtraj_top.n_atoms} atoms")
    print(f"  MDTraj trajectory: {mdtraj_traj.n_frames} frames")

except ImportError:
    print("  [SKIP] MDTraj not installed")


# =============================================================================
# MDAnalysis Interface Tests
# =============================================================================
section("12. MDAnalysis Interface")

try:
    import MDAnalysis as mda

    mda_universe = mol.to_mdanalysis_universe()
    test("to_mdanalysis_universe()", mda_universe.atoms.n_atoms == mol.num_atoms)

    mda_traj_universe = traj.to_mdanalysis_universe(mol)
    test("traj.to_mdanalysis_universe()",
         mda_traj_universe.trajectory.n_frames == traj.nframe)

    print(f"\n  MDAnalysis universe: {mda_universe.atoms.n_atoms} atoms")
    print(f"  MDAnalysis trajectory: {mda_traj_universe.trajectory.n_frames} frames")

except ImportError:
    print("  [SKIP] MDAnalysis not installed")


# =============================================================================
# Scikit-learn Integration Tests
# =============================================================================
section("13. Scikit-learn Integration")

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.manifold import TSNE

    # Prepare features from CA coordinates
    ca_indices = [i for i, n in enumerate(mol.atom_name)
                  if ''.join(n).strip() == "CA"]
    features = traj.coords[:, ca_indices, :].reshape(traj.nframe, -1).astype(np.float32)

    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(features)
    test("StandardScaler", scaled_features.shape == features.shape)

    # Quick t-SNE with minimum iterations for testing
    n_samples = min(500, traj.nframe)  # Limit samples for speed
    tsne = TSNE(n_components=2, random_state=42,
                perplexity=min(30.0, n_samples - 1), max_iter=250)
    embedding = tsne.fit_transform(scaled_features[:n_samples])
    test("t-SNE embedding", embedding.shape == (n_samples, 2))

    print(f"\n  Features shape: {features.shape}")
    print(f"  t-SNE embedding shape: {embedding.shape}")

except ImportError:
    print("  [SKIP] scikit-learn not installed")


# =============================================================================
# PyTorch Integration Tests
# =============================================================================
section("14. PyTorch Integration")

try:
    import torch

    # Convert trajectory to PyTorch tensor
    ca_indices = [i for i, n in enumerate(mol.atom_name)
                  if ''.join(n).strip() == "CA"]
    features = traj.coords[:, ca_indices, :].reshape(traj.nframe, -1).astype(np.float32)

    tensor = torch.tensor(features, dtype=torch.float32)
    test("PyTorch tensor creation", tensor.shape == (traj.nframe, len(ca_indices) * 3))

    print(f"\n  PyTorch tensor shape: {tensor.shape}")
    print(f"  PyTorch version: {torch.__version__}")

except ImportError:
    print("  [SKIP] PyTorch not installed")


# =============================================================================
# Error Handling Tests
# =============================================================================
section("15. Error Handling")

# Test custom exceptions exist and work
try:
    raise GenesisValidationError("test error")
except GenesisError:
    test("GenesisValidationError inherits GenesisError", True)

try:
    raise GenesisFortranError("test error", code=1, stderr_output="test stderr")
except GenesisError as e:
    test("GenesisFortranError has code attribute", hasattr(e, 'code'))
    test("GenesisFortranError has stderr_output attribute", hasattr(e, 'stderr_output'))


# =============================================================================
# Summary
# =============================================================================
section("Test Summary")

total = tests_passed + tests_failed
print(f"\n  Tests passed: {tests_passed}/{total}")
print(f"  Tests failed: {tests_failed}/{total}")

if tests_failed == 0:
    print("\n" + "="*60)
    print("  ALL TESTS PASSED!")
    print("="*60)
    sys.exit(0)
else:
    print("\n" + "="*60)
    print(f"  {tests_failed} TEST(S) FAILED")
    print("="*60)
    sys.exit(1)
