#!/usr/bin/env python
"""
Comprehensive integration test script for GENESIS Python interface (genepie).

This script tests the main functionality demonstrated in demo.ipynb,
including molecular loading, trajectory processing, and various analyses.

Usage:
    python -m genepie.tests.test_integration

Requirements:
    - chignolin test data (download with: python -m genepie.tests.download_test_data)
    - genepie package installed
"""
# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------

import os
import sys
import tempfile
import numpy as np
import warnings

from .conftest import CHIGNOLIN_PDB, CHIGNOLIN_PSF, CHIGNOLIN_DCD

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


def run_tests():
    """Main test runner function."""
    global tests_passed, tests_failed
    tests_passed = 0
    tests_failed = 0

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
        return 1

    # Check input files exist
    test("chignolin.pdb exists", CHIGNOLIN_PDB.exists())
    test("chignolin.psf exists", CHIGNOLIN_PSF.exists())
    test("chignolin.dcd exists", CHIGNOLIN_DCD.exists())

    if not all([CHIGNOLIN_PDB.exists(), CHIGNOLIN_PSF.exists(), CHIGNOLIN_DCD.exists()]):
        print("\nFATAL: Required input files not found.")
        print("       Run: python -m genepie.tests.download_test_data")
        return 1


    # =============================================================================
    # SMolecule Tests
    # =============================================================================
    section("2. SMolecule - Loading Molecular Information")

    mol_allatom = SMolecule.from_file(pdb=CHIGNOLIN_PDB, psf=CHIGNOLIN_PSF)
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

    # Test 1: CA atom selection
    trajs_array_ca, mol_ca = genesis_exe.crd_convert(
        molecule=mol_allatom,
        trj_files=[str(CHIGNOLIN_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="an:CA",
        fitting_selection="an:CA",
        fitting_method="TR+ROT",
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
        trj_files=[str(CHIGNOLIN_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="segid:PROA",
        fitting_selection="segid:PROA",
        fitting_method="TR+ROT",
        pbc_correct="NO"
    )
    traj = trajs_array[0]
    test("PROA selection - trajectory loaded", traj.nframe > 0)
    test("PROA selection - natom matches mol", traj.natom == mol.num_atoms)

    print(f"  PROA trajectory: {traj.nframe} frames, {traj.natom} atoms")

    # Test 3: Heavy atoms selection
    trajs_heavy, mol_heavy = genesis_exe.crd_convert(
        molecule=mol_allatom,
        trj_files=[str(CHIGNOLIN_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="heavy",
        fitting_selection="heavy",
        fitting_method="TR+ROT",
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

    # Get CA atom indices for residues 1-4
    ca_indices = genesis_exe.selection(mol, "an:CA")
    # Take first 4 CA atoms for testing
    if len(ca_indices) >= 4:
        # Create distance pairs (1-2, 2-3), angle triplets (1-2-3), torsion quads (1-2-3-4)
        distance_pairs = np.array([[ca_indices[0], ca_indices[1]]], dtype=np.int32)
        angle_triplets = np.array([[ca_indices[0], ca_indices[1], ca_indices[2]]], dtype=np.int32)
        torsion_quadruplets = np.array([[ca_indices[0], ca_indices[1], ca_indices[2], ca_indices[3]]], dtype=np.int32)

        analysis_results = genesis_exe.trj_analysis(
            trajs=traj,
            distance_pairs=distance_pairs,
            angle_triplets=angle_triplets,
            torsion_quadruplets=torsion_quadruplets
        )

        test("distance result shape", analysis_results.distance.shape[0] == traj.nframe)
        test("angle result shape", analysis_results.angle.shape[0] == traj.nframe)
        test("torsion result shape", analysis_results.torsion.shape[0] == traj.nframe)
        test("distance values reasonable", 2.0 < analysis_results.distance.mean() < 6.0,
             f"mean = {analysis_results.distance.mean():.3f}")

        print(f"\n  Distance mean: {analysis_results.distance.mean():.3f} A")
        print(f"  Angle mean: {analysis_results.angle.mean():.3f} deg")
        print(f"  Torsion mean: {analysis_results.torsion.mean():.3f} deg")
    else:
        print("  [SKIP] Not enough CA atoms for trj_analysis test")


    # =============================================================================
    # rg_analysis Tests
    # =============================================================================
    section("5. rg_analysis - Radius of Gyration")

    rg_result = genesis_exe.rg_analysis(
        molecule=mol, trajs=traj,
        analysis_selection="an:CA",
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
        analysis_selection="an:CA",
        fitting_selection="an:CA",
        fitting_method="TR+ROT",
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

    # Create msd_data with time column (diffusion_analysis expects column 0 = time)
    nsteps = msd_result.msd.shape[0]
    time_col = np.arange(nsteps, dtype=np.float64).reshape(-1, 1)
    msd_data = np.hstack([time_col, msd_result.msd])

    # Use start_step as integer (20% of data)
    start_step = max(1, int(nsteps * 0.2))
    diffusion_result = genesis_exe.diffusion_analysis(
        msd_data=msd_data,
        time_step=1.0,
        start_step=start_step
    )

    test("diffusion result not None", diffusion_result is not None)
    test("diffusion coefficients available",
         diffusion_result.diffusion_coefficients is not None)

    print(f"\n  Diffusion coeffs: {diffusion_result.diffusion_coefficients}")


    # =============================================================================
    # wham_analysis Tests
    # =============================================================================
    section("10. wham_analysis - Free Energy")

    # Create test CV data in temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        cv_dir = os.path.join(tmpdir, "cv_wham_data")
        os.makedirs(cv_dir)
        num_cv_frames = 200
        rows = np.arange(1, num_cv_frames + 1)
        np.random.seed(42)  # For reproducibility
        for i, loc in enumerate([1.0, 2.0, 3.0], 1):
            values = np.random.normal(loc=loc, scale=0.3, size=num_cv_frames)
            data = np.column_stack((rows, values))
            np.savetxt(os.path.join(cv_dir, f"window_{i}.dat"), data, fmt="%d %.5f")

        cv_file_template = os.path.join(cv_dir, "window_{:d}.dat")
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
        return 0
    else:
        print("\n" + "="*60)
        print(f"  {tests_failed} TEST(S) FAILED")
        print("="*60)
        return 1


def main():
    sys.exit(run_tests())


if __name__ == "__main__":
    main()
