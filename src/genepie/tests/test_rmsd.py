# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
import subprocess
import sys
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_rmsd_no_fitting():
    """Test RMSD analysis without fitting (direct comparison to reference)."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol = genesis_exe.crd_convert(
            mol,
            trj_files=[str(BPTI_DCD)],
            trj_format="DCD",
            trj_type="COOR+BOX",
            selection="all",
    )
    _ = subset_mol

    for t in trajs:
        # Test RMSD without fitting (fitting_selection=None)
        result = genesis_exe.rmsd_analysis(
            mol, t,
            analysis_selection="an:CA",
            fitting_selection=None,  # No fitting
            ana_period=1,
            mass_weighted=False,
        )

        # Validate results
        assert result.rmsd is not None, "RMSD result should not be None"
        assert len(result.rmsd) > 0, "RMSD result should have at least one frame"
        assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
        assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

        print(f"RMSD no fitting (n={len(result.rmsd)}): "
              f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")


def test_rmsd_with_fitting():
    """Test RMSD analysis with TR+ROT fitting."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol = genesis_exe.crd_convert(
            mol,
            trj_files=[str(BPTI_DCD)],
            trj_format="DCD",
            trj_type="COOR+BOX",
            selection="all",
    )
    _ = subset_mol

    for t in trajs:
        # Test RMSD with TR+ROT fitting
        result = genesis_exe.rmsd_analysis(
            mol, t,
            analysis_selection="sid:BPTI and an:CA",
            fitting_selection="sid:BPTI and an:CA",  # Enable fitting
            fitting_method="TR+ROT",
            ana_period=1,
            mass_weighted=False,
        )

        # Validate results
        assert result.rmsd is not None, "RMSD with fitting should not be None"
        assert len(result.rmsd) > 0, "RMSD should have at least one frame"
        assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
        assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

        print(f"RMSD with fitting (n={len(result.rmsd)}): "
              f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")


def test_rmsd_lazy_unified_api():
    """Test RMSD analysis using unified lazy API (crd_convert(lazy=True) + rmsd_analysis)."""
    import numpy as np

    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Create lazy trajectory via crd_convert(lazy=True)
    lazy_trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )

    lazy_traj = lazy_trajs[0]
    assert lazy_traj.is_lazy, "Trajectory should be lazy"
    assert lazy_traj.lazy_dcd_file is not None, "lazy_dcd_file should be set"
    assert lazy_traj.lazy_trj_type == 2, "lazy_trj_type should be COOR+BOX (2)"

    # Call rmsd_analysis with lazy trajectory (unified API)
    result = genesis_exe.rmsd_analysis(
        subset_mol,
        lazy_traj,
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=1,
        mass_weighted=False,
    )

    # Validate results
    assert result.rmsd is not None, "RMSD result should not be None"
    assert len(result.rmsd) > 0, "RMSD result should have at least one frame"
    assert all(r >= 0 for r in result.rmsd), "RMSD values should be non-negative"
    assert all(r < 50.0 for r in result.rmsd), "RMSD values should be reasonable (< 50 A)"

    print(f"RMSD lazy unified API (n={len(result.rmsd)}): "
          f"min={min(result.rmsd):.5f}, max={max(result.rmsd):.5f}")


def test_rmsd_lazy_vs_memory():
    """Compare lazy and memory-based RMSD to ensure identical results."""
    import numpy as np

    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)

    # Memory-based RMSD
    trajs_mem, mol_mem = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=False,
    )
    result_mem = genesis_exe.rmsd_analysis(
        mol_mem,
        trajs_mem[0],
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=1,
        mass_weighted=False,
    )

    # Lazy-based RMSD
    trajs_lazy, mol_lazy = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        lazy=True,
    )
    result_lazy = genesis_exe.rmsd_analysis(
        mol_lazy,
        trajs_lazy[0],
        analysis_selection="an:CA",
        fitting_selection=None,
        ana_period=1,
        mass_weighted=False,
    )

    # Compare results
    assert len(result_mem.rmsd) == len(result_lazy.rmsd), \
        f"Frame count mismatch: memory={len(result_mem.rmsd)}, lazy={len(result_lazy.rmsd)}"

    for i, (mem_val, lazy_val) in enumerate(zip(result_mem.rmsd, result_lazy.rmsd)):
        assert np.isclose(mem_val, lazy_val, rtol=1e-4, atol=1e-6), \
            f"Frame {i}: memory={mem_val}, lazy={lazy_val}"

    print(f"Memory vs Lazy comparison passed: {len(result_mem.rmsd)} frames")
    print(f"  Memory: min={min(result_mem.rmsd):.5f}, max={max(result_mem.rmsd):.5f}")
    print(f"  Lazy:   min={min(result_lazy.rmsd):.5f}, max={max(result_lazy.rmsd):.5f}")


def _run_test_in_subprocess(test_name: str) -> bool:
    """Run a single test function in isolated subprocess to avoid Fortran state issues."""
    code = f'''
import sys
# Handle package imports when run as subprocess
if __name__ == "__main__":
    import pathlib
    pkg_dir = pathlib.Path("{__file__}").resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))

from genepie.tests.test_rmsd import {test_name}
{test_name}()
'''
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=120
    )

    if result.returncode == 0:
        # Print stdout (test output)
        if result.stdout:
            print(result.stdout, end='')
        return True
    else:
        print(f"stderr: {result.stderr}" if result.stderr else "")
        return False


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")

    # Run each test in separate subprocess to avoid Fortran global state accumulation
    # This prevents crashes on macOS ARM64 when multiple tests run in same process
    tests = [
        "test_rmsd_no_fitting",
        "test_rmsd_with_fitting",
        "test_rmsd_lazy_unified_api",
        "test_rmsd_lazy_vs_memory",
    ]

    failed = []
    for test_name in tests:
        if _run_test_in_subprocess(test_name):
            print(f"\n{test_name}: PASSED")
        else:
            print(f"\n{test_name}: FAILED")
            failed.append(test_name)

    if failed:
        raise RuntimeError(f"Tests failed: {', '.join(failed)}")

    print("\nAll RMSD tests passed!")


if __name__ == "__main__":
    main()
