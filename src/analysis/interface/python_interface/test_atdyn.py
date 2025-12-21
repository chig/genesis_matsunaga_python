#!/usr/bin/env python3
"""
Test script for atdyn Python interface.
Uses regression test data files.
"""

import sys
import os

# Add path for development
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from python_interface import genesis_exe


def test_atdyn_min_glycam():
    """Test energy minimization with AMBER/glycam system."""
    print("=" * 60)
    print("Testing run_atdyn_min (glycam system)")
    print("=" * 60)

    # Get path to test files
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
    test_dir = os.path.join(base_dir, 'tests', 'regression_test', 'build', 'glycam')

    prmtopfile = os.path.join(test_dir, 'glycam.top')
    ambcrdfile = os.path.join(test_dir, 'glycam.rst')
    rstfile = os.path.join(test_dir, 'rst')

    print(f"  prmtopfile: {prmtopfile}")
    print(f"  ambcrdfile: {ambcrdfile}")
    print(f"  rstfile: {rstfile}")

    if not os.path.exists(prmtopfile):
        print(f"Test files not found: {prmtopfile}")
        print("SKIPPED")
        return True

    try:
        # Run minimization via Python interface (same params as reference test)
        result = genesis_exe.run_atdyn_min(
            prmtopfile=prmtopfile,
            ambcrdfile=ambcrdfile,
            rstfile=rstfile,
            forcefield="AMBER",
            electrostatic="PME",
            switchdist=12.0,
            cutoffdist=12.0,
            pairlistdist=14.0,
            pme_alpha=0.34,
            pme_ngrid_x=64,
            pme_ngrid_y=64,
            pme_ngrid_z=64,
            pme_nspline=4,
            dispersion_corr="epress",
            method="SD",
            nsteps=20,
            eneout_period=2,
            nbupdate_period=4,
            rigid_bond=False,
            boundary_type="PBC",
            box_size_x=69.5294360,
            box_size_y=68.0597930,
            box_size_z=56.2256950,
        )

        print(f"Minimization completed successfully!")
        print(f"  Number of atoms: {result.final_coords.shape[1]}")
        print(f"  Energy shape: {result.energies.shape}")
        print(f"  Final energy (total): {result.energies[0, 0]:.4f}")
        print(f"  Final gradient: {result.final_gradient:.6f}")
        print(f"  Converged: {result.converged}")

        # Check that energy is reasonable (should be negative for this system)
        if result.energies[0, 0] > 0:
            print("ERROR: Total energy should be negative for this system")
            return False

        print("PASSED")
        return True

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    results = []

    results.append(("run_atdyn_min (glycam)", test_atdyn_min_glycam()))

    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)

    all_passed = True
    for name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    sys.exit(0 if all_passed else 1)
