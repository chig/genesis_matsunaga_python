# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_avecrd_analysis():
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(BPTI_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
    )

    _ = subset_mol

    for t in trajs:
        d = genesis_exe.avecrd_analysis(
            mol, t,
            selection_group=["an:CA"],
            fitting_method="TR+ROT",
            fitting_atom=1,
            check_only=False,
            num_iterations=5,
            analysis_atom=1,
        )
        print(d.pdb)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    try:
        test_avecrd_analysis()
        print("\n✓ test_avecrd: PASSED")
    except Exception as e:
        print(f"\n✗ test_avecrd: FAILED - {e}")
        raise


if __name__ == "__main__":
    main()
