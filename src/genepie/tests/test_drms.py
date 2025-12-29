# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..ctrl_files import TrajectoryParameters
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_drms_analysis():
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol =  genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = str(BPTI_DCD),
                    md_step = 10,
                    mdout_period = 1,
                    ana_period = 1,
                    repeat = 1,
                    ),
                ],
            trj_format = "DCD",
            trj_type = "COOR+BOX",
            trj_natom = 0,
            selection_group = ["all", ],
            fitting_method = "NO",
            fitting_atom = 1,
            check_only = False,
            pbc_correct = "NO",
    ) 

    _ = subset_mol

    try: 
        for t in trajs:
            drms, = genesis_exe.drms_analysis(
                    mol, t,
                    selection_group = ["an: CA", ],
                    check_only = False,
                    contact_groups = 1,
                    ignore_hydrogen  = False,
                    two_states       = False,
                    avoid_bonding    = True,
                    exclude_residues = 4,
                    minimum_distance = 1.0,
                    maximum_distance = 6.0,
                    pbc_correct      = False,
                    verbose          = True,
                    )
            print(drms)
    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    try:
        test_drms_analysis()
        print("\n✓ test_drms: PASSED")
    except Exception as e:
        print(f"\n✗ test_drms: FAILED - {e}")
        raise


if __name__ == "__main__":
    main()
