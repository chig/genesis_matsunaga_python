# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
from .conftest import RALP_PDB, RALP_PSF, RALP_DCD
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_hb_analysis_Count_atom():
    mol = SMolecule.from_file(pdb=RALP_PDB, psf=RALP_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=[str(RALP_DCD)],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
        centering=True,
        centering_selection="all",
        center_coord=(0.0, 0.0, 0.0),
        rename_res=["HSE HIS", "HSD HIS"],
    )

    _ = subset_mol

    for t in trajs:
        d = genesis_exe.hb_analysis(
            mol, t,
            selection_group=["sid:PROA",
                             "resname:DPPC & (an:O11 | an:O12 | an:O13 | an:O14)"],
            check_only=False,
            output_type="Count_atom",
            solvent_list="DPPC",
            analysis_atom=1,
            target_atom=2,
            boundary_type="PBC",
            hb_distance=3.4,
            dha_angle=120.0,
            hda_angle=30.0,
        )
        print(d, flush=True)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_hb_analysis_Count_atom()


if __name__ == "__main__":
    main()
