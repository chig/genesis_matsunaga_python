# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent))
    __package__ = "genepie"
# --------------------------------------------
import os
import pathlib
from .s_molecule import SMolecule
from . import genesis_exe


def test_hb_analysis_Count_snap():
    pdb_path = pathlib.Path("RALP_DPPC_run.pdb")
    psf_path = pathlib.Path("RALP_DPPC.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=["RALP_DPPC_run.dcd"],
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
            output_type="Count_Snap",
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
    test_hb_analysis_Count_snap()


if __name__ == "__main__":
    main()
