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


def test_rmsd_analysis():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    trajs, subset_mol = genesis_exe.crd_convert(
        mol,
        trj_files=["BPTI_run.dcd"],
        trj_format="DCD",
        trj_type="COOR+BOX",
        selection="all",
    )

    _ = subset_mol

    for t in trajs:
        d = genesis_exe.rmsd_analysis(
            mol, t,
            selection_group=["sid:BPTI and an:CA"],
            fitting_method="TR+ROT",
            fitting_atom=1,
            check_only=False,
            analysis_atom=1,
        )
        print(d.rmsd, flush=True)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_rmsd_analysis()


if __name__ == "__main__":
    main()
