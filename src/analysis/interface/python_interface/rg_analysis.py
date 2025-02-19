import os
import ctypes
import pathlib
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
import genesis_exe


def test_rg_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    rg_analysis_ctrl_path = pathlib.Path("test_rg_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                rg = genesis_exe.rg_analysis(mol, t, 1, rg_analysis_ctrl_path)
                print(rg)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_rg_analysis()


if __name__ == "__main__":
    main()
