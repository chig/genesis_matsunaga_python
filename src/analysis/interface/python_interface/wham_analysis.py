import os
import ctypes
import pathlib
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
import genesis_exe


def test_wham_analysis():
    ctrl_path = pathlib.Path("test_wham_analysis_inp")
    genesis_exe.wham_analysis(ctrl_path)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_wham_analysis()


if __name__ == "__main__":
    main()
