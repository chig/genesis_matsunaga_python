import os
import ctypes
import pathlib
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
import genesis_exe


def main():
    if os.path.exists("out"):
        os.remove("out")
    # test_mbar_analysis_umbrella_block()


if __name__ == "__main__":
    main()
