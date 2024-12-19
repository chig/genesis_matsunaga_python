import ctypes
import os
from s_molecule import SMoleculeC, c2py_s_molecule


def define_prototypes(lib: ctypes.CDLL):
    # 関数のプロトタイプを定義
    lib.define_molecule_from_pdb.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(SMoleculeC),
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_void_p),
            ctypes.POINTER(ctypes.c_void_p)]
    lib.define_molecule_from_pdb.restype = None

    lib.deallocate_s_molecule_c.argtypes = [
            ctypes.POINTER(SMoleculeC)]
    lib.deallocate_s_molecule_c.restype = None

def test():
    # ライブラリをロード
    lib_name = 'libpython_interface.so'
    lib_dir = os.path.join(os.path.dirname(__file__), '.libs')
    lib_path = os.path.join(lib_dir, lib_name)

    if not os.path.exists(lib_path):
        raise FileNotFoundError(f"Library file {lib_name} not found in {lib_dir}")

    lib = ctypes.CDLL(lib_path)

    define_prototypes(lib)

    # 関数を呼び出す
    pdb_filename = b"molecule.pdb"
    mol_c = SMoleculeC()
    num_atoms = ctypes.c_int()
    atom_names_ptr = ctypes.c_void_p()
    atom_coords_ptr = ctypes.c_void_p()

    lib.define_molecule_from_pdb(
            pdb_filename,
            ctypes.byref(mol_c),
            ctypes.byref(num_atoms),
            ctypes.byref(atom_names_ptr),
            ctypes.byref(atom_coords_ptr))
    mol = c2py_s_molecule(mol_c)

    # 結果を処理する
    print("num_atoms = ", mol.num_atoms)
    for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
        print(mol.atom_no[i], mol.segment_name[i], mol.atom_name[i])

    # メモリを解放する
    lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


if __name__ == "__main__":
    test()

