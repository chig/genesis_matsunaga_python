import ctypes


class SMolecule:
    num_deg_freedom: int
    num_atoms: int
    num_bonds: int
    num_enm_bonds: int
    num_angles: int
    num_dihedrals: int
    num_impropers: int
    num_cmaps: int
    num_residues: int
    num_molecules: int
    num_segments: int
    shift_origin: bool
    special_hydrogen: bool
    total_charge: float
    atom_no: list[int]
    segment_name: list[str]
    segment_no: list[int]
    residue_no: list[int]
    residue_c_no: list[int]
    residue_name: list[str]
    atom_name: list[str]
    atom_cls_name: list[str]
    atom_cls_no: list[int]
    charge: list[float]
    mass: list[float]
    inv_mass: list[float]
    imove: list[int]
    stokes_radius: list[float]
    inv_stokes_radius: list[float]
    chain_id: list[str]
    atom_coord: list[list[float]]
    atom_occupancy: list[float]
    atom_temp_factor: list[float]
    atom_velocity: list[list[float]]
    light_atom_name: list[bool]
    light_atom_mass: list[bool]
    molecule_no: list[int]
    bond_list: list[list[int]]
    enm_list: list[list[int]]
    angl_list: list[list[int]]
    dihe_list: list[list[int]]
    impr_list: list[list[int]]
    cmap_list: list[list[int]]
    molecule_atom_no: list[int]
    molecule_mass: list[float]
    molecule_name: list[str]
    atom_refcoord: list[list[float]]
    atom_fitcoord: list[list[float]]
    num_pc_modes: int
    pc_mode: list[float]
    fep_topology: int
    num_hbonds_singleA: int
    num_hbonds_singleB: int
    num_atoms_fep: list[int]
    num_bonds_fep: list[int]
    num_angles_fep: list[int]
    num_dihedrals_fep: list[int]
    num_impropers_fep: list[int]
    num_cmaps_fep: list[int]
    bond_list_fep: list[list[list[int]]]
    angl_list_fep: list[list[list[int]]]
    dihe_list_fep: list[list[list[int]]]
    impr_list_fep: list[list[list[int]]]
    cmap_list_fep: list[list[list[int]]]
    id_singleA: list[int]
    id_singleB: list[int]
    fepgrp: list[int]
    fepgrp_bond: list[list[int]]
    fepgrp_angl: list[list[list[int]]]
    fepgrp_dihe: list[list[list[list[int]]]]
    fepgrp_cmap: list[int]


class SMoleculeC(ctypes.Structure):
     _fields_ = [("num_deg_freedom", ctypes.c_int),
                 ("num_atoms", ctypes.c_int),
                 ("num_bonds", ctypes.c_int),
                 ("num_enm_bonds", ctypes.c_int),
                 ("num_angles", ctypes.c_int),
                 ("num_dihedrals", ctypes.c_int),
                 ("num_impropers", ctypes.c_int),
                 ("num_cmaps", ctypes.c_int),
                 ("num_residues", ctypes.c_int),
                 ("num_molecules", ctypes.c_int),
                 ("num_segments", ctypes.c_int),
                 ("shift_origin", ctypes.c_bool),
                 ("special_hydrogen", ctypes.c_bool),
                 ("total_charge", ctypes.c_double),
                 ("atom_no", ctypes.c_void_p),
                 ("segment_name", ctypes.c_void_p),
                 ("segment_no", ctypes.c_void_p),
                 ("residue_no", ctypes.c_void_p),
                 ("residue_c_no", ctypes.c_void_p),
                 ("residue_name", ctypes.c_void_p),
                 ("atom_name", ctypes.c_void_p),
                 ("atom_cls_name", ctypes.c_void_p),
                 ("atom_cls_no", ctypes.c_void_p),
                 ("charge", ctypes.c_void_p),
                 ("mass", ctypes.c_void_p),
                 ("inv_mass", ctypes.c_void_p),
                 ("imove", ctypes.c_void_p),
                 ("stokes_radius", ctypes.c_void_p),
                 ("inv_stokes_radius", ctypes.c_void_p),
                 ("chain_id", ctypes.c_void_p),
                 ("atom_coord", ctypes.c_void_p),
                 ("atom_occupancy", ctypes.c_void_p),
                 ("atom_temp_factor", ctypes.c_void_p),
                 ("atom_velocity", ctypes.c_void_p),
                 ("light_atom_name", ctypes.c_void_p),
                 ("light_atom_mass", ctypes.c_void_p),
                 ("molecule_no", ctypes.c_void_p),
                 ("bond_list", ctypes.c_void_p),
                 ("enm_list", ctypes.c_void_p),
                 ("angl_list", ctypes.c_void_p),
                 ("dihe_list", ctypes.c_void_p),
                 ("impr_list", ctypes.c_void_p),
                 ("cmap_list", ctypes.c_void_p),
                 ("molecule_atom_no", ctypes.c_void_p),
                 ("molecule_mass", ctypes.c_void_p),
                 ("molecule_name", ctypes.c_void_p),
                 ("atom_refcoord", ctypes.c_void_p),
                 ("atom_fitcoord", ctypes.c_void_p),
                 ("num_pc_modes", ctypes.c_int),
                 ("pc_mode", ctypes.c_void_p),
                 ("fep_topology", ctypes.c_int),
                 ("num_hbonds_singleA", ctypes.c_int),
                 ("num_hbonds_singleB", ctypes.c_int),
                 ("num_atoms_fep", ctypes.c_void_p),
                 ("num_bonds_fep", ctypes.c_void_p),
                 ("num_angles_fep", ctypes.c_void_p),
                 ("num_dihedrals_fep", ctypes.c_void_p),
                 ("num_impropers_fep", ctypes.c_void_p),
                 ("num_cmaps_fep", ctypes.c_void_p),
                 ("bond_list_fep", ctypes.c_void_p),
                 ("angl_list_fep", ctypes.c_void_p),
                 ("dihe_list_fep", ctypes.c_void_p),
                 ("impr_list_fep", ctypes.c_void_p),
                 ("cmap_list_fep", ctypes.c_void_p),
                 ("id_singleA", ctypes.c_void_p),
                 ("id_singleB", ctypes.c_void_p),
                 ("fepgrp", ctypes.c_void_p),
                 ("fepgrp_bond", ctypes.c_void_p),
                 ("fepgrp_angl", ctypes.c_void_p),
                 ("fepgrp_dihe", ctypes.c_void_p),
                 ("fepgrp_cmap", ctypes.c_void_p),
                 ("nbnd_fep_max", ctypes.c_int),
                 ("nangl_fep_max", ctypes.c_int),
                 ("ndihe_fep_max", ctypes.c_int),
                 ("nimpr_fep_max", ctypes.c_int),
                 ("ncmap_fep_max", ctypes.c_int),
                 ("size_id_singleA", ctypes.c_int),
                 ("size_id_singleB", ctypes.c_int),
                 ("size_fepgrp", ctypes.c_int),
                 ]


def c2py_s_molecule(src: SMoleculeC) -> SMolecule:
    dst = SMolecule()
    dst.num_deg_freedom = src.num_deg_freedom
    dst.num_atoms     = src.num_atoms
    dst.num_bonds     = src.num_bonds
    dst.num_enm_bonds = src.num_enm_bonds
    dst.num_angles    = src.num_angles
    dst.num_dihedrals = src.num_dihedrals
    dst.num_impropers = src.num_impropers
    dst.num_cmaps     = src.num_cmaps
    dst.num_residues  = src.num_residues
    dst.num_molecules = src.num_molecules
    dst.num_segments  = src.num_segments
    dst.shift_origin  = src.shift_origin
    dst.special_hydrogen = src.special_hydrogen
    dst.total_charge  = src.total_charge
    dst.atom_no = c2py_int_array(src.atom_no, dst.num_atoms)
    dst.segment_name = c2py_str_array(src.segment_name, dst.num_atoms, 4)
    dst.segment_no = c2py_int_array(src.segment_no, dst.num_atoms)
    dst.residue_no = c2py_int_array(src.residue_no, dst.num_atoms)
    dst.residue_c_no = c2py_int_array(src.residue_c_no, dst.num_atoms)
    dst.residue_name = c2py_str_array(src.residue_name, dst.num_atoms, 6)
    dst.atom_name = c2py_str_array(src.atom_name, dst.num_atoms, 4)
    dst.atom_cls_name = c2py_str_array(src.atom_cls_name, dst.num_atoms, 6)
    dst.atom_cls_no = c2py_int_array(src.atom_cls_no, dst.num_atoms)
    dst.charge = c2py_double_array(src.charge, dst.num_atoms)
    dst.mass = c2py_double_array(src.mass, dst.num_atoms)
    dst.inv_mass = c2py_double_array(src.inv_mass, dst.num_atoms)
    dst.imove = c2py_int_array(src.imove, dst.num_atoms)
    dst.stokes_radius = c2py_double_array(src.stokes_radius, dst.num_atoms)
    dst.inv_stokes_radius = c2py_double_array(src.inv_stokes_radius, dst.num_atoms)
    return dst


def c2py_int_array(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_int))
        return list(int(ptr[i]) for i in range(0, size))
    else:
        return ()


def c2py_double_array(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_double))
        return list(float(ptr[i]) for i in range(0, size))
    else:
        return ()


def c2py_str_array(src: ctypes.c_void_p, size_array: int, size_str: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_char))
        return list(c2py_str(ptr[i * size_str], size_str) for i in range(0, size_array))
    else:
        return ()


def c2py_str(src: ctypes.c_void_p, size: int):
    if src:
        ptr = ctypes.cast(src, ctypes.POINTER(ctypes.c_char))
        return ''.join(str(ptr[i].decode('utf8')) for i in range(0, size))
    return ' ' * size
