import pathlib
import unittest
import mdtraj
import numpy as np
from ctrl_files import TrajectoryParameters
import genesis_exe
from s_molecule import SMolecule
from s_trajectories import STrajectories, STrajectoriesArray


def create_traj_by_genesis(pdb_path, dcd_path) \
        -> tuple[STrajectoriesArray, SMolecule]:
    mol = SMolecule.from_file(pdb=pdb_path)
    trajs = genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = dcd_path,
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
    return (trajs, mol)


class TestMDTraj(unittest.TestCase):

    def test_from_mdtraj_trajectory(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        trj_path = pathlib.Path("BPTI_run.dcd")
        mdt = mdtraj.load(trj_path, top=pdb_path)

        trj, mol = STrajectories.from_mdtraj_trajectory(mdt)
        gtrajs, gmol = create_traj_by_genesis(pdb_path, trj_path)
        with trj, gtrajs:
            self.assertEqual(trj.natom, gtrajs[0].natom)
            self.assertEqual(trj.nframe, gtrajs[0].nframe)
            self.assertTrue(np.allclose(trj.coords, gtrajs[0].coords,
                                        rtol=1e-4, atol=1e-4))
            self.assertTrue(np.allclose(trj.pbc_boxes, gtrajs[0].pbc_boxes,
                                        rtol=1e-4, atol=1e-4))

            self.assertEqual(mol.num_deg_freedom, gmol.num_deg_freedom)
            self.assertEqual(mol.num_atoms, gmol.num_atoms)
            self.assertEqual(mol.num_bonds, gmol.num_bonds)
            self.assertEqual(mol.num_enm_bonds, gmol.num_enm_bonds)
            self.assertEqual(mol.num_angles, gmol.num_angles)
            self.assertEqual(mol.num_dihedrals, gmol.num_dihedrals)
            self.assertEqual(mol.num_impropers, gmol.num_impropers)

    def test_to_mdtraj_trajectory(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        trj_path = pathlib.Path("BPTI_run.dcd")
        strajs, smol = create_traj_by_genesis(pdb_path, trj_path)
        with strajs:
            for t in strajs:
                _ = t.to_mdtraj_trajectory(smol)

if __name__ == "__main__":
    unittest.main()
