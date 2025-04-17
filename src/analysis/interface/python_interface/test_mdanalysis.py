import pathlib
import unittest
import numpy as np
from s_trajectories import STrajectories, STrajectoriesArray
import MDAnalysis as mda
from ctrl_files import TrajectoryParameters
from s_molecule import SMolecule
import genesis_exe


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


class TestMdAnalysis(unittest.TestCase):

    def test_from_mdanalysis_universe(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        trj_path = pathlib.Path("BPTI_run.dcd")
        uni = mda.Universe(pdb_path, trj_path)
        trj, mol = STrajectories.from_mdanalysis_universe(uni)
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

    def test_to_mdanalysis_universe(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        psf_path = pathlib.Path("BPTI_ionize.psf")
        mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
        with genesis_exe.crd_convert(
                mol,
                traj_params = [
                    TrajectoryParameters(
                        trjfile = "BPTI_run.dcd",
                        md_step = 10,
                        mdout_period = 1,
                        ana_period = 1,
                        repeat = 1,
                        ),
                    ],
                trj_format = "DCD",
                trj_type = "COOR+BOX",
                trj_natom = 0,
                selection_group = ["an:CA", ],
                fitting_method = "TR+ROT",
                fitting_atom = 1,
                check_only = False,
                centering = True,
                centering_atom = 1,
                center_coord = (0.0, 0.0, 0.0),
                pbc_correct = "molecule",
                rename_res = ["HSE HIS", "HSE HIS"],
                ) as trajs:
            for t in trajs:
                _ = t.to_mdanalysis_universe(mol)
                # uni.atoms.write("test.pdb")

if __name__ == "__main__":
    unittest.main()
