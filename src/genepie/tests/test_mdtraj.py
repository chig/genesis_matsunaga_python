# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import unittest
import mdtraj
from ..s_trajectories import STrajectories
from ..custom_test_case import CustomTestCase


class TestMDTraj(CustomTestCase):

    def test_from_mdtraj_trajectory(self):
        mdt = mdtraj.load(self.TRJ_PATH, top=self.PDB_PATH)

        trj, mol = STrajectories.from_mdtraj_trajectory(mdt)
        gtrajs, gmol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH)
        # New API returns List[STrajectories], not context manager
        self.assertAlmostEqual(gtrajs[0], trj)

    def test_to_mdtraj_trajectory(self):
        strajs, smol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        # New API returns List[STrajectories], iterate directly
        for t in strajs:
            mdt = t.to_mdtraj_trajectory(smol)
            gt, gm = STrajectories.from_mdtraj_trajectory(mdt)
            self.assertAlmostEqual(t, gt)


if __name__ == "__main__":
    unittest.main()
