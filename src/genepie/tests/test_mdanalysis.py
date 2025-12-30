# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import unittest
from ..s_trajectories import STrajectories
import MDAnalysis as mda
from ..custom_test_case import CustomTestCase


class TestMdAnalysis(CustomTestCase):

    def test_from_mdanalysis_universe(self):
        uni = mda.Universe(self.PDB_PATH, self.TRJ_PATH)
        uni.atoms.guess_bonds()
        trj, mol = STrajectories.from_mdanalysis_universe(uni)
        gtrajs, gmol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        # New API returns List[STrajectories], not context manager
        self.assertAlmostEqual(gtrajs[0], trj)

    def test_to_mdanalysis_universe(self):
        gtrajs, gmol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        # New API returns List[STrajectories], iterate directly
        for t in gtrajs:
            uni = t.to_mdanalysis_universe(gmol)
            gt, gm = STrajectories.from_mdanalysis_universe(uni)
            self.assertAlmostEqual(t, gt)


if __name__ == "__main__":
    unittest.main()
