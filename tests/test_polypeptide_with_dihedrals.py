import unittest
from itertools import islice
from pathlib import Path

from Bio.PDB import PDBParser
from numpy.ma.testutils import assert_equal

from polypeptide_with_dihedrals import PolypeptideWithDihedrals


class TestPolypeptideWithDihedrals(unittest.TestCase):
    def setUp(self):
        base_dir = Path(__file__).resolve().parent.parent  # go up from tests/ to project root
        pdb_path = base_dir / "examples" / "9lil.pdb"
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("example", pdb_path)
        self.polypeptide = self.structure[0]['A'].child_list[11:20]
        self.polypeptide_with_dihedrals = PolypeptideWithDihedrals()
        for res in self.polypeptide:
            self.polypeptide_with_dihedrals.append(res)
        phi_psi_list = self.polypeptide_with_dihedrals.get_phi_psi_list()
        self.assertIsInstance(phi_psi_list, dict)
        self.second_res_with_angles = next(islice(phi_psi_list.items(), 1, 2))  # first residue does not contain the phi/psi angles
        self.second_res = self.second_res_with_angles[0]
        self.second_res_angles = self.second_res_with_angles[1]

    def test_polypeptide(self):
        self.assertGreater(len(self.polypeptide_with_dihedrals), 0)  # add assertion here
        self.assertIsInstance(self.polypeptide_with_dihedrals, list)

    def test_get_phi_psi_list(self):
        #assert the angles list contains two elements: phi and psi
        self.assertEqual(len(self.second_res_angles), 2)
        #assert the phi angle is of float type
        self.assertIsInstance(self.second_res_angles[0], float)

    def test_assign_secondary_structure(self):
        self.polypeptide_with_dihedrals.assign_secondary_structure()
        ss_map = self.polypeptide_with_dihedrals.secondary_structure_map
        self.assertIsInstance(ss_map, dict)
        self.assertIn(ss_map[self.second_res], ('U', 'alpha', 'beta', 'b_par', 'b_anti'))


if __name__ == '__main__':
    unittest.main()
