import unittest

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from residue_analyzer import ResidueAnalyzer, ResidueNotFoundError


class TestResidueAnalyzer(unittest.TestCase):
    def setUp(self):

        # Create a new Structure
        structure = Structure("example")

        # Create a Model and add it to the Structure
        model = Model(0)
        structure.add(model)

        # Create a Chain and add it to the Model
        chain = Chain('A')
        model.add(chain)

        # Create a Residue and add it to the Chain
        # We'll use a simple "ALA" (Alanine) residue with a few atoms
        self.residue = Residue((' ', 1, ' '), 'ALA', ' ')
        chain.add(self.residue)

        # Create Atoms for the residue
        # For simplicity, we are adding only the CA (alpha-carbon) atom here
        atom_ca = Atom('CA', np.array([0.0, 0.0, 0.0]), 1.0, 1.0, ' ', 'CA', 1, 'CA', 1.0)

        self.residue.add(atom_ca)

        # You can add more residues or atoms similarly
        # Add another residue, for example "GLY" (Glycine)
        self.residue2 = Residue((' ', 2, ' '), 'GLY', ' ')
        atom_n = Atom('N', np.array([1.0, 0.0, 0.0]), 1.0, 1.0, ' ', 'N', 1, 'N', 1.0)
        self.residue2.add(atom_n)

        self.residue_analyzer = ResidueAnalyzer(residue_number=1, chain_name = 'A',  structure=structure)
        self.residue_analyzer2 = ResidueAnalyzer(residue_number=2, chain_name = 'A',  structure=structure)

    def test_found_res_in_structure(self):
        self.assertEqual(self.residue_analyzer.find_res_in_structure(), self.residue)

    def test_not_found_res_in_structure(self):
        # Test if the custom error is raised when residue is not found
        with self.assertRaises(ResidueNotFoundError):
            self.residue_analyzer2.find_res_in_structure()

if __name__ == '__main__':
    unittest.main()
