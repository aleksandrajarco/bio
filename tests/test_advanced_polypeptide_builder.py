import unittest
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.Entity import Entity
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Structure import Structure

from advanced_polypeptide_builder import (
    DihedralPPBuilder,
)  # Assuming your builder is imported like this
from polypeptide_with_dihedrals import PolypeptideWithDihedrals


class TestDihedralPPBuilder(unittest.TestCase):
    def setUp(self):

        base_dir = (
            Path(__file__).resolve().parent.parent
        )  # go up from tests/ to project root
        pdb_path = base_dir / "examples" / "9lil.pdb"

        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("example", pdb_path)

    def test_build_peptides(self):
        # Initialize PPBuilder
        ppb = DihedralPPBuilder()

        # Build polypeptides
        polypeptides = ppb.build_peptides(self.structure)

        # Assertions
        self.assertIsInstance(polypeptides, list)
        self.assertEqual(len(polypeptides), 2)

        self.assertIsInstance(polypeptides[0], PolypeptideWithDihedrals)

    def test_first_peptide_residues(self):
        peptides = DihedralPPBuilder().build_peptides(self.structure)
        self.assertGreater(len(peptides[0]), 0)  # Make sure it contains residues

    def test_build_peptides_empty_structure(self):
        empty_structure = Structure("empty_structure")
        builder = DihedralPPBuilder()
        with self.assertRaises(PDBException) as context:
            builder.build_peptides(empty_structure)

        self.assertIn("Structure contains no models", str(context.exception))


if __name__ == "__main__":
    unittest.main()
