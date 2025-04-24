import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch
from numpy import pi

from polypeptide_with_dihedrals import PolypeptideWithDihedrals
from utils import (
    unify_angle,
    calculate_dihedral_angle,
    plot_ramachandran,
    export_phi_psi_to_csv,
    generate_ss_summary,
)


class TestUtils(unittest.TestCase):
    def setUp(self):
        # Mock atoms with mock vectors
        self.atom1 = MagicMock()
        self.atom2 = MagicMock()
        self.atom3 = MagicMock()
        self.atom4 = MagicMock()

        self.mock_vector = MagicMock()
        for atom in [self.atom1, self.atom2, self.atom3, self.atom4]:
            atom.get_vector.return_value = self.mock_vector

        self.residue = MagicMock()
        self.residue.get_resname.return_value = "ALA"
        self.residue.get_id.return_value = (" ", 1, " ")
        # dummy phi_psi_map
        self.phi_psi_map = {self.residue: [0.123456, 0.987654]}
        # dummy ss_map
        self.ss_map = {self.residue: "U"}

    def test_unify_angle(self):
        self.assertEqual(unify_angle(230), -130)
        self.assertEqual(unify_angle(-170), -170)
        self.assertEqual(unify_angle(-200), 160)

    @patch("utils.vectors.calc_dihedral", return_value=pi / 2)
    def test_calculate_dihedral_angle(self, mock_calc):
        result = calculate_dihedral_angle(
            self.atom1, self.atom2, self.atom3, self.atom4
        )
        self.assertAlmostEqual(result, 90.0)
        mock_calc.assert_called_once()  # optional, to verify interaction

    def test_ramachandran_plot_saves_file(self):
        # Dummy phi/psi angle data (should be valid floats)
        dummy_angles = [
            (-60.0, -45.0),
            (-120.0, 135.0),
            (-135.0, 140.0),
            (-75.0, 160.0),
        ]

        # Create a temporary file path
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_file:
            output_path = tmp_file.name

        try:
            # Call the function
            plot_ramachandran(dummy_angles, output_path)

            # Assert file was created
            self.assertTrue(os.path.exists(output_path))
            self.assertGreater(os.path.getsize(output_path), 0)  # Not an empty file

        finally:
            # Cleanup
            if os.path.exists(output_path):
                os.remove(output_path)

    def test_export_phi_psi_to_csv(self):

        # Create a temporary file path
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp_file:
            output_path = tmp_file.name

        try:
            export_phi_psi_to_csv(self.phi_psi_map, output_path)
            self.assertTrue(os.path.exists(output_path))
            self.assertGreater(os.path.getsize(output_path), 0)  # Not an empty file

            # Check that the file contains expected data (example: first residue's angles)
            with open(output_path, "r") as file:
                file_content = file.read()

            # Assert that the content includes phi and psi angles for residue1
            self.assertIn("residue_id", file_content)
            self.assertIn("residue_name", file_content)
            self.assertIn("0.123456", file_content)
            self.assertIn("0.987654", file_content)

        finally:
            # Cleanup
            if os.path.exists(output_path):
                os.remove(output_path)

    def test_generate_ss_summary(self):
        self.polypeptide = PolypeptideWithDihedrals()
        self.polypeptide.secondary_structure_map = {
            "residue1": "alpha",
            "residue2": "beta",
            "residue3": "b_par",
            "residue4": "b_anti",
            "residue5": "U",
        }
        # Call the function
        summary = generate_ss_summary(self.ss_map)

        # Check the output
        self.assertIn("alpha", summary)
        self.assertIn("beta", summary)
        self.assertIn("b_par", summary)
        self.assertIn("b_anti", summary)
        self.assertIn("U", summary)


if __name__ == "__main__":
    unittest.main()
