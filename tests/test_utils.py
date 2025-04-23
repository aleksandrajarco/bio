import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch
from numpy import pi

from utils import unify_angle, calculate_dihedral_angle, plot_ramachandran


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
            (-75.0, 160.0)
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

if __name__ == "__main__":
    unittest.main()
