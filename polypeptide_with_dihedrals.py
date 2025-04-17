from Bio.PDB import calc_dihedral
from Bio.PDB.Polypeptide import Polypeptide
from utils import unify_angle


class PolypeptideWithDihedrals(Polypeptide):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.secondary_structure_map = {}

    def get_phi_psi_list(self):
        """Return a dictionary of phi/psi angles per residue."""
        angles = {}
        for i, res in enumerate(self):
            try:
                n = res['N'].get_vector()
                ca = res['CA'].get_vector()
                c = res['C'].get_vector()
            except KeyError:
                angles[res] = (None, None)
                res.xtra["PHI"] = res.xtra["PSI"] = None
                continue

            phi = psi = None

            if i > 0:
                try:
                    cp = self[i - 1]['C'].get_vector()
                    phi = calc_dihedral(cp, n, ca, c)
                except KeyError:
                    pass

            if i < len(self) - 1:
                try:
                    nn = self[i + 1]['N'].get_vector()
                    psi = calc_dihedral(n, ca, c, nn)
                except KeyError:
                    pass

            phi_deg = unify_angle(phi * 180) if phi is not None else None
            psi_deg = unify_angle(psi * 180) if psi is not None else None
            angles[res] = (phi_deg, psi_deg)
            res.xtra["PHI"], res.xtra["PSI"] = phi_deg, psi_deg

        return angles

    def assign_secondary_structure(self):
        """Define secondary structure based on phi/psi angles."""
        phi_psi_dict = self.get_phi_psi_list()

        for residue, (phi, psi) in phi_psi_dict.items():
            secondary_structure = 'U'
            if phi is not None and psi is not None:
                if -145 < phi < -135 and 130 < psi < 140:
                    secondary_structure = 'b_anti'
                elif -125 < phi < -115 and 110 < psi < 120:
                    secondary_structure = 'b_par'
                elif -180 < phi < -100 and 100 < psi < 180:
                    secondary_structure = 'beta'
                elif -135 < phi < -45 and -50 < psi < 40:
                    secondary_structure = 'alpha'

            self.secondary_structure_map[residue] = secondary_structure

        return self.secondary_structure_map
