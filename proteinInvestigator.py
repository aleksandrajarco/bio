#!/usr/bin/env python3
"""
Protein Secondary Structure Analyzer

Created on Jan 21, 2014
Author: ola
"""

import sys
import argparse
from math import pi
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Polypeptide import PPBuilder, Polypeptide
from Bio.PDB.vectors import calc_dihedral

parser = PDBParser()


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-pdb", dest="pdb_file", required=True, help="Path to PDB file")
    return parser.parse_args()


def unify_angle(angle):
    """Normalize angle to [-180, 180]"""
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360
    return angle


class Polypeptide2(Polypeptide):
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


class _PPBuilder2(PPBuilder):
    def build_peptides(self, entity, aa_only=1):
        """Build and return a list of Polypeptide2 objects."""
        level = entity.get_level()

        if level == "S":
            chains = entity[0].get_list()
        elif level == "M":
            chains = entity.get_list()
        elif level == "C":
            chains = [entity]
        else:
            raise PDBException("Entity must be Structure, Model, or Chain.")

        peptides = []

        for chain in chains:
            iterator = iter(chain)
            try:
                prev_res = next(iterator)
                while not self._accept(prev_res, aa_only):
                    prev_res = next(iterator)
            except StopIteration:
                continue

            pp = None
            for next_res in iterator:
                if (self._accept(prev_res, aa_only) and self._accept(next_res, aa_only)
                        and self._is_connected(prev_res, next_res)):
                    if pp is None:
                        pp = Polypeptide2()
                        pp.append(prev_res)
                        peptides.append(pp)
                    pp.append(next_res)
                else:
                    pp = None
                prev_res = next_res

        return peptides


class ProteinInvestigation:
    def __init__(self, res_nr, chain_name, pdb_file):
        self.resNr = res_nr
        self.chainName = chain_name
        self.pdbFile = pdb_file
        self.structure = parser.get_structure(pdb_file, pdb_file)
        self.chain = None
        self.residue = None

    def if_chain_in_structure(self):
        self.chain = [c for c in self.structure[0] if c.id == self.chainName]
        return bool(self.chain)

    def if_structure_in_file(self):
        return bool([r for r in self.structure.get_residues()])

    def find_res_in_structure(self):
        if not self.if_structure_in_file() or not self.if_chain_in_structure():
            print('Cannot find residue in structure.')
            sys.exit()

        try:
            self.residue = [
                r for r in self.structure.get_residues()
                if r.id[1] == int(self.resNr)
            ][0]
        except IndexError:
            print('Residue not found.')
            sys.exit()

        return self.residue

    def calc_dihedral(self, atom1, atom2, atom3, atom4):
        self.find_res_in_structure()
        angle = calc_dihedral(
            atom1.get_vector(),
            atom2.get_vector(),
            atom3.get_vector(),
            atom4.get_vector()
        ) * 180 / pi
        return angle

    def define_ss(self, polypeptide):
        """Define secondary structure based on phi/psi angles."""
        phi_psi_dict = polypeptide.get_phi_psi_list()

        for residue, (phi, psi) in phi_psi_dict.items():
            ss = 'U'
            if phi is not None and psi is not None:
                if -145 < phi < -135 and 130 < psi < 140:
                    ss = 'b_anti'
                elif -125 < phi < -115 and 110 < psi < 120:
                    ss = 'b_par'
                elif -180 < phi < -100 and 100 < psi < 180:
                    ss = 'beta'
                elif -135 < phi < -45 and -50 < psi < 40:
                    ss = 'alpha'
            print(residue.resname, residue.id[1], ss, phi, psi)


def main():
    # args = parse_args()
    pdb_file = 'examples/9lil.pdb'  # Hardcoded for now
    chain_name = 'A'
    res_nr = 238

    structure = parser.get_structure(pdb_file, pdb_file)
    chain_a = [c for c in structure.get_chains() if c.id == chain_name][0]

    ppb = _PPBuilder2()
    peptide = ppb.build_peptides(chain_a, aa_only=1)[0]

    PI = ProteinInvestigation(res_nr, chain_name, pdb_file)
    PI.define_ss(peptide)


if __name__ == '__main__':
    main()
