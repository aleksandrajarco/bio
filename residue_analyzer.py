#!/usr/bin/env python3
"""
Protein Secondary Structure Analyzer

Created on Jan 21, 2014
Author: ola
"""

import sys
import argparse
from Bio.PDB import PDBParser
from advanced_polypeptide_builder import DihedralPPBuilder

class ResidueAnalyzer:
    def __init__(self, residue_number, chain_name, pdb_file):
        self.residue_number = residue_number
        self.chain_name = chain_name
        self.pdb_file = pdb_file
        #self.structure = parser.get_structure(pdb_file, pdb_file)
        self.chain = None
        self.residue = None

    def if_chain_in_structure(self):
        self.chain = [c for c in self.pdb_file[0] if c.id == self.chain_name]
        return bool(self.chain)

    def find_res_in_structure(self):
        if not self.if_chain_in_structure():
            print('Cannot find residue in structure.')
            sys.exit()

        try:
            self.residue = [
                r for r in self.pdb_file.get_residues()
                if r.id[1] == int(self.residue_number)
            ][0]
        except IndexError:
            print('Residue not found.')
            sys.exit()

        return self.residue
