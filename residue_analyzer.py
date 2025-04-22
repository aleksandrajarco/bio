#!/usr/bin/env python3
"""
Protein Secondary Structure Analyzer

Created on Jan 21, 2014
Author: ola
"""

import sys

class ChainNotFoundError(Exception):
    """Custom exception for when a chain is not found in the structure."""
    def __init__(self, message="Cannot find residue in structure."):
        self.message = message
        super().__init__(self.message)

class ResidueNotFoundError(Exception):
    """Custom exception for when a residue is not found in the structure."""
    def __init__(self, message="Residue not found."):
        self.message = message
        super().__init__(self.message)

class ResidueAnalyzer:
    def __init__(self, residue_number, chain_name, structure):
        self.residue_number = residue_number
        self.chain_name = chain_name
        self.structure = structure
        self.chain = chain_name
        self.residue = None

    def if_chain_in_structure(self):
        self.chain = [c for c in self.structure[0] if c.id == self.chain_name]
        return bool(self.chain)

    def find_res_in_structure(self):
        # Check if the chain exists in the structure
        if not self.if_chain_in_structure():
            raise ChainNotFoundError("Cannot find residue in structure.")

        try:
            # Try to find the residue in the structure
            self.residue = [
                r for r in self.structure.get_residues()
                if r.id[1] == int(self.residue_number)
            ][0]
        except IndexError:
            raise ResidueNotFoundError("Residue not found.")

        return self.residue

