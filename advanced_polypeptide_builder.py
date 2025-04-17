from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Polypeptide import PPBuilder

from polypeptide_with_dihedrals import PolypeptideWithDihedrals

class DihedralPPBuilder(PPBuilder):
    def build_peptides(self, entity, aa_only=1):
        """Build and return a list of PolypeptideWithDihedrals objects."""
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
                        pp = PolypeptideWithDihedrals()
                        pp.append(prev_res)
                        peptides.append(pp)
                    pp.append(next_res)
                else:
                    pp = None
                prev_res = next_res

        return peptides