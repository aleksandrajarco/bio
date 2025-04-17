import argparse

from Bio.PDB import PDBParser

from advanced_polypeptide_builder import DihedralPPBuilder
from residue_analyzer import ResidueAnalyzer

parser = PDBParser()

def parse_args():
    """Parse command-line arguments."""
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg_parser.add_argument("-pdb", dest="pdb_file", required=True, help="Path to PDB file")
    arg_parser.add_argument("-res", dest="residue_number", required=False, help="Residue number")
    arg_parser.add_argument("-chain", dest="chain_name", required=False,  default = 'A',  help="Chain name")
    return arg_parser.parse_args()

def main():
    args = parse_args()
    pdb_file = args.pdb_file
    chain_name = args.chain_name
    structure = parser.get_structure(pdb_file, pdb_file)
    chain_a = [c for c in structure.get_chains() if c.id == chain_name][0]

    ppb = DihedralPPBuilder()
    peptide = ppb.build_peptides(chain_a, aa_only=1)[0]
    ss_map = peptide.assign_secondary_structure()
    residue_analyzer = ResidueAnalyzer(238, 'A', structure)
    residue = residue_analyzer.find_res_in_structure()
    print(residue)

    #print(ss_map)

if __name__ == '__main__':
    main()
