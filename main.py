import argparse

from Bio.PDB import PDBParser

import utils
from advanced_polypeptide_builder import DihedralPPBuilder
from utils import plot_ramachandran

parser = PDBParser()


def parse_args():
    """Parse command-line arguments."""
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    arg_parser.add_argument(
        "-pdb", dest="pdb_file", required=True, help="Path to PDB file"
    )
    arg_parser.add_argument(
        "-res", dest="residue_number", required=False, help="Residue number"
    )
    arg_parser.add_argument(
        "-chain", dest="chain_name", required=False, default="A", help="Chain name"
    )

    arg_parser.add_argument(
        "-output",
        dest="ramachandran_output",
        required=False,
        default="ramachandran_plot.png",
        help="Output plot image name",
    )
    arg_parser.add_argument(
        "-plot", dest="plot", action="store_true", help="Plot Ramachandran plot"
    )
    arg_parser.add_argument("-export_csv", help="Path to export phi/psi angles to CSV")
    arg_parser.add_argument(
        "-ss_summary", action="store_true", help="Print secondary structure summary"
    )

    return arg_parser.parse_args()


def main():
    args = parse_args()
    pdb_file = args.pdb_file
    chain_name = args.chain_name
    structure = parser.get_structure(pdb_file, pdb_file)
    chain = [c for c in structure.get_chains() if c.id == chain_name][0]

    ppb = DihedralPPBuilder()

    peptides = ppb.build_peptides(chain, aa_only=1)

    phi_psi_data = {}
    for peptide in peptides:
        phi_psi_data.update(peptide.get_phi_psi_list())

    # create and save ramachandran plot
    if args.plot:
        plot_ramachandran(phi_psi_data.values(), args.ramachandran_output)

    # export to csv
    if args.export_csv:
        utils.export_phi_psi_to_csv(phi_psi_data, args.export_csv)

    # write secondary structure summary
    if args.ss_summary:
        ss_map = {}
        for peptide in peptides:
            ss_map.update(peptide.assign_secondary_structure())
        summary = utils.generate_ss_summary(ss_map)
        print(summary)


if __name__ == "__main__":
    main()
