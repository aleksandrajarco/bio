# Protein Secondary Structure Analyzer

This project provides tools for analyzing protein secondary structures by computing backbone dihedral angles (ϕ and ψ) from PDB files. Leveraging Biopython, it assigns secondary structure elements (e.g., α-helix, β-sheet) based on these angles.

## Repository Overview

Key features:
- Parses PDB files and extracts chains and residues
- Computes phi/psi angles using backbone atoms
- Assigns secondary structure classes (alpha, beta, b_anti, b_par, unstructured)
- Built on top of Bio.PDB with modular, extensible architecture

## Getting Started

### Prerequisites

- Python 3.8 or higher
- [Biopython](https://biopython.org/)
- matplotlib >= 3.4.0

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/aleksandrajarco/bio.git
   cd bio
   ```

2. Create a virtual environment (optional but recommended):

   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```	
3. Install the required packages:

   ```bash
   pip install -r requirements.txt
   ```
   
### Usage
## Example CLI Commands:
1. Generate Ramachandran Plot for a Chain
This command generates a Ramachandran plot for a given PDB file and chain.

```bash
python main.py -pdb path/to/your/pdb_file.pdb -chain A -output plot.png
```
Arguments:

-pdb: Path to the PDB file.

-chain: The chain to analyze.

-output: Optional. File name to save the plot to (default is ramachandran_plot.png).

2. Export Phi/Psi Angles to CSV
This command exports the Phi/Psi angles of the peptide chain to a CSV file.

```bash
python main.py -pdb path/to/your/pdb_file.pdb -chain A -export_csv
```
Arguments:

-pdb: Path to the PDB file.

-chain: The chain to analyze.

-export_csv: path to CSV file where the phi/psi angles together with residue information will be exported..

3. Display a Summary of Secondary Structures
This command outputs a summary of secondary structures in the given PDB file and chain.

```bash
python main.py -pdb path/to/your/pdb_file.pdb -chain A -ss_summary
```
Arguments:

-pdb: Path to the PDB file.

-chain: The chain to analyze.

-ss_summary: Flag to display secondary structure summary.

## Features
Dihedral Angle Calculation: Computes backbone ϕ and ψ angles for residues.

Secondary Structure Assignment: Classifies residues into secondary structures using angle thresholds.

Visualization: Generate plots to visualize secondary structure distribution and Ramachandran plots.

## Future Enhancements
Integration with Data Analysis Libraries: Consider integrating the tool with libraries like pandas or scipy for data analysis or better handling of larger datasets.

More Plotting Options: Besides Ramachandran plots, providing options for other plots such as secondary structure distributions might be a good enhancement.

Handling Multiple Chains: Add a feature for handling multiple chains in the PDB file for batch processing.

3D Plotting of Secondary Structures: This could be a more advanced visualization feature where secondary structure elements (alpha-helices, beta-sheets) are visualized in 3D.

