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
- [matplotlib](https://matplotlib.org/) (optional, for plotting)

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
   pip install requirements.txt
   ```
   
### Usage
To analyze a protein structure:

```bash
python main.py -pdb examples/9lil.pdb -chain A
```
This will:

Parse the specified PDB file.

Compute the ϕ and ψ angles for each residue.

Assign secondary structure elements based on these angles.

Output the results to the console.

## Features
Dihedral Angle Calculation: Computes backbone ϕ and ψ angles for residues.

Secondary Structure Assignment: Classifies residues into secondary structures using angle thresholds.

Visualization: (Planned) Generate plots to visualize secondary structure distribution and Ramachandran plots.

## Future Enhancements
Implement visualization tools for secondary structure distribution.

Develop Ramachandran plot generation.

Expand secondary structure classification criteria.

Integrate with molecular visualization tools like PyMOL or NGLView.


