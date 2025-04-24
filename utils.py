from collections import Counter

from numpy import pi
from Bio.PDB import vectors
import matplotlib.pyplot as plt


def plot_ramachandran(phi_psi_angles, output_file="ramachandran_plot.png", title=None):
    """
    Saves a Ramachandran plot of phi/psi angles.

    Parameters:
    - phi_psi_angles: List of tuples (phi, psi). Can contain None values.
    - output_path: Path to save the PNG image.
    - title: Optional chart title.
    """
    # Filter out None values
    filtered_angles = [
        (phi, psi) for phi, psi in phi_psi_angles if phi is not None and psi is not None
    ]

    if not filtered_angles:
        raise ValueError("No valid phi/psi angles to plot.")

    phi_vals, psi_vals = zip(*filtered_angles)

    plt.figure(figsize=(6, 6))
    plt.scatter(phi_vals, psi_vals, alpha=0.7, c="blue", edgecolors="k")
    plt.xlabel("Phi (ϕ) angles")
    plt.ylabel("Psi (ψ) angles")
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)
    plt.axhline(0, color="gray", linewidth=0.5)
    plt.axvline(0, color="gray", linewidth=0.5)

    if title:
        plt.title(title)

    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def unify_angle(angle):
    """Normalize angle to [-180, 180]"""
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360
    return angle


def calculate_dihedral_angle(atom1, atom2, atom3, atom4):
    angle = (
        vectors.calc_dihedral(
            atom1.get_vector(),
            atom2.get_vector(),
            atom3.get_vector(),
            atom4.get_vector(),
        )
        * 180
        / pi
    )
    return angle


def export_phi_psi_to_csv(phi_psi_data, output_path):
    with open(output_path, "w") as f:
        f.write("residue_id,residue_name,phi,psi\n")
        for res, angles in phi_psi_data.items():
            res.id = res.id[1]
            res_name = res.resname
            phi, psi = angles
            f.write(f"{res.id},{res_name},{phi},{psi}\n")
        # for residue_id, res_name, phi, psi in phi_psi_data:
        # f.write(f"{residue_id},{res_name},{phi},{psi}\n")


def generate_ss_summary(ss_map):
    counts = Counter(ss_map.values())
    summary = "Secondary Structure Summary:\n"
    for ss_type in ("alpha", "beta", "b_par", "b_anti", "U"):
        summary += f"{ss_type:}: {counts.get(ss_type, 0)}\n"
    return summary
