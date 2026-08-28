#!/usr/bin/env python3
"""
Supplementary Figure S3 — Ramachandran analysis of the Boltz-2-predicted
structures of CsTx2 and CsTx3 only (Cs32.1 and CsTx1 have no structural
model: Cs32.1 was not carried forward past purification, and CsTx1's
C-terminus was never resolved by Edman/MALDI sequencing).

Uses the best-confidence Boltz-2 sample per peptide (model_0.pdb, selected by
scripts/02_boltz2_cofolding/03_evaluate_boltz.py) and a PROCHECK-style four-zone
background (Morris et al. 1992; Lovell et al. 2003).

Run (from the repository root):
    python scripts/05_figures_and_statistics/S3_render_figure_ramachandran.py

Output:
    generated/figures/FigureS3_ramachandran.png / .pdf / .tiff
"""
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MplPath
from Bio.PDB import PDBParser, PPBuilder
from PIL import Image

FIGURES_DIR = "generated/figures"
PREDICTIONS_DIR = "generated/predictions"

PEPTIDES = {
    "CsTx2": (os.path.join(PREDICTIONS_DIR, "CsTx2", "CsTx2_model_0.pdb"), "#2a78d6"),
    "CsTx3": (os.path.join(PREDICTIONS_DIR, "CsTx3", "CsTx3_model_0.pdb"), "#eb6834"),
}


def _polygon(*points):
    closed = list(points) + [points[0]]
    return MplPath(closed)


MOST_FAVORED = [
    _polygon((-155, -65), (-155, 20), (-100, 52), (-60, 52), (-40, 20), (-40, -65)),
    _polygon((-160, 110), (-160, 180), (-50, 180), (-50, 110)),
    _polygon((-160, -180), (-160, -140), (-100, -140), (-100, -180)),
    _polygon((30, -40), (30, 65), (90, 65), (90, -40)),
]
ADDITIONALLY_ALLOWED = [
    _polygon((-175, -85), (-175, 38), (-25, 38), (-25, -85)),
    _polygon((-175, 80), (-175, 180), (-25, 180), (-25, 80)),
    _polygon((-175, -180), (-175, -110), (-40, -110), (-40, -180)),
    _polygon((-90, 50), (-90, 120), (-35, 120), (-35, 50)),
    _polygon((-175, 45), (-175, 115), (-90, 115), (-90, 45)),
    _polygon((25, -55), (25, 80), (105, 80), (105, -55)),
]
GENEROUSLY_ALLOWED = [
    _polygon((-180, -100), (-180, 50), (-10, 50), (-10, -100)),
    _polygon((-180, 60), (-180, 180), (0, 180), (0, 60)),
    _polygon((-180, -180), (-180, -95), (0, -95), (0, -180)),
    _polygon((15, -70), (15, 100), (120, 100), (120, -70)),
    _polygon((-115, 25), (-115, 150), (-10, 150), (-10, 25)),
]


def draw_ramachandran_regions(ax):
    for poly in GENEROUSLY_ALLOWED:
        v = poly.vertices
        ax.fill(v[:, 0], v[:, 1], color="#b749e6", alpha=0.35, zorder=0)
    for poly in ADDITIONALLY_ALLOWED:
        v = poly.vertices
        ax.fill(v[:, 0], v[:, 1], color="#00e5e5", alpha=0.55, zorder=1)
    for poly in MOST_FAVORED:
        v = poly.vertices
        ax.fill(v[:, 0], v[:, 1], color="#1df892", alpha=0.70, zorder=2)


def procheck_zone(phi, psi):
    point = np.array([[phi, psi]])
    if any(p.contains_points(point)[0] for p in MOST_FAVORED):
        return "most_favored"
    if any(p.contains_points(point)[0] for p in ADDITIONALLY_ALLOWED):
        return "additional_allowed"
    if any(p.contains_points(point)[0] for p in GENEROUSLY_ALLOWED):
        return "generously_allowed"
    return "disallowed"


def compute_phi_psi(structure):
    builder = PPBuilder()
    result = []
    for peptide in builder.build_peptides(structure):
        for residue, (phi, psi) in zip(peptide, peptide.get_phi_psi_list()):
            if phi is not None and psi is not None:
                result.append((residue.get_resname(), math.degrees(phi), math.degrees(psi)))
    return result


def main():
    parser = PDBParser(QUIET=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))

    plt.rcParams.update({
        "font.size": 10, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
        "text.color": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    })

    stats = {}
    for ax, label, (name, (pdb_file, color)) in zip(axes, ("A", "B"), PEPTIDES.items()):
        structure = parser.get_structure(name, pdb_file)
        phi_psi = compute_phi_psi(structure)
        draw_ramachandran_regions(ax)
        phis = [p for _, p, _ in phi_psi]
        psis = [p for _, _, p in phi_psi]
        zones = [procheck_zone(p, s) for _, p, s in phi_psi]
        n_favored = sum(1 for z in zones if z == "most_favored")
        n_allowed = sum(1 for z in zones if z in ("most_favored", "additional_allowed"))
        stats[name] = (n_favored, n_allowed, len(phi_psi))
        ax.scatter(phis, psis, c=color, s=70, edgecolors="black", linewidths=0.6, zorder=3, alpha=0.9)
        for resname, phi, psi in phi_psi:
            ax.annotate(resname[:1], (phi, psi), fontsize=6, alpha=0.75, xytext=(3, 3),
                        textcoords="offset points")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-180, 181, 60))
        ax.axhline(0, color="gray", lw=0.5, zorder=1)
        ax.axvline(0, color="gray", lw=0.5, zorder=1)
        ax.set_xlabel(r"$\phi$ (degrees)")
        ax.set_ylabel(r"$\psi$ (degrees)")
        ax.set_aspect("equal")
        # Panel letter only — the caption identifies which panel is which
        # peptide, matching the convention used in every composited figure.
        ax.text(0.02, 0.98, label, transform=ax.transAxes, fontsize=18, fontweight="bold",
                va="top", ha="left", color="black")

    legend_elements = [
        mpatches.Patch(facecolor="#1df892", alpha=0.85, label="Most favored"),
        mpatches.Patch(facecolor="#00e5e5", alpha=0.85, label="Additional allowed"),
        mpatches.Patch(facecolor="#b749e6", alpha=0.85, label="Generously allowed"),
        mpatches.Patch(facecolor="white", edgecolor="gray", label="Disallowed"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.02), fontsize=9)
    fig.tight_layout()

    os.makedirs(FIGURES_DIR, exist_ok=True)
    out_png = os.path.join(FIGURES_DIR, "FigureS3_ramachandran.png")
    out_pdf = os.path.join(FIGURES_DIR, "FigureS3_ramachandran.pdf")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    out_tiff = os.path.join(FIGURES_DIR, "FigureS3_ramachandran.tiff")
    Image.open(out_png).save(out_tiff, compression="tiff_lzw", dpi=(600, 600))
    print(f"wrote {out_png}, {out_pdf}, {out_tiff}")
    for name, (n_favored, n_allowed, n_total) in stats.items():
        print(f"  {name}: {n_favored}/{n_total} most favored, {n_allowed}/{n_total} allowed")


if __name__ == "__main__":
    main()
