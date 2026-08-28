#!/usr/bin/env python3
"""
Step 5b — Figure 5: composite the co-folding pose render (panel A) with the
pose-convergence-vs-ipTM scatter plot (panel B, from the Boltz-2 co-folding
statistics) into the final two-panel figure.

Run 13_render_figure5_pose.cxc in ChimeraX first, and run
scripts/05_figures_and_statistics/07_make_cofolding_convergence_figure.py to
produce generated/figures/fig_cofolding_convergence_vs_iptm.png (see the
top-level README).

Run (from the repository root):
    python scripts/02_boltz2_cofolding/14_composite_figure5_pose.py

Output:
    generated/figures/Figure5_composite.png / .tiff
"""
import os
from PIL import Image, ImageChops
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FIGURES_DIR = "generated/figures"


def autocrop(img, padding=20):
    bg = Image.new(img.mode, img.size, (255, 255, 255))
    diff = ImageChops.difference(img.convert("RGB"), bg.convert("RGB"))
    bbox = diff.getbbox()
    if bbox:
        bbox = (max(0, bbox[0] - padding), max(0, bbox[1] - padding),
                min(img.width, bbox[2] + padding), min(img.height, bbox[3] + padding))
        return img.crop(bbox)
    return img


def main():
    img_a = autocrop(Image.open(os.path.join(FIGURES_DIR, "Figure5_panelA_pose.png")))
    img_b = autocrop(Image.open(os.path.join(FIGURES_DIR, "fig_cofolding_convergence_vs_iptm.png")))

    fig = plt.figure(figsize=(14, 7), dpi=600)
    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1, 1.1], wspace=0.04)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    for ax, img, label in [(ax_a, img_a, "A"), (ax_b, img_b, "B")]:
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.02, 0.98, label, transform=ax.transAxes,
                 fontsize=18, fontweight="bold", va="top", ha="left", color="black")

    os.makedirs(FIGURES_DIR, exist_ok=True)
    out_png = os.path.join(FIGURES_DIR, "Figure5_composite.png")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"wrote {out_png}")

    out_tiff = os.path.join(FIGURES_DIR, "Figure5_composite.tiff")
    Image.open(out_png).save(out_tiff, compression="tiff_lzw", dpi=(600, 600))
    print(f"wrote {out_tiff}")


if __name__ == "__main__":
    main()
