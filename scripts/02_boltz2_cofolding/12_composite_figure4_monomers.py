#!/usr/bin/env python3
"""
Step 4b — Figure 4: composite the three ChimeraX panels produced by
11_render_figure4_monomers.cxc into the final two-row figure
(A/B stacked monomers, C superposition spanning the bottom row).

Run 11_render_figure4_monomers.cxc in ChimeraX first.

Run (from the repository root):
    python scripts/02_boltz2_cofolding/12_composite_figure4_monomers.py

Output:
    generated/figures/Figure4_composite.png / .tiff
"""
import os
from PIL import Image, ImageChops
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

FIGURES_DIR = "generated/figures"


def content_bbox(img, padding=20):
    bg = Image.new(img.mode, img.size, (255, 255, 255))
    diff = ImageChops.difference(img.convert("RGB"), bg.convert("RGB"))
    bbox = diff.getbbox()
    if bbox is None:
        return (0, 0, img.width, img.height)
    return (max(0, bbox[0] - padding), max(0, bbox[1] - padding),
            min(img.width, bbox[2] + padding), min(img.height, bbox[3] + padding))


def autocrop(img, padding=20):
    return img.crop(content_bbox(img, padding))


def crop_common_scale(imgs, padding=20):
    """Crop each image to the same window size, centered on its own content,
    so relative molecular scale between panels is preserved. An independent
    per-panel autocrop makes ChimeraX's per-model auto-fit view invisible —
    e.g. CsTx3's C-terminal extension gives it a looser bounding box than
    CsTx2's, which would otherwise render CsTx3 smaller once both panels are
    placed in an equal-sized grid cell."""
    bboxes = [content_bbox(im, padding) for im in imgs]
    w = max(b[2] - b[0] for b in bboxes)
    h = max(b[3] - b[1] for b in bboxes)
    out = []
    for im, (x0, y0, x1, y1) in zip(imgs, bboxes):
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        left = max(0, min(im.width - w, int(cx - w / 2)))
        top = max(0, min(im.height - h, int(cy - h / 2)))
        out.append(im.crop((left, top, left + w, top + h)))
    return out


def main():
    img_a, img_b = crop_common_scale([
        Image.open(os.path.join(FIGURES_DIR, "Figure4_panelA_CsTx2.png")),
        Image.open(os.path.join(FIGURES_DIR, "Figure4_panelB_CsTx3.png")),
    ])
    img_c = autocrop(Image.open(os.path.join(FIGURES_DIR, "Figure4_panelC_superposition.png")))

    fig = plt.figure(figsize=(12, 10), dpi=600)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.04, wspace=0.04)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    for ax, img, label in [(ax_a, img_a, "A"), (ax_b, img_b, "B"), (ax_c, img_c, "C")]:
        ax.imshow(img)
        ax.axis("off")
        ax.text(0.02, 0.98, label, transform=ax.transAxes,
                 fontsize=18, fontweight="bold", va="top", ha="left", color="black")

    os.makedirs(FIGURES_DIR, exist_ok=True)
    out_png = os.path.join(FIGURES_DIR, "Figure4_composite.png")
    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"wrote {out_png}")

    out_tiff = os.path.join(FIGURES_DIR, "Figure4_composite.tiff")
    Image.open(out_png).save(out_tiff, compression="tiff_lzw", dpi=(600, 600))
    print(f"wrote {out_tiff}")


if __name__ == "__main__":
    main()
