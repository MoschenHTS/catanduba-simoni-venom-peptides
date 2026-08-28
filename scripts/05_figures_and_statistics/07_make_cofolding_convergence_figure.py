#!/usr/bin/env python3
"""
Boltz-2 co-folding finding: ipTM does not discriminate real complexes from
scrambled-toxin negatives, but ensemble pose-convergence does. Scatter: x=best
ipTM, y=pose convergence %, two groups (scramble null vs real complex), same
categorical color convention as the MD figures (04_make_figures.py) for the 4
real complexes; neutral gray for the null.

DATA PROVENANCE. The values below are literals, unlike the membrane-MD figures
which read results/*.csv (01_collect_data.py). That difference is deliberate: the
MD numbers changed when replicas were extended n=3 -> n=8, so they needed a
drift-proof single source; the Boltz-2 co-folding campaign is complete and its
inputs are frozen. They are not unverified, though - every value here was
regenerated from the raw prediction ensembles and matched exactly against
results/boltz2_cofolding_summary.txt (committed). To re-verify at any time:

    python scripts/02_boltz2_cofolding/08_aggregate_cofolding_statistics.py   # from the repository root

and compare against results/boltz2_cofolding_summary.txt's "MULTI-SEED
REPRODUCIBILITY AND SCRAMBLE-NULL STATISTICS" section. If you ever re-run
Boltz-2 with different seeds or samples, update these literals from that
output rather than editing them by hand.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# Where the figure is written. Defaults to generated/figures/ at the repository
# root; override with CSIMONI_FIGDIR to write elsewhere.
FIGDIR = os.environ.get("CSIMONI_FIGDIR", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "generated", "figures"))
COLOR = {"CsTx2–Kv4.2": "#2a78d6", "CsTx2–Kv4.3": "#1baf7a",
         "CsTx3–Kv4.2": "#eda100", "CsTx3–Kv4.3": "#008300"}
NULL_COLOR = "#8a8a86"

# scramble null: seed42 (original) + 5 independent shuffles (shuf1-5)
null_iptm = [0.92, 0.93, 0.92, 0.93, 0.89, 0.88]
null_conv = [22, 8, 54, 10, 6, 6]

# real complexes: multi-seed means (seeds 15/13/1/999), best ipTM mean +- sd,
# convergence mean +- sd
real = {
    "CsTx2–Kv4.2": dict(iptm=(0.87, 0.01), conv=(72, 4)),
    "CsTx2–Kv4.3": dict(iptm=(0.89, 0.01), conv=(58, 6)),
    "CsTx3–Kv4.2": dict(iptm=(0.89, 0.01), conv=(84, 3)),
    "CsTx3–Kv4.3": dict(iptm=(0.91, 0.01), conv=(72, 5)),
}
# Reference line = the null MAXIMUM (54%), not mean+2SD. mean+2SD assumes an
# approximately Gaussian null; ddof also flips it between 52% (population) and 55.2%
# (sample) for this n=6 - an artifact that shouldn't decide a verdict. The null max is
# ddof-independent and is what "all four real complexes exceed the entire scramble null"
# actually means.
null_max = np.max(null_conv)

mpl.rcParams.update({
    "font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e5e4df", "grid.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

fig, ax = plt.subplots(figsize=(7.2, 5.8))

ax.axhspan(0, null_max, color=NULL_COLOR, alpha=0.08, zorder=0)
ax.axhline(null_max, color=NULL_COLOR, lw=1.2, ls="--", zorder=1)
ax.text(0.933, null_max + 2, f"scramble-null maximum\n({null_max:.0f}%, n=6)",
        fontsize=8, color=NULL_COLOR, style="italic", ha="right", va="bottom")

ax.scatter(null_iptm, null_conv, s=90, color=NULL_COLOR, edgecolor="white",
           linewidth=1.2, zorder=3, label="scrambled toxin (n=6, negative control)")

for name, d in real.items():
    ax.errorbar(d["iptm"][0], d["conv"][0], xerr=d["iptm"][1], yerr=d["conv"][1],
                fmt="o", ms=11, color=COLOR[name], ecolor=COLOR[name], elinewidth=1.3,
                capsize=3, markeredgecolor="white", markeredgewidth=1.2, zorder=4)
    ax.annotate(name, (d["iptm"][0], d["conv"][0]), xytext=(8, 6),
                textcoords="offset points", fontsize=9.5, fontweight="bold",
                color=COLOR[name])

ax.set_xlim(0.855, 0.935)
ax.set_ylim(0, 92)
ax.set_xlabel("best ipTM (higher = more \"confident\")")
ax.set_ylabel("ensemble pose convergence (% of 50 samples)")
# Title removed - the finding it stated ("ipTM cannot separate real complexes from
# scrambled negatives, but convergence can") is already the Figure 6 caption verbatim;
# keeping it here too was pure duplication.
# Legend anchored below the axes (axes-fraction coords, not figure coords) rather than
# inside the plot at loc="lower left", which risked overlapping the scrambled-null points.
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False, fontsize=9,
          handles=[plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=NULL_COLOR,
                               markersize=9, label='scrambled toxin (n=6, negative control)'),
                   plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#52514e",
                              markersize=9, label='real complex (mean ± sd, 4 seeds)')])

fig.tight_layout()
os.makedirs(FIGDIR, exist_ok=True)
fig.savefig(f"{FIGDIR}/fig_cofolding_convergence_vs_iptm.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{FIGDIR}/fig_cofolding_convergence_vs_iptm.pdf", bbox_inches="tight")
print("wrote fig_cofolding_convergence_vs_iptm.{png,pdf}")
print(f"null max = {null_max:.1f}%")
