#!/usr/bin/env python3
"""
Publication figures for the membrane MD campaign.
Fig 1: ligand RMSD vs time, 4 complexes x 8 replicas (2x2 small multiples).
Fig 2: MM-GBSA binding-energy comparison, 4 complexes, replica points + mean.
Colors: fixed categorical order (validated CVD-safe set, dataviz skill reference
palette), one color per complex, consistent across both figures.
Extended 2026-08 from 3 -> 8 replicas/complex for statistical power (see
03_power_analysis.py); individual per-replica legend dropped (too many series),
replaced by a bold mean trace over the 8 thin per-replica lines.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from load_tables import mmgbsa_total

# Project root. Set CSIMONI_BASE to point at wherever the primary MD data (not
# shipped in this repository - see the top-level README) has been placed.
CSIMONI_BASE = os.environ.get("CSIMONI_BASE", "/path/to/this/repo")

# CSIMONI_MD_BASE overrides the directory holding <system>/ data directly; unset,
# falls back to CSIMONI_BASE/md.
BASE = os.environ.get("CSIMONI_MD_BASE", os.path.join(CSIMONI_BASE, "md"))

# Where figures are written. Defaults to generated/figures/ at the repository root;
# override with CSIMONI_FIGDIR to write elsewhere.
FIGDIR = os.environ.get("CSIMONI_FIGDIR", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "generated", "figures"))
COMPLEXES = ["cstx2_kv42", "cstx3_kv42", "cstx2_kv43", "cstx3_kv43"]
LABELS = {"cstx2_kv42": "CsTx2–Kv4.2", "cstx3_kv42": "CsTx3–Kv4.2",
          "cstx2_kv43": "CsTx2–Kv4.3", "cstx3_kv43": "CsTx3–Kv4.3"}
COLOR = {"cstx2_kv42": "#2a78d6", "cstx2_kv43": "#1baf7a",
         "cstx3_kv42": "#eda100", "cstx3_kv43": "#008300"}

mpl.rcParams.update({
    "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e5e4df", "grid.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

# ---------- Figure 1: RMSD vs time, 2x2 small multiples ----------
fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)

for ax, sysname in zip(axes.flat, COMPLEXES):
    c = COLOR[sysname]
    traces = []
    for rep in range(1, 9):
        f = f"{BASE}/{sysname}/prod_rep{rep}_rmsd.dat"
        d = np.loadtxt(f)
        t, rmsd = d[:, 0], d[:, 1]
        ax.plot(t, rmsd, color=c, alpha=0.35, lw=0.9, zorder=2)
        traces.append(rmsd)
    mean_trace = np.mean(np.vstack(traces), axis=0)
    ax.plot(t, mean_trace, color=c, alpha=1.0, lw=2.0, zorder=3, label="8-replica mean")
    # Small in-axes corner label (not a title bar) identifies each complex; thin/bold
    # line meaning and the 8-replica-mean legend now live in the Figure 7 caption text
    # instead of a per-subplot legend repeated identically 4 times.
    ax.text(0.03, 0.95, LABELS[sysname], transform=ax.transAxes, fontsize=10,
            fontweight="bold", color=c, va="top", ha="left")
    ax.set_ylim(0, 18)

for ax in axes[-1, :]:
    ax.set_xlabel("time (ns)")
for ax in axes[:, 0]:
    ax.set_ylabel("ligand RMSD (Å)")

# Suptitle removed - "Toxin pose stability during 200 ns membrane MD (8 replicas/complex)"
# duplicates the Figure 7 caption verbatim.
fig.tight_layout()
os.makedirs(FIGDIR, exist_ok=True)
fig.savefig(f"{FIGDIR}/fig_md_rmsd_stability.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{FIGDIR}/fig_md_rmsd_stability.pdf", bbox_inches="tight")
print("wrote fig_md_rmsd_stability.{png,pdf}")

# ---------- Figure 2: MM-GBSA comparison ----------
# dTOTAL per replica, read from results/md_mmgbsa_binding_energy.csv (01_collect_data.py),
# which parses mmpbsa/repN_FINAL_RESULTS_MMPBSA.dat. n=8/complex, extended 2026-08
# from n=3 for statistical power; see 03_power_analysis.py.
_gbsa_by_label = mmgbsa_total()
mmgbsa = {sysname: _gbsa_by_label[LABELS[sysname].replace("–", "-")]
          for sysname in COMPLEXES}

fig2, ax2 = plt.subplots(figsize=(6.5, 5))
x = np.arange(len(COMPLEXES))
means = [np.mean(mmgbsa[s]) for s in COMPLEXES]
sds = [np.std(mmgbsa[s], ddof=1) for s in COMPLEXES]
colors = [COLOR[s] for s in COMPLEXES]

ax2.axhline(0, color="#52514e", lw=1, zorder=1)
ax2.bar(x, means, yerr=sds, capsize=4, color=colors, alpha=0.85,
        edgecolor="white", linewidth=1.5, width=0.6, zorder=2,
        error_kw=dict(ecolor="#0b0b0b", elinewidth=1.2, capthick=1.2))
# individual replica points (show the raw spread honestly)
rng = np.random.default_rng(0)
for xi, s in zip(x, COMPLEXES):
    jitter = rng.uniform(-0.12, 0.12, size=len(mmgbsa[s]))
    ax2.scatter(xi + jitter, mmgbsa[s], color="#0b0b0b", s=20, zorder=3,
                edgecolor="white", linewidth=0.6, alpha=0.85)

ax2.set_xticks(x)
ax2.set_xticklabels([LABELS[s] for s in COMPLEXES], fontsize=10)
ax2.set_ylabel("MM-GBSA ΔG  (kcal/mol, relative ranking only)")
# Title removed - "Binding-energy comparison (single-trajectory MM-GBSA, n=8/complex)"
# duplicates the Figure 7 caption verbatim.
ax2.text(0.02, 0.02, "more favorable ↓", transform=ax2.transAxes,
          fontsize=8, color="#52514e", style="italic")
fig2.tight_layout()
fig2.savefig(f"{FIGDIR}/fig_mmgbsa_comparison.png", dpi=300, bbox_inches="tight")
fig2.savefig(f"{FIGDIR}/fig_mmgbsa_comparison.pdf", bbox_inches="tight")
print("wrote fig_mmgbsa_comparison.{png,pdf}")
