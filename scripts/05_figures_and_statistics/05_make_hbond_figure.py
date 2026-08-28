#!/usr/bin/env python3
"""
Trp4(toxin)-Thr22(VSD) BACKBONE-carbonyl H-bond occupancy, n=8/complex.

Nonzero in all 16 Kv4.3 replicas and exactly 0.00% in all 16 Kv4.2 replicas - a clean
subtype split. Two Kv4.3 replicas are only marginally nonzero (5.2% and 3.0%) and should
not be read as evidence of a persistent bond there; the other 14 sit at 27-89%.

Plots the backbone-carbonyl acceptor ONLY. The same residue pair also interacts via the
Thr22 side-chain hydroxyl (OG1) - a different interaction, correctly shown here as 0 where
it is the only one present. The side-chain variant is NOT Kv4.2-exclusive: it appears at
>=10% occupancy in 2/16 Kv4.2 replicas (28%, 42%) AND 2/16 Kv4.3 replicas (both alongside a
substantial backbone bond in the same replica), plus several more sub-10% traces on both
channels.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from load_tables import trp4_thr22_backbone

# Where the figure is written. Defaults to generated/figures/ at the repository
# root; override with CSIMONI_FIGDIR to write elsewhere.
FIGDIR = os.environ.get("CSIMONI_FIGDIR", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "generated", "figures"))
# Occupancy (%) of TOX:TRP4:N -- VSD:THR22:O (backbone carbonyl), read from
# results/md_hbond_occupancy.csv (01_collect_data.py -> per-replica *_hbonds.csv).
# Previously transcribed by hand from the batch logs, which only printed pairs at
# >=10% occupancy - that silently rounded two genuinely small values down to 0
# (CsTx2-Kv4.3 rep4 = 5.2%, CsTx3-Kv4.3 rep5 = 3.0%). Now read at full precision.
_bb = trp4_thr22_backbone()
DATA = {"CsTx2–Kv4.2": _bb["CsTx2-Kv4.2"], "CsTx3–Kv4.2": _bb["CsTx3-Kv4.2"],
        "CsTx2–Kv4.3": _bb["CsTx2-Kv4.3"], "CsTx3–Kv4.3": _bb["CsTx3-Kv4.3"]}
COLOR = {"CsTx2–Kv4.2": "#2a78d6", "CsTx2–Kv4.3": "#1baf7a",
         "CsTx3–Kv4.2": "#eda100", "CsTx3–Kv4.3": "#008300"}

mpl.rcParams.update({
    "font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
    "text.color": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e5e4df", "grid.linewidth": 0.6,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

fig, ax = plt.subplots(figsize=(6.8, 5))
names = list(DATA.keys())
x = np.arange(len(names))
means = [np.mean(DATA[n]) for n in names]
colors = [COLOR[n] for n in names]

ax.bar(x, means, color=colors, alpha=0.85, edgecolor="white", linewidth=1.5, width=0.6, zorder=2)
for xi, n in zip(x, names):
    jitter = np.linspace(-0.1, 0.1, len(DATA[n]))
    ax.scatter(xi + jitter, DATA[n], color="#0b0b0b", s=20, zorder=3,
               edgecolor="white", linewidth=0.6, alpha=0.85)

ax.set_xticks(x)
ax.set_xticklabels(names, fontsize=10)
ax.set_ylabel("Trp4(N)–Thr22(O) H-bond occupancy (% of trajectory)")
ax.set_title("A Kv4.3-specific backbone anchor, absent from Kv4.2 (n=8/complex)\n(toxin Trp4 backbone amide to VSD Thr22 backbone carbonyl)",
              fontsize=11.5, fontweight="bold")
ax.set_ylim(-3, 100)
fig.tight_layout()
os.makedirs(FIGDIR, exist_ok=True)
fig.savefig(f"{FIGDIR}/fig_trp4_thr22_hbond.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{FIGDIR}/fig_trp4_thr22_hbond.pdf", bbox_inches="tight")
print("wrote fig_trp4_thr22_hbond.{png,pdf}")
