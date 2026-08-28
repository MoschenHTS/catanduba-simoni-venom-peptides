#!/usr/bin/env python3
"""
Step 2 — Figure 2: multiple sequence alignment of the four C. simoni
peptides against their closest theraphotoxin homologs in UniProtKB/Swiss-Prot.

Pipeline:
  1. Reference panel = the top-10 distinct homologs ranked by pairwise
     % identity to CsTx3 over ungapped columns (hits originally parsed from
     generated/blast_results/*.xml, produced by 01_blast_homolog_search.py;
     mature chains from UniProt, isoform duplicates collapsed).
  2. MAFFT L-INS-i alignment (BLOSUM62).
  3. Property-coloured cells, ICK cysteine framework (CI-CVI) with
     disulfide connectivity, accessions, source species, and % identity
     to CsTx3.
  4. CsTx1's C-terminus (~223 Da) is unresolved -> rendered as a hatched
     block, excluded from the identity calculation.

Requires MAFFT on PATH (https://mafft.cbrc.jp/alignment/software/).

Run (from the repository root):
    python scripts/01_sequence_analysis/02_msa_alignment_figure.py

Output:
    generated/figures/Figure2_alignment.png / .pdf
    generated/blast_results/figure2_alignment.fasta
"""
import subprocess
import tempfile
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from Bio.Align import PairwiseAligner, substitution_matrices

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,  # editable text in vector output
})

# ---------------------------------------------------------------- % identity
CSTX3 = "GCRWMFGACKTTADCCKALACVGTCIWDGTYGN"
_pa = PairwiseAligner()
_pa.substitution_matrix = substitution_matrices.load("BLOSUM62")
_pa.open_gap_score = -10
_pa.extend_gap_score = -0.5
_pa.mode = "global"


def pairwise_identity(seq, ref=CSTX3):
    """Pairwise %% identity to `ref`, over ungapped aligned columns only."""
    a = _pa.align(ref, seq)[0]
    ta, tb = str(a[0]), str(a[1])
    identical = sum(1 for x, y in zip(ta, tb) if x == y and x != "-")
    ungapped = sum(1 for x, y in zip(ta, tb) if x != "-" and y != "-")
    return 100 * identical / ungapped if ungapped else 0


def abbreviate_species(species):
    parts = species.split()
    return parts[0][0] + ". " + " ".join(parts[1:]) if len(parts) > 1 else species


# ---------------------------------------------------------------------- data
QUERIES = [  # name, sequence, species, C-terminus unresolved
    ("Cs32.1", "GCRWMFGACKTTADCCKALACVGTCIWDGS",    "Catanduba simoni", False),
    ("CsTx1",  "GCRWMFGACKTTADCCKALACVGTCIWDG",     "Catanduba simoni", True),
    ("CsTx2",  "GCRWMFGACKTTADCCKALACVGTCIWDGTYG",  "Catanduba simoni", False),
    ("CsTx3",  "GCRWMFGACKTTADCCKALACVGTCIWDGTYGN", "Catanduba simoni", False),
]

REFERENCE_HOMOLOGS = [  # name, UniProt accession, mature sequence, species
    ("SGTx1",   "P56855", "TCRYLFGGCKTTADCCKHLACRSDGKYCAWDGTF",    "Stromatopelma calceatum griseipes"),
    ("HmTx2",   "P0DOC5", "ECRYLFGGCKTTADCCKHLGCRTDLYYCAWDGTF",    "Heteroscodra maculata"),
    ("JZTX-42", "B1P1A2", "ECRWMFGGCTTDSDCCEHLGCRWEKPSWCAWDGTVRK", "Chilobrachys guangxiensis"),
    ("HaTx2",   "P56853", "ECRYLFGGCKTTADCCKHLGCKFRDKYCAWDFTFS",   "Grammostola rosea"),
    ("VaTx1",   "P0C244", "SECRWFMGGCDSTLDCCKHLSCKMGLYYCAWDGTF",   "Psalmopoeus cambridgei"),
    ("GxTX-2",  "P84837", "ECRKMFGGCSVDSDCCAHLGCKPTLKYCAWDGT",     "Chilobrachys guangxiensis"),
    ("ScTx1",   "P60991", "DCTRMFGACRRDSDCCPHLGCKPTSKYCAWDGTI",    "Stromatopelma calceatum"),
    ("Eo1a",    "P0DW95", "DCRWFLGGCSKDSDCCKHLACRIDGYIKYCAWDGTF",  "Encyocratella olivacea"),
    ("Pmu1a",   "P0DQZ2", "ECRWFWGGCNNDADCCKHLECKRKWPHICLWDGTFT",  "Pterinochilus murinus"),
    ("HaTx1",   "P56852", "ECRYLFGGCKTTSDCCKHLGCKFRDKYCAWDFTFS",   "Grammostola rosea"),
]
REFERENCE_HOMOLOGS.sort(key=lambda r: pairwise_identity(r[2]), reverse=True)  # objective rank order

# ----------------------------------------------------------------- alignment
def run_mafft(named_seqs):
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as f:
        for name, s in named_seqs:
            f.write(f">{name}\n{s}\n")
        path = f.name
    out = subprocess.run(
        ["mafft", "--localpair", "--maxiterate", "1000", "--quiet", path],
        capture_output=True, text=True,
    ).stdout
    os.unlink(path)
    alignment, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            alignment[name] = ""
        elif name:
            alignment[name] += line.strip()
    return alignment


names = [q[0] for q in QUERIES] + [r[0] for r in REFERENCE_HOMOLOGS]
seq_by_name = {q[0]: q[1] for q in QUERIES} | {r[0]: r[2] for r in REFERENCE_HOMOLOGS}
species_by_name = {q[0]: q[2] for q in QUERIES} | {r[0]: r[3] for r in REFERENCE_HOMOLOGS}
accession_by_name = {r[0]: r[1] for r in REFERENCE_HOMOLOGS}
unresolved_by_name = {q[0]: q[3] for q in QUERIES}

alignment = run_mafft([(n, seq_by_name[n]) for n in names])
aln_len = len(next(iter(alignment.values())))

# ------------------------------------------------------- cysteine framework
cys_columns = [i for i in range(aln_len) if alignment["CsTx3"][i] == "C"]
cstx1_last_col = max(i for i in range(aln_len) if alignment["CsTx1"][i] != "-")
query_last_col = max(
    i for n in ["Cs32.1", "CsTx2", "CsTx3"]
    for i in range(aln_len) if alignment[n][i] != "-"
)

# -------------------------------------------------- ClustalW-style consensus
STRONG_GROUPS = ["STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"]
WEAK_GROUPS = ["CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK",
               "NEQHRK", "FVLIM", "HFY"]


def consensus_symbol(column):
    residues = [c for c in column if c != "-"]
    if len(residues) < len(column):
        return " "
    unique = set(residues)
    if len(unique) == 1:
        return "*"
    if any(unique <= set(g) for g in STRONG_GROUPS):
        return ":"
    if any(unique <= set(g) for g in WEAK_GROUPS):
        return "."
    return " "


consensus = [
    consensus_symbol("".join(alignment[n][i] for n in names))
    for i in range(aln_len)
]

# --------------------------------------------------------- residue colours
PROPERTY_COLOR = {
    **{c: "#E6B800" for c in "C"},        # cysteine    - gold
    **{c: "#2C6FBF" for c in "KRH"},      # positive    - blue
    **{c: "#CC3333" for c in "DE"},       # negative    - red
    **{c: "#33994D" for c in "FYWILMV"},  # hydrophobic - green
    **{c: "#E08A1E" for c in "STNQ"},     # polar       - orange
    **{c: "#9AA0A6" for c in "AG"},       # small       - grey
}


def residue_color(c):
    return PROPERTY_COLOR.get(c, "#FFFFFF")


# ------------------------------------------------------------------ layout
CELL_W, CELL_H = 0.42, 0.50
ROW_GAP = 0.12
ROW_PITCH = CELL_H + ROW_GAP
BLOCK_GAP = 0.42  # extra gap between query block and reference block
n_queries, n_rows = len(QUERIES), len(names)


def row_y(r):
    y = (n_rows - 1 - r) * ROW_PITCH
    return y + BLOCK_GAP if r < n_queries else y


align_top = (n_rows - 1) * ROW_PITCH + BLOCK_GAP + CELL_H
ref_top = (n_rows - n_queries) * ROW_PITCH - ROW_GAP
gap_mid = ref_top + (BLOCK_GAP + ROW_GAP) / 2

x_left, x_right = -2.8, aln_len * CELL_W + 5.6
y_bottom, y_top = -1.95, align_top + 1.45
scale = 6.85 / (x_right - x_left)  # 174 mm = Springer full-page (double-column) width
fig, ax = plt.subplots(figsize=((x_right - x_left) * scale, (y_top - y_bottom) * scale))
ax.set_xlim(x_left, x_right)
ax.set_ylim(y_bottom, y_top)
ax.set_aspect("equal")
ax.axis("off")

# cysteine numerals above the alignment; disulfide arcs above the numerals
for k, ci in enumerate(cys_columns):
    ax.text(ci * CELL_W + CELL_W / 2, align_top + 0.08, ["I", "II", "III", "IV", "V", "VI"][k],
            ha="center", va="bottom", fontsize=6, color="#8a6d00", fontweight="bold")
arc_top = align_top + 0.50
for a, b in [(0, 3), (1, 4), (2, 5)]:
    xa = cys_columns[a] * CELL_W + CELL_W / 2
    xb = cys_columns[b] * CELL_W + CELL_W / 2
    ax.add_patch(mpatches.FancyArrowPatch(
        (xa, arc_top), (xb, arc_top), connectionstyle="arc3,rad=-0.11",
        arrowstyle="-", lw=1.0, color="#8a6d00"))

# column headers
header_y = align_top + 0.30
ax.text(-0.3, header_y, "Toxin", fontsize=7, fontweight="bold", va="center", ha="right")
for x, label in [(aln_len * CELL_W + 0.3, "Acc."),
                  (aln_len * CELL_W + 2.3, "% id"),
                  (aln_len * CELL_W + 3.6, "Species")]:
    ax.text(x, header_y, label, fontsize=7, fontweight="bold", va="center",
            style="italic" if label == "Species" else "normal")

# alignment cells
for r, name in enumerate(names):
    y = row_y(r)
    seq = alignment[name]
    for i in range(aln_len):
        if name == "CsTx1" and cstx1_last_col < i <= query_last_col:  # unresolved block
            ax.add_patch(Rectangle((i * CELL_W, y), CELL_W, CELL_H, facecolor="#E5E5E5",
                                    edgecolor="#B0B0B0", hatch="////", lw=0.4))
            continue
        c = seq[i]
        if c == "-":
            continue
        ax.add_patch(Rectangle((i * CELL_W, y), CELL_W, CELL_H, facecolor=residue_color(c),
                                edgecolor="white", lw=0.5))
        ax.text(i * CELL_W + CELL_W / 2, y + CELL_H / 2, c, ha="center", va="center",
                fontsize=5.8, color="white", fontweight="bold")

# row labels + right-hand columns
for r, name in enumerate(names):
    y = row_y(r) + CELL_H / 2
    is_query = name in unresolved_by_name
    ax.text(-0.3, y, name, fontsize=6.8, va="center", ha="right",
            fontweight="bold" if is_query else "normal")
    ax.text(aln_len * CELL_W + 3.6, y, abbreviate_species(species_by_name[name]),
            fontsize=6.0, va="center", style="italic", color="#222")
    if not is_query:
        ax.text(aln_len * CELL_W + 0.3, y, accession_by_name[name], fontsize=6, va="center", color="#333")
        ax.text(aln_len * CELL_W + 2.3, y, f"{pairwise_identity(seq_by_name[name]):.1f}",
                fontsize=6.2, va="center")

# consensus symbol row (symbol key belongs in the manuscript caption)
consensus_y = -0.42
ax.text(-0.3, consensus_y, "conservation", fontsize=6.0, style="italic", color="#666",
        va="center", ha="right")
for i in range(aln_len):
    if consensus[i] != " ":
        ax.text(i * CELL_W + CELL_W / 2, consensus_y, consensus[i], ha="center", va="center",
                fontsize=7, fontweight="bold", color="#333")

# divider between query block and reference block
ax.plot([-2.6, aln_len * CELL_W + 5.3], [gap_mid, gap_mid], lw=0.7, color="#999", clip_on=False)

# colour legend
legend_row1 = [("Cysteine (C)", "#E6B800"), ("Positive (K,R,H)", "#2C6FBF"),
               ("Negative (D,E)", "#CC3333"), ("Hydrophobic (F,Y,W,I,L,M,V)", "#33994D")]
legend_row2 = [("Polar (S,T,N,Q)", "#E08A1E"), ("Small (A,G)", "#9AA0A6"),
               ("Unresolved", "#E5E5E5")]
char_w, pad = 0.17, 1.2
col_max_len = [max(len(legend_row1[c][0]), len(legend_row2[c][0]) if c < len(legend_row2) else 0)
               for c in range(4)]
col_x = [0.0]
for c in range(3):
    col_x.append(col_x[-1] + 0.62 + col_max_len[c] * char_w + pad)
for ly, row in [(-1.10, legend_row1), (-1.68, legend_row2)]:
    for ci, (label, color) in enumerate(row):
        x = col_x[ci]
        is_unresolved = label == "Unresolved"
        ax.add_patch(Rectangle((x, ly), 0.46, 0.40, facecolor=color,
                                edgecolor="#999" if is_unresolved else "white",
                                hatch="////" if is_unresolved else None, lw=0.5))
        ax.text(x + 0.62, ly + 0.20, label, fontsize=6.0, va="center")

# ------------------------------------------------------------------- output
os.makedirs("generated/figures", exist_ok=True)
plt.savefig("generated/figures/Figure2_alignment.png", dpi=600, bbox_inches="tight", pad_inches=0.03)
plt.savefig("generated/figures/Figure2_alignment.pdf", bbox_inches="tight", pad_inches=0.03)
w_in, h_in = fig.get_size_inches()
print(f"wrote generated/figures/Figure2_alignment.png/.pdf  (canvas {w_in * 25.4:.0f} x {h_in * 25.4:.0f} mm)")

os.makedirs("generated/blast_results", exist_ok=True)
with open("generated/blast_results/figure2_alignment.fasta", "w") as f:
    for n in names:
        f.write(f">{n}\n{alignment[n]}\n")

print("\n% identity to CsTx3 (pairwise, ungapped):")
for r in REFERENCE_HOMOLOGS:
    print(f"  {r[0]:8} {r[1]:8} {pairwise_identity(r[2]):5.1f}%")
