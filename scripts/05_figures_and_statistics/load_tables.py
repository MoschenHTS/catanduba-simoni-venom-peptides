#!/usr/bin/env python3
"""
Shared loader for the derived tables produced by 01_collect_data.py.

Exists so 02_stats_test.py, 03_power_analysis.py, 04_make_figures.py and 05_make_hbond_figure.py
all read the SAME numbers from the SAME files, instead of each carrying its own
hand-transcribed copy (which is how two display-cutoff artifacts and three
independent copies of the MM-GBSA values ended up in the tree).

If the tables are missing, run:  python 01_collect_data.py
"""
import os
import csv
from collections import defaultdict

TABLES = os.environ.get("CSIMONI_TABLES", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "results"))

# Fixed presentation order, used by every downstream figure/table.
COMPLEXES = ["CsTx2-Kv4.2", "CsTx3-Kv4.2", "CsTx2-Kv4.3", "CsTx3-Kv4.3"]


def _load(filename, column):
    """-> {complex_label: [rep1..rep8 values]}, ordered by replica number."""
    path = os.path.join(TABLES, filename)
    if not os.path.exists(path):
        raise SystemExit(
            f"missing derived table: {path}\nRun 'python 01_collect_data.py' first.")
    per = defaultdict(dict)
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            per[row["complex"]][int(row["replica"])] = float(row[column])
    return {cx: [per[cx][r] for r in sorted(per[cx])] for cx in COMPLEXES}


def rmsd_last50():
    """Mean ligand RMSD (A) over the last 50 ns of each replica."""
    return _load("md_rmsd_paddle_contact.csv", "rmsd_last50ns_A")


def mmgbsa_total():
    """MM-GBSA dTOTAL (kcal/mol). Enthalpy-only, relative ranking ONLY - not an
    absolute binding affinity."""
    return _load("md_mmgbsa_binding_energy.csv", "dTOTAL")


def trp4_thr22_backbone():
    """Trp4(N)-Thr22(O) BACKBONE H-bond occupancy (%). Distinct from the side-chain
    OG1 variant - see trp4_thr22_sidechain()."""
    return _load("md_hbond_occupancy.csv","trp4_thr22_backbone_O_pct")


def trp4_thr22_sidechain():
    """Trp4(N)-Thr22(OG1) SIDE-CHAIN H-bond occupancy (%)."""
    return _load("md_hbond_occupancy.csv","trp4_thr22_sidechain_OG1_pct")
