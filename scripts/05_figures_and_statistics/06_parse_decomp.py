#!/usr/bin/env python3
"""
Parse FINAL_DECOMP_MMPBSA.dat (Total Energy Decomposition section) for all 4
complexes x 8 replicas (extended 2026-08 from n=3), average across replicas, rank
residues by NON-BONDED
contribution (vdW+Elec+PolarSolv+NonPolarSolv - excluding the noisy per-residue
'Internal'/bonded term, which is only meaningful summed over the whole complex,
not attributable to individual residues).
"""
import os
import glob, re, numpy as np
from collections import defaultdict

# Project root. Set CSIMONI_BASE to point at wherever the primary MD data (not
# shipped in this repository - see the top-level README) has been placed.
CSIMONI_BASE = os.environ.get("CSIMONI_BASE", "/path/to/this/repo")

# CSIMONI_MD_BASE overrides the directory holding <system>/ data directly (used by
# the reorganized repo, where it is not literally "<root>/md"); unset, falls back
# to CSIMONI_BASE/md, matching the original working-tree layout.
BASE = os.environ.get("CSIMONI_MD_BASE", os.path.join(CSIMONI_BASE, "md"))
COMPLEXES = ["cstx2_kv42", "cstx3_kv42", "cstx2_kv43", "cstx3_kv43"]


def parse_file(path):
    """Return {residue_label: nonbonded_total} for the 'Total Energy
    Decomposition' section (skip the later Sidechain/Backbone sections)."""
    lines = open(path).read().splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("DELTAS:"))
    out = {}
    for l in lines[start:]:
        if l.startswith("Sidechain Energy") or l.startswith("Backbone Energy"):
            break
        if not l.startswith(("L:", "R:")):
            continue
        parts = l.split(",")
        label = parts[0]  # e.g. "L:A:TRP:4"
        vals = list(map(float, parts[1:]))
        # columns: Internal(avg,sd,sem) vdW(avg,sd,sem) Elec(avg,sd,sem)
        #          PolarSolv(avg,sd,sem) NonPolarSolv(avg,sd,sem) TOTAL(avg,sd,sem)
        vdw, elec, polar, nonpolar = vals[3], vals[6], vals[9], vals[12]
        out[label] = vdw + elec + polar + nonpolar  # non-bonded total (excl. Internal)
    return out


for cx in COMPLEXES:
    per_rep = []
    for rep in (1, 2, 3, 4, 5, 6, 7, 8):
        f = f"{BASE}/{cx}/mmpbsa/rep{rep}_FINAL_DECOMP_MMPBSA.dat"
        try:
            per_rep.append(parse_file(f))
        except FileNotFoundError:
            print(f"[{cx}] rep{rep} decomp file missing, skip")
    if not per_rep:
        continue
    all_labels = set().union(*[set(d) for d in per_rep])
    avg = {lbl: np.mean([d.get(lbl, np.nan) for d in per_rep]) for lbl in all_labels}
    # split toxin (L:) vs VSD (R:)
    tox = sorted([(l, v) for l, v in avg.items() if l.startswith("L:")], key=lambda x: x[1])
    vsd = sorted([(l, v) for l, v in avg.items() if l.startswith("R:")], key=lambda x: x[1])
    print(f"\n{'='*70}\n{cx}  (non-bonded per-residue dG, kcal/mol, mean of 8 reps; more negative = more favorable)\n{'='*70}")
    print("  TOXIN residues (most favorable first):")
    for l, v in tox[:6]:
        _, ch, resn, resi = l.split(":")
        print(f"    {resn}{resi:<4} {v:+7.2f}")
    print("  VSD residues (most favorable first):")
    for l, v in vsd[:6]:
        _, ch, resn, resi = l.split(":")
        print(f"    {resn}{resi:<4} {v:+7.2f}")
    print("  Most UNFAVORABLE (either side):")
    worst = sorted(avg.items(), key=lambda x: -x[1])[:4]
    for l, v in worst:
        _, ch, resn, resi = l.split(":")
        side = "TOX" if l.startswith("L:") else "VSD"
        print(f"    [{side}] {resn}{resi:<4} {v:+7.2f}")
