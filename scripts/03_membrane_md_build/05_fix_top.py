#!/usr/bin/env python3
"""Deterministically rewrite the [ molecules ] block from a .gro (contiguous
protein/POPC/SOL ordering, before ions are added). Usage: 05_fix_top.py gro top"""
import sys
gro, top = sys.argv[1:3]
L = open(gro).read().splitlines(); n = int(L[1]); atoms = L[2:2+n]
resn = [a[5:10].strip() for a in atoms]
nP = sum(r == "POPC" for r in resn) // 134
nS = sum(r == "SOL" for r in resn) // 3
t = open(top).read().splitlines(); out = []
for line in t:
    out.append(line)
    if line.strip().replace(" ", "").lower() == "[molecules]":
        out += ["; Compound        #mols", "Protein_chain_A     1",
                f"POPC              {nP}", f"SOL               {nS}"]
        break
open(top, "w").write("\n".join(out) + "\n")
print(f"fix_top: Protein 1, POPC {nP}, SOL {nS}")
