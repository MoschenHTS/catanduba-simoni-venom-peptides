#!/usr/bin/env python3
"""Write index.ndx with the two coupling groups the mdps expect:
   MEMBRANE = protein + POPC ; SOLV = water + ions (everything else).
Also emits Protein / POPC / SOLU convenience groups.
Usage: 06_make_index.py system.gro index.ndx n_prot_atoms"""
import sys
GRO, NDX = sys.argv[1:3]
NPROT = int(sys.argv[3])
L = open(GRO).read().splitlines()
n = int(L[1]); atoms = L[2:2+n]
resn = [a[5:10].strip() for a in atoms]

groups = {"System": list(range(1, n+1)),
          "Protein": list(range(1, NPROT+1)),
          "POPC": [i+1 for i in range(n) if resn[i] == "POPC"],
          "MEMBRANE": [i+1 for i in range(n) if i < NPROT or resn[i] == "POPC"],
          "SOLV": [i+1 for i in range(n) if i >= NPROT and resn[i] != "POPC"]}
with open(NDX, "w") as f:
    for name, idx in groups.items():
        f.write(f"[ {name} ]\n")
        for j in range(0, len(idx), 15):
            f.write(" ".join(map(str, idx[j:j+15])) + "\n")
print("index.ndx groups:", {k: len(v) for k, v in groups.items()})
assert len(groups["MEMBRANE"]) + len(groups["SOLV"]) == n, "groups must partition system"
