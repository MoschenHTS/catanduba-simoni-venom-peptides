#!/usr/bin/env python3
"""Delete waters that gmx solvate placed inside the lipid hydrophobic core
(|z-midplane| < CORE_HALF) and that are NOT lining the protein (crevice water,
kept if within KEEP_NEAR of a protein atom). Updates the .top SOL count.
Usage: 04_clean_core_water.py in.gro out.gro system.top n_prot_atoms"""
import sys, numpy as np
GRO_IN, GRO_OUT, TOP = sys.argv[1:4]
NPROT = int(sys.argv[4])
MID = 4.19; CORE_HALF = 1.4; KEEP_NEAR = 0.5

L = open(GRO_IN).read().splitlines()
n = int(L[1]); atoms = L[2:2+n]; box = L[2+n]
xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms])
resn = [a[5:10].strip() for a in atoms]
prot = xyz[:NPROT]

# SOL molecules = 3 consecutive atoms (OW, HW1, HW2), contiguous block(s)
sol_idx = [i for i in range(n) if resn[i] == "SOL"]
assert len(sol_idx) % 3 == 0
remove = set()
kept_sol = 0
for k in range(0, len(sol_idx), 3):
    ow = sol_idx[k]
    z = xyz[ow, 2]
    if abs(z - MID) < CORE_HALF:
        d = np.sqrt(((xyz[ow] - prot)**2).sum(1)).min()
        if d > KEEP_NEAR:
            remove.update(sol_idx[k:k+3]); continue
    kept_sol += 1

keep = [i for i in range(n) if i not in remove]
with open(GRO_OUT, "w") as f:
    f.write(L[0] + "\n%d\n" % len(keep))
    for i in keep:
        f.write(atoms[i] + "\n")
    f.write(box + "\n")
print(f"core-water cleanup: removed {len(remove)//3} waters, kept {kept_sol} SOL")

# update .top SOL count (last SOL line)
t = open(TOP).read().splitlines()
for i in range(len(t)-1, -1, -1):
    if t[i].split() and t[i].split()[0] == "SOL":
        t[i] = f"SOL               {kept_sol}"; break
open(TOP, "w").write("\n".join(t) + "\n")
print("updated", TOP, "SOL ->", kept_sol)
