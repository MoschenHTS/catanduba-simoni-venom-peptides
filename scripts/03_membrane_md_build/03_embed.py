#!/usr/bin/env python3
"""
Insert the oriented protein into the equilibrated POPC bilayer:
 - remove POPC molecules that clash with the protein (any lipid atom < CLASH_LIP)
 - remove bilayer waters that overlap the protein (< CLASH_WAT)
 - box z extended to BOXZ nm for water headroom above the toxin
Writes <system>/embedded.gro + <system>/system.top (protein + kept POPC + kept SOL).

Usage: 03_embed.py <system_dir> <n_prot_atoms>
"""
import sys
import numpy as np

SYS, NPROT = sys.argv[1], int(sys.argv[2])
PROT = f"{SYS}/prot_oriented.gro"
BILA = "bilayer/start.gro"
OUTG = f"{SYS}/embedded.gro"
OUTT = f"{SYS}/system.top"
BOXZ = 10.0
CLASH_LIP = 0.35
CLASH_WAT = 0.28
POPC_N = 134
SOL_N = 3


def read_gro(path):
    L = open(path).read().splitlines()
    n = int(L[1]); atoms = L[2:2 + n]
    xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms])
    return L[0], atoms, xyz, L[2 + n]


_, patoms, pxyz, _ = read_gro(PROT)
_, batoms, bxyz, bbox = read_gro(BILA)
assert len(patoms) == NPROT, f"expected {NPROT} protein atoms, prot_oriented.gro has {len(patoms)}"

resn = [a[5:10].strip() for a in batoms]
popc_idx = [i for i, r in enumerate(resn) if r == "POPC"]
sol_idx = [i for i, r in enumerate(resn) if r == "SOL"]
assert len(popc_idx) % POPC_N == 0 and len(sol_idx) % SOL_N == 0

def min_dist_to_prot(coords):
    out = np.empty(len(coords))
    for s in range(0, len(coords), 2000):
        c = coords[s:s+2000]
        d = np.sqrt(((c[:, None, :] - pxyz[None, :, :])**2).sum(-1)).min(1)
        out[s:s+2000] = d
    return out

keepP = []
pc = bxyz[popc_idx]
pd = min_dist_to_prot(pc).reshape(-1, POPC_N)
for m in range(pd.shape[0]):
    if pd[m].min() >= CLASH_LIP:
        keepP.extend(popc_idx[m*POPC_N:(m+1)*POPC_N])
nP = len(keepP) // POPC_N

keepS = []
sc = bxyz[sol_idx]
sd = min_dist_to_prot(sc).reshape(-1, SOL_N)
for m in range(sd.shape[0]):
    if sd[m].min() >= CLASH_WAT:
        keepS.extend(sol_idx[m*SOL_N:(m+1)*SOL_N])
nS = len(keepS) // SOL_N

print(f"[{SYS}] POPC: kept {nP}/200 (removed {200-nP})")
print(f"[{SYS}] SOL : kept {nS}/{len(sol_idx)//SOL_N} (removed {len(sol_idx)//SOL_N - nS})")

out = patoms + [batoms[i] for i in keepP] + [batoms[i] for i in keepS]
bx = bbox.split()
with open(OUTG, "w") as f:
    f.write(f"{SYS} in POPC\n%d\n" % len(out))
    for a in out:
        f.write(a + "\n")
    f.write(f"{float(bx[0]):10.5f}{float(bx[1]):10.5f}{BOXZ:10.5f}\n")
print("wrote", OUTG, "(", len(out), "atoms )")

top = open(f"{SYS}/topol.top").read()
# Relative include: grompp resolves #include against the directory of the including
# .top file, so "../bilayer/POPC.itp" works from md/<system>/system.top regardless of
# where the project tree lives. This previously wrote an ABSOLUTE path, which made the
# generated topologies non-portable and - worse - silently readable from a stale copy
# of the tree if one existed at the old location. NOTE: pdb2gmx separately bakes three
# absolute includes (forcefield.itp, tip3p.itp, ions.itp) into topol.top; the committed
# topologies therefore still contain absolute paths. Run 13_normalize_top_paths.py to
# convert them before using these topologies on another machine.
popc_inc = '#include "../bilayer/POPC.itp"\n'
top = top.replace("[ system ]", popc_inc + "\n[ system ]", 1)
top = top.rstrip() + f"\nPOPC              {nP}\nSOL               {nS}\n"
open(OUTT, "w").write(top)
print("wrote", OUTT, "-> [molecules]: Protein_chain_A 1 | POPC", nP, "| SOL", nS)
