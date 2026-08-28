#!/usr/bin/env python3
"""
Orient a CsTx-Kv4 complex for membrane insertion and position it in the bilayer.
- membrane normal = 1st principal axis of the VSD TM-helix-core CA cloud -> rotate to z.
- flip so the TOXIN (extracellular gating-modifier) sits on +z.
- centre xy on the bilayer centre; put the VSD TM-core CA centre at the bilayer midplane.
Writes <system>/prot_oriented.gro and prints diagnostics to validate BEFORE inserting.

Usage: 02_orient_place.py <system_dir> <n_tox_res> <vsd_type: Kv42|Kv43>
"""
import numpy as np, sys

SYS, N_TOX_RES, VSD_TYPE = sys.argv[1], int(sys.argv[2]), sys.argv[3]
GRO = f"{SYS}/prot.gro"
OUT = f"{SYS}/prot_oriented.gro"
BOX_XY = 7.962        # bilayer box (nm)
MID_Z = 4.19          # bilayer midplane (nm)

# local (1-indexed within VSD) S1-S4 TM-helix boundaries, from UniProt TRANSMEM
# annotations (Kv4.2 Q9NZV8 res185-307; Kv4.3 Q9UK17 res183-306).
LOCAL_TM = {
    "Kv42": [(1, 22), (43, 65), (73, 97), (104, 123)],
    "Kv43": [(1, 22), (42, 64), (72, 95), (102, 124)],
}[VSD_TYPE]
TM_RANGES = [(N_TOX_RES + a - 1, N_TOX_RES + b - 1) for a, b in LOCAL_TM]

lines = open(GRO).read().splitlines()
title, natom = lines[0], int(lines[1])
atoms = lines[2:2+natom]
box = lines[2+natom]

resid_seq = []
xyz = np.zeros((natom, 3)); name = []; prev = None; ridx = -1
for i, l in enumerate(atoms):
    rnum = l[0:5]
    if rnum != prev:
        ridx += 1; prev = rnum
    resid_seq.append(ridx)
    name.append(l[10:15].strip())
    xyz[i] = [float(l[20:28]), float(l[28:36]), float(l[36:44])]
resid_seq = np.array(resid_seq)
is_ca = np.array([n == "CA" for n in name])
is_tox = resid_seq < N_TOX_RES
tm_res = np.zeros(natom, bool)
for a, b in TM_RANGES:
    tm_res |= (resid_seq >= a) & (resid_seq <= b)
vsd_ca = is_ca & tm_res

P = xyz[vsd_ca]; Pc = P - P.mean(0)
evals, evecs = np.linalg.eigh(Pc.T @ Pc)
normal = evecs[:, np.argmax(evals)]
z = np.array([0, 0, 1.0]); v = np.cross(normal, z); s = np.linalg.norm(v); c = normal @ z
if s < 1e-8:
    R = np.eye(3) if c > 0 else np.diag([1, -1, -1.0])
else:
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.eye(3) + vx + vx @ vx * ((1 - c) / s**2)
xyz = xyz @ R.T

if xyz[is_tox].mean(0)[2] < xyz[vsd_ca].mean(0)[2]:
    xyz = xyz @ np.diag([1, -1, -1.0]).T
    print(">> flipped so toxin is extracellular (+z)")

shift = np.array([BOX_XY/2, BOX_XY/2, MID_Z]) - xyz[vsd_ca].mean(0)
xyz = xyz + shift

vz = xyz[vsd_ca][:, 2]; tz = xyz[is_tox][:, 2]
print(f"[{SYS}] VSD TM-core z-span : {vz.min():.2f} .. {vz.max():.2f} nm (span {vz.max()-vz.min():.2f}; want ~3-4)")
print(f"[{SYS}] toxin z-range      : {tz.min():.2f} .. {tz.max():.2f} nm  COM {tz.mean():.2f}")
print(f"[{SYS}] protein z-extent   : {xyz[:,2].min():.2f} .. {xyz[:,2].max():.2f} nm")
print(f"[{SYS}] xy extent          : x {xyz[:,0].min():.2f}-{xyz[:,0].max():.2f}, y {xyz[:,1].min():.2f}-{xyz[:,1].max():.2f} (box {BOX_XY})")

with open(OUT, "w") as f:
    f.write(title + "\n" + f"{natom}\n")
    for i, l in enumerate(atoms):
        f.write(f"{l[:20]}{xyz[i,0]:8.3f}{xyz[i,1]:8.3f}{xyz[i,2]:8.3f}\n")
    f.write(box + "\n")
print("wrote", OUT)
