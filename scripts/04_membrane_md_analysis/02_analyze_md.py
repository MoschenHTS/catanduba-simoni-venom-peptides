#!/usr/bin/env python3
"""
MD pose-stability analysis for a CsTx-Kv4 replica (run under boltz_env, needs
MDAnalysis). Run from inside the system dir (e.g. md/cstx3_kv43/). Uses npt2.gro
(not the .tpr - GROMACS 2026 tpx version 138 is too new for MDAnalysis's TPR
parser) as topology, since only resid/resname/name are needed (no bonds).

IMPORTANT: resid is NOT unique across chains - it RESETS to 1 for the VSD
(pdb2gmx numbers each chain from 1). Toxin = protein.residues[:n_tox_res],
VSD = protein.residues[n_tox_res:], selected by RESIDUE INDEX, not resid.

For each replica:
  * ligand RMSD over time (toxin CA, after superposing on VSD TM-core CA) ->
    does the pose drift / dissociate, or stay put?
  * persistent contacts: % of frames the S4/paddle region is within 4.5 A of
    the toxin (gating charges included)
  * R3 - VSD acidic salt bridge occupancy over the trajectory (the Boltz-2 co-folding
    pharmacophore hit) as a time series, not just a single-model snapshot
Writes a plain-text report + a time-series table (time, RMSD, contact, salt bridge).

Usage: 02_analyze_md.py <replica e.g. prod_rep1> <n_tox_res: 32|33> <vsd_type: Kv42|Kv43>
"""
import sys, numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align

REP = sys.argv[1] if len(sys.argv) > 1 else "prod_rep1"
N_TOX_RES = int(sys.argv[2]) if len(sys.argv) > 2 else 32
VSD_TYPE = sys.argv[3] if len(sys.argv) > 3 else "Kv42"
TOP = "npt2.gro"
XTC = f"{REP}.xtc"

# local (1-indexed within VSD) S1-S4 boundaries - same table as 02_orient_place.py
LOCAL_TM = {
    "Kv42": [(1, 22), (43, 65), (73, 97), (104, 123)],
    "Kv43": [(1, 22), (42, 64), (72, 95), (102, 124)],
}[VSD_TYPE]
PADDLE_START = {"Kv42": 100, "Kv43": 98}[VSD_TYPE]   # S4 + a bit of S3b

u = mda.Universe(TOP, XTC)
prot = u.select_atoms("protein")
residues = prot.residues
tox_res = residues[:N_TOX_RES]
vsd_res = residues[N_TOX_RES:]
tox = tox_res.atoms
vsd = vsd_res.atoms
tox_ca = tox.select_atoms("name CA")
vsd_ca = vsd.select_atoms("name CA")

def local_atoms(lo, hi):
    return vsd_res[lo-1:hi].atoms
tm_ca = sum((local_atoms(a, b) for a, b in LOCAL_TM), vsd_res[0:0].atoms).select_atoms("name CA")
paddle = local_atoms(PADDLE_START, LOCAL_TM[-1][1])
print(f"[{REP}] toxin residues: {len(tox_res)} | VSD residues: {len(vsd_res)} | "
      f"TM-core CA: {len(tm_ca)} | paddle atoms: {len(paddle)}")
assert len(tox_res) == N_TOX_RES and len(tm_ca) > 20, "residue split looks wrong - check topology"

print(f"=== {REP}: {len(u.trajectory)} frames, {u.trajectory.totaltime/1000:.0f} ns ===")

acidic = vsd.select_atoms("(resname ASP GLU) and (name OD1 or name OD2 or name OE1 or name OE2)")
r3 = tox_res[2].atoms.select_atoms("name NH1 or name NH2 or name NE")  # residue index 2 = resid 3 (R3)
assert tox_res[2].resname == "ARG", f"expected R3 (Arg), got {tox_res[2].resname}"

u.trajectory[0]
ref_tox = tox_ca.positions.copy()
ref_tm_com = tm_ca.center_of_mass()
ref_tm_pos = tm_ca.positions - ref_tm_com

times, lrmsd, contact_frac, sb_frac = [], [], [], []
for ts in u.trajectory:
    times.append(ts.time / 1000)
    tm_com = tm_ca.center_of_mass()
    mobile = tm_ca.positions - tm_com
    R, _ = align.rotation_matrix(mobile, ref_tm_pos)
    mob_tox = (tox_ca.positions - tm_com) @ R.T + ref_tm_com
    lrmsd.append(float(np.sqrt(((mob_tox - ref_tox) ** 2).sum(1).mean())))
    dmat = np.linalg.norm(tox.positions[:, None, :] - paddle.positions[None, :, :], axis=2)
    contact_frac.append(float((dmat < 4.5).any()))
    if len(r3) and len(acidic):
        dsb = np.linalg.norm(r3.positions[:, None, :] - acidic.positions[None, :, :], axis=2)
        sb_frac.append(float((dsb < 4.0).any()))
    else:
        sb_frac.append(0.0)

times = np.array(times); lrmsd = np.array(lrmsd)
print(f"\nligand RMSD (toxin, superposed on VSD TM-core):")
print(f"  mean {lrmsd.mean():.2f} A | first 20ns {lrmsd[times<20].mean():.2f} A | "
      f"last 50ns {lrmsd[times>times.max()-50].mean():.2f} A | max {lrmsd.max():.2f} A")
print(f"paddle contact occupancy: {100*np.mean(contact_frac):.0f}% of frames")
print(f"R3-acidic salt bridge occupancy: {100*np.mean(sb_frac):.0f}% of frames")
half = len(times)//2
print(f"stability: first-half RMSD std {lrmsd[:half].std():.2f} A | second-half std {lrmsd[half:].std():.2f} A")

out = f"{REP}_rmsd.dat"
np.savetxt(out, np.column_stack([times, lrmsd, contact_frac, sb_frac]),
           header="time_ns ligand_RMSD_A paddle_contact salt_bridge_R3")
print("wrote", out)
