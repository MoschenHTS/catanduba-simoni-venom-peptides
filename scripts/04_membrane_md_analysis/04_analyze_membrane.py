#!/usr/bin/env python3
"""
Membrane partitioning: does the toxin's hydrophobic face (W4,M5,F6,W27 - the
Boltz-2 co-folding pharmacophore) actually approach/insert into the lipid acyl tails, as
the manuscript's Discussion claims ("amphipathic surface... facilitates
partitioning into the membrane")? Tests this directly against a cationic
control (R3) which should NOT partition into the hydrophobic core.

Uses the FULL (lipid-containing) trajectory - npt2.gro topology + prod_repN.xtc
(NOT the protein-only mmpbsa/ trajectories). Same residue-index gotcha as
02_analyze_md.py: resid resets per chain, so select by index not resid.

Usage: 04_analyze_membrane.py <system_dir> <replica e.g. prod_rep1> <n_tox_res>
"""
import sys, numpy as np
import MDAnalysis as mda

SYS, REP, N_TOX_RES = sys.argv[1], sys.argv[2], int(sys.argv[3])
import os
os.chdir(SYS)

u = mda.Universe("npt2.gro", f"{REP}.xtc")
prot = u.select_atoms("protein")
tox_res = prot.residues[:N_TOX_RES]

hydrophobic = (tox_res[3].atoms + tox_res[4].atoms + tox_res[5].atoms + tox_res[26].atoms)  # W4,M5,F6,W27 (0-idx 3,4,5,26)
cationic = tox_res[2].atoms  # R3 (0-idx 2)
assert tox_res[3].resname == "TRP" and tox_res[4].resname == "MET" and tox_res[5].resname == "PHE" \
    and tox_res[26].resname == "TRP", "hydrophobic-face residue check failed"
assert tox_res[2].resname == "ARG", "R3 check failed"

tails = u.select_atoms("resname POPC and (name C2* or name C3*)")
popc_p = u.select_atoms("resname POPC and name P")
print(f"[{SYS}/{REP}] hydrophobic-face atoms: {len(hydrophobic)} | R3 atoms: {len(cationic)} | "
      f"lipid tail atoms: {len(tails)} | POPC P atoms: {len(popc_p)}")

times, hydro_mindist, cat_mindist, hydro_z, cat_z, leaflet_z = [], [], [], [], [], []
for ts in u.trajectory:
    times.append(ts.time / 1000)
    tp = tails.positions
    hp = hydrophobic.positions
    cp = cationic.positions
    dh = np.sqrt(((hp[:, None, :] - tp[None, :, :]) ** 2).sum(-1)).min()
    dc = np.sqrt(((cp[:, None, :] - tp[None, :, :]) ** 2).sum(-1)).min()
    hydro_mindist.append(dh)
    cat_mindist.append(dc)
    hydro_z.append(hp[:, 2].mean())
    cat_z.append(cp[:, 2].mean())
    pz = popc_p.positions[:, 2]
    leaflet_z.append(pz[pz > pz.mean()].mean())  # upper leaflet (toxin side) mean P z

times = np.array(times); hydro_mindist = np.array(hydro_mindist); cat_mindist = np.array(cat_mindist)
hydro_z = np.array(hydro_z); cat_z = np.array(cat_z); leaflet_z = np.array(leaflet_z)

INSERT_CUT = 5.0  # A, "close to/inserted in tail region"
print(f"hydrophobic face (W4/M5/F6/W27) min-dist to lipid tails: "
      f"mean {hydro_mindist.mean():.2f} A | <{ INSERT_CUT:.0f}A in {100*np.mean(hydro_mindist<INSERT_CUT):.0f}% of frames")
print(f"R3 (cationic control) min-dist to lipid tails:            "
      f"mean {cat_mindist.mean():.2f} A | <{ INSERT_CUT:.0f}A in {100*np.mean(cat_mindist<INSERT_CUT):.0f}% of frames")
print(f"hydrophobic face z rel. to upper leaflet P: mean {(hydro_z-leaflet_z).mean():+.2f} A "
      f"(negative = below/inside headgroup plane)")
print(f"R3 z rel. to upper leaflet P:               mean {(cat_z-leaflet_z).mean():+.2f} A")

out = f"{REP}_membrane.dat"
np.savetxt(out, np.column_stack([times, hydro_mindist, cat_mindist, hydro_z - leaflet_z, cat_z - leaflet_z]),
           header="time_ns hydrophobic_mindist_A cationic_mindist_A hydrophobic_z_rel_leaflet cationic_z_rel_leaflet")
print("wrote", out)
