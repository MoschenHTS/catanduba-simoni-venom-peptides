#!/usr/bin/env python3
"""
Two checks per replica, protein-only, full 200ns trajectory:
  (1) Interface hydrogen bonds (toxin<->VSD), beyond the single R3 salt bridge
      already tracked in 02_analyze_md.py - occupancy per donor/acceptor pair.
  (2) Toxin secondary-structure stability (DSSP): does the ICK beta-hairpin
      fold (claimed in the manuscript: two antiparallel beta-strands) persist
      throughout MD, or fray/unfold?

Same index-based chain split as 02_analyze_md.py/04_analyze_membrane.py (resid resets
per chain in this topology - classify H-bonds by ATOM INDEX, not resid/segid).

Usage: 06_analyze_hbonds_ss.py <system_dir> <replica e.g. prod_rep1> <n_tox_res>
"""
import sys, os, numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.dssp import DSSP

SYS, REP, N_TOX_RES = sys.argv[1], sys.argv[2], int(sys.argv[3])
os.chdir(SYS)

STEP = 10  # subsample: 400 frames across the full 200ns (0.76s/frame for HBA -> ~5min/replica)

u = mda.Universe("npt2.gro", f"{REP}.xtc")
prot = u.select_atoms("protein")
tox_res = prot.residues[:N_TOX_RES]
n_tox_atoms = len(tox_res.atoms)
tox_atom_ix = set(tox_res.atoms.ix)   # global atom indices belonging to the toxin

# ---------------- (1) interface H-bonds ----------------
# No charge info in this .gro-based Universe (same reason we can't read the .tpr
# directly - see 02_analyze_md.py notes), so guess_hydrogens/guess_acceptors (which
# need atom charges) fail. Use bond-based pairing instead: guess covalent bonds,
# then give explicit N/O/S-based donor/acceptor/hydrogen selections.
u.atoms.guess_bonds()
hba = HydrogenBondAnalysis(u, d_a_cutoff=3.5, d_h_a_angle_cutoff=150, update_selections=False,
                            donors_sel="protein and (name N* or name O* or name S*)",
                            hydrogens_sel="protein and name H*",
                            acceptors_sel="protein and (name O* or name N* or name S*)")
hba.run(step=STEP, verbose=False)
hb = hba.results.hbonds   # columns: frame, donor_ix, hydrogen_ix, acceptor_ix, dist, angle
n_frames = len(range(0, len(u.trajectory), STEP))

pair_count = {}
for row in hb:
    d_ix, a_ix = int(row[1]), int(row[3])
    d_tox, a_tox = d_ix in tox_atom_ix, a_ix in tox_atom_ix
    if d_tox == a_tox:
        continue  # both toxin or both VSD -> intramolecular, skip
    d_atom, a_atom = u.atoms[d_ix], u.atoms[a_ix]
    key = (f"{'TOX' if d_tox else 'VSD'}:{d_atom.resname}{d_atom.resid}:{d_atom.name}",
           f"{'TOX' if a_tox else 'VSD'}:{a_atom.resname}{a_atom.resid}:{a_atom.name}")
    pair_count[key] = pair_count.get(key, 0) + 1

print(f"=== [{SYS}/{REP}] interface H-bonds (>=10% occupancy of {n_frames} frames) ===")
persistent = sorted(pair_count.items(), key=lambda x: -x[1])
for (d, a), c in persistent:
    occ = 100 * c / n_frames
    if occ >= 10:
        print(f"  {d} -- {a}   {occ:.0f}%")
if not any(100 * c / n_frames >= 10 for c in pair_count.values()):
    print("  (none >= 10% occupancy)")

# Machine-readable dump of EVERY pair (no >=10% cutoff) so downstream figures/tables
# never have to parse stdout. The occupancies used to live only in the batch-run logs,
# which meant a replica skipped by an idempotent re-run left no recoverable record
# (this bit us for cstx2_kv42/prod_rep1) - and absence there is indistinguishable from
# a real 0%. Full precision here; round at presentation time, not at storage time.
hb_out = f"{REP}_hbonds.csv"
with open(hb_out, "w") as fh:
    fh.write("donor,acceptor,n_frames_present,n_frames_total,occupancy_pct\n")
    for (d, a), c in persistent:
        fh.write(f"{d},{a},{c},{n_frames},{100 * c / n_frames:.6f}\n")
print("wrote", hb_out)

# ---------------- (2) toxin secondary structure stability ----------------
dssp = DSSP(prot, guess_hydrogens=True).run(step=STEP, verbose=False)
ss = dssp.results.dssp  # shape (n_frames, n_residues), single-letter SS codes
tox_ss = ss[:, :N_TOX_RES]

# ICK beta-hairpin: fraction of frames each toxin residue is assigned beta-strand ('E')
is_beta = (tox_ss == "E")
beta_frac_per_res = is_beta.mean(axis=0)
overall_beta_frac = is_beta.mean()
print(f"\n=== [{SYS}/{REP}] toxin secondary structure (DSSP, {N_TOX_RES} residues) ===")
print(f"overall beta-strand occupancy (any toxin residue, avg over trajectory): {100*overall_beta_frac:.0f}%")
print("per-residue beta-strand occupancy (residues with >=30%):")
for i, f in enumerate(beta_frac_per_res):
    if f >= 0.30:
        print(f"  res{i+1} ({tox_res[i].resname}): {100*f:.0f}%")

out = f"{REP}_ss.dat"
np.savetxt(out, is_beta.astype(int), fmt="%d",
           header=f"per-frame beta-strand assignment (1/0), columns=toxin residues 1-{N_TOX_RES}")
print("wrote", out)
