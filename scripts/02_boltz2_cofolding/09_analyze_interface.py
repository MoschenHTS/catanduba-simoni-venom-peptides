#!/usr/bin/env python3
"""
Deep interface analysis across the full 50-model ensemble (run under boltz_env).

For each complex it reports, as fractions of the 50 samples (i.e. how REPRODUCIBLE
each contact is, not a single cherry-picked model):
  * per-residue contact occupancy (peptide and VSD), VSD residues flagged S4/paddle
  * persistent salt bridges (toxin acidic <-> VSD basic, and vice versa)
  * the gating-modifier PHARMACOPHORE: which cationic (R/K/H) and which hydrophobic/
    aromatic (W/F/Y/M/L/I/V) toxin residues bury at the interface, with occupancy.

The pharmacophore is computed identically for the CsTx complexes and for our two
validated reference gating modifiers — SGTx1 (closest homolog) and ProTx-II (the
solved-structure 6N4I benchmark) — so CsTx is compared to the toxins we actually
modelled, not to an external prototype.

Heavy-atom contacts at 4.5 A; salt bridge = acidic O within 4.0 A of basic N.
"""
import glob, os, re, json
import numpy as np
import warnings; warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.IUPACData import protein_letters_3to1 as _T31

ROOT = os.environ.get("BOLTZ_PRED_ROOT", "out/boltz_results_inputs/predictions")
COMPLEXES = ["CsTx2_Kv42", "CsTx2_Kv43", "CsTx3_Kv42", "CsTx3_Kv43"]
REFS = ["SGTx1_Kv42", "ProTx2_Nav17VSD2"]   # validated references (NOT hanatoxin)
CONTACT = 4.5
SB = 4.0
S4_START = {"Kv42": 100, "Kv43": 100}       # renumbered VSD res >= -> S4/paddle
ACID = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
BASE = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}
BASIC_AA = {"ARG", "LYS", "HIS"}
HYDRO_AA = {"TRP", "PHE", "TYR", "MET", "LEU", "ILE", "VAL"}
p = PDBParser(QUIET=True)


def one(resname):
    return _T31.get(resname.capitalize(), "X")


def vsd_of(name):
    return "Kv42" if "Kv42" in name else "Kv43" if "Kv43" in name else None


def models(name):
    d = os.path.join(ROOT, name)
    return sorted(glob.glob(os.path.join(d, f"{name}_model_*.pdb")),
                  key=lambda x: int(re.search(r"model_(\d+)", x).group(1)))


def chains(path):
    m = p.get_structure("m", path)[0]
    A = [r for r in m["A"] if is_aa(r)] if "A" in m else []
    B = [r for r in m["B"] if is_aa(r)] if "B" in m else []
    return A, B


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"])


def res_contacts(A, B):
    cA, cB = set(), set()
    Bh = [(r, heavy(r)) for r in B]
    for ra in A:
        ha = heavy(ra)
        for rb, hb in Bh:
            if (np.linalg.norm(ha[:, None] - hb[None], axis=2) < CONTACT).any():
                cA.add((ra.id[1], ra.resname)); cB.add((rb.id[1], rb.resname))
    return cA, cB


def salt_bridges(A, B):
    out = []
    def charged(res, table):
        return [a.coord for a in res if a.name in table.get(res.resname, [])]
    for ra in A:
        ra_ac, ra_ba = charged(ra, ACID), charged(ra, BASE)
        for rb in B:
            rb_ac, rb_ba = charged(rb, ACID), charged(rb, BASE)
            if ra_ac and rb_ba and (np.linalg.norm(np.array(ra_ac)[:, None]-np.array(rb_ba)[None], axis=2) < SB).any():
                out.append((f"{ra.resname}{ra.id[1]}", f"{rb.resname}{rb.id[1]}", "tox(-)/vsd(+)"))
            if ra_ba and rb_ac and (np.linalg.norm(np.array(ra_ba)[:, None]-np.array(rb_ac)[None], axis=2) < SB).any():
                out.append((f"{ra.resname}{ra.id[1]}", f"{rb.resname}{rb.id[1]}", "tox(+)/vsd(-)"))
    return out


def analyse(name):
    vsd = vsd_of(name)
    s4 = S4_START.get(vsd)
    ms = models(name)
    n = len(ms)
    if not n:
        print(f"  {name}: no models"); return
    pepocc, vsdocc, sbocc = {}, {}, {}
    for path in ms:
        A, B = chains(path)
        cA, cB = res_contacts(A, B)
        for r in cA: pepocc[r] = pepocc.get(r, 0) + 1
        for r in cB: vsdocc[r] = vsdocc.get(r, 0) + 1
        for sb in set(salt_bridges(A, B)):
            sbocc[sb] = sbocc.get(sb, 0) + 1
    print(f"\n## {name}  (n={n} models)")
    pep = sorted(pepocc.items(), key=lambda x: -x[1])
    print("  peptide interface (residue: % of models):")
    print("    " + "  ".join(f"{one(nm)}{ri}:{100*c//n}%" for (ri, nm), c in pep if c >= n*0.3))
    vsd_sorted = sorted(vsdocc.items(), key=lambda x: -x[1])
    print("  VSD interface (residue: %;" + (" * = S4/paddle):" if s4 else "):"))
    print("    " + "  ".join(f"{one(nm)}{ri}{'*' if s4 and ri>=s4 else ''}:{100*c//n}%"
                              for (ri, nm), c in vsd_sorted if c >= n*0.3))
    if s4:
        s4hit = [c for (ri, nm), c in vsdocc.items() if ri >= s4]
        m = max(s4hit) if s4hit else 0
        print(f"  --> S4/paddle engaged in {m}/{n} models ({100*m//n}% of ensemble)")
    if sbocc:
        print("  persistent salt bridges (>=20% of models):")
        for (pr, vr, kind), c in sorted(sbocc.items(), key=lambda x: -x[1]):
            if c >= n*0.2:
                print(f"    {pr} -- {vr}  [{kind}]  {100*c//n}%")
    # generic gating-modifier pharmacophore (cationic + hydrophobic), per toxin
    basics = sorted([((ri, nm), c) for (ri, nm), c in pepocc.items() if nm in BASIC_AA], key=lambda x: -x[1])
    hydro = sorted([((ri, nm), c) for (ri, nm), c in pepocc.items() if nm in HYDRO_AA], key=lambda x: -x[1])
    print("  pharmacophore | cationic: " +
          (" ".join(f"{one(nm)}{ri}:{100*c//n}%" for (ri, nm), c in basics if c >= n*0.3) or "none>30%"))
    print("                | hydrophobic: " +
          (" ".join(f"{one(nm)}{ri}:{100*c//n}%" for (ri, nm), c in hydro if c >= n*0.3) or "none>30%"))


print("=" * 80)
print("DEEP INTERFACE ANALYSIS — contact occupancy across the 50-model ensemble")
print("=" * 80)
print("Pharmacophore = cationic (R/K/H) + hydrophobic (W/F/Y/M/L/I/V) toxin residues")
print("that bury at the interface, computed identically for CsTx and for our two")
print("validated references: SGTx1 (homolog) and ProTx-II (6N4I benchmark).")

print("\n" + "-" * 30 + " CsTx COMPLEXES " + "-" * 30)
for c in COMPLEXES:
    analyse(c)

print("\n" + "-" * 26 + " REFERENCE GATING MODIFIERS " + "-" * 26)
for c in REFS:
    analyse(c)
