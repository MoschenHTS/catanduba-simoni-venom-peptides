#!/usr/bin/env python3
"""
Fold / disulfide validation for the monomers. Run under boltz_env python.

Two independent checks:

1. TM-align each peptide monomer to the SGTx1 NMR structure (PDB 1LA4). High
   TM-score on SGTx1_control validates that the pipeline reproduces the known ICK
   fold; CsTx2/CsTx3 TM-scores show how close their predicted fold is to that
   experimental ICK template.  (needs `tmtools`; skips gracefully if absent.)

2. DE-NOVO disulfide connectivity on the UNCONSTRAINED monomers (CsTx2_free,
   CsTx3_free). The constrained runs satisfy the imposed bonds almost by
   construction, so they cannot validate the connectivity. Here we measure all
   SG-SG distances and report whether the ICK pattern (2-16,9-21,15-25) or the
   rival defensin pattern (2-25,9-16,15-21) is what Boltz forms WITHOUT being
   told. ICK winning de novo is the actual evidence for the assignment.
"""
import glob, os, re, json, itertools
import numpy as np
import warnings; warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

PRED = os.environ.get("BOLTZ_PRED_ROOT", "out/boltz_results_inputs/predictions")
# Default "ref/" preserves the original working-tree layout unchanged. Set to
# <repo>/inputs/reference_structures to run against the reorganized/portable repo instead
# (see CSIMONI_TIER1_ROOT in 08_aggregate_cofolding_statistics.py for the same pattern).
REF_ROOT = os.environ.get("BOLTZ_REF_ROOT", "ref")
LA4 = f"{REF_ROOT}/1LA4.pdb"
ICK = [(2, 16), (9, 21), (15, 25)]
DEFENSIN = [(2, 25), (9, 16), (15, 21)]
p = PDBParser(QUIET=True)


def best_model(name):
    d = os.path.join(PRED, name)
    pdbs = glob.glob(os.path.join(d, f"{name}_model_*.pdb"))
    if not pdbs:
        return None
    def conf(pp):
        idx = re.search(r"model_(\d+)", pp).group(1)
        cj = os.path.join(d, f"confidence_{name}_model_{idx}.json")
        return json.load(open(cj)).get("confidence_score", 0) if os.path.exists(cj) else 0
    return max(pdbs, key=conf)


def ca_seq(struct, chain=None):
    from Bio.Data.IUPACData import protein_letters_3to1 as t31
    m = struct[0]
    ch = m[chain] if chain else next(iter(m))
    coords, seq = [], ""
    for r in ch:
        if is_aa(r) and "CA" in r:
            coords.append(r["CA"].coord)
            seq += t31.get(r.resname.capitalize(), "X")
    return np.array(coords), seq


def sg_dists(struct):
    m = struct[0]; ch = next(iter(m))
    sg = {r.id[1]: r["SG"].coord for r in ch if is_aa(r) and "SG" in r}
    def dd(pairs):
        return [(a, b, float(np.linalg.norm(sg[a] - sg[b])))
                for a, b in pairs if a in sg and b in sg]
    return dd, sg


def tmalign(name):
    try:
        from tmtools import tm_align
    except ImportError:
        return None
    bm = best_model(name)
    if not bm or not os.path.exists(LA4):
        return None
    c1, s1 = ca_seq(p.get_structure("m", bm))
    c2, s2 = ca_seq(p.get_structure("ref", LA4))  # 1LA4 model 1, chain A (SGTx1)
    r = tm_align(c1.astype(np.float64), c2.astype(np.float64), s1, s2)
    return r.tm_norm_chain2, r.rmsd  # normalized by 1LA4 length


print("=" * 70)
print("FOLD VALIDATION — TM-align vs 1LA4 (SGTx1 NMR) + de-novo disulfides")
print("=" * 70)

print("\n[1] TM-align to 1LA4 (SGTx1 ICK reference):")
for name in ["SGTx1_control", "CsTx2", "CsTx3", "CsTx2_free", "CsTx3_free"]:
    r = tmalign(name)
    if r is None:
        print(f"  {name:16} (skipped — tmtools missing or no model)")
    else:
        tm, rmsd = r
        flag = " <- fold pipeline validated" if name == "SGTx1_control" and tm > 0.5 else ""
        print(f"  {name:16} TM-score={tm:.3f}  RMSD={rmsd:.2f} A{flag}")

print("\n[2] De-novo disulfide connectivity (UNCONSTRAINED monomers):")
for name in ["CsTx2_free", "CsTx3_free"]:
    bm = best_model(name)
    if not bm:
        print(f"  {name:16} (no model)"); continue
    dd, sg = sg_dists(p.get_structure("m", bm))
    ick = dd(ICK); defn = dd(DEFENSIN)
    ick_ok = sum(1 for *_, d in ick if d < 2.5)
    def_ok = sum(1 for *_, d in defn if d < 2.5)
    verdict = ("ICK confirmed de novo" if ick_ok == 3 else
               "DEFENSIN pattern!" if def_ok == 3 else "ambiguous/partial")
    print(f"  {name:16} ICK pairs<2.5A: {ick_ok}/3  |  defensin pairs<2.5A: {def_ok}/3"
          f"  -> {verdict}")
    print(f"    {'  '.join(f'{a}-{b}:{d:.1f}' for a,b,d in ick)}")

print("\nInterpretation: SGTx1_control TM>0.5 confirms the monomer pipeline. If the")
print("UNCONSTRAINED CsTx*_free independently form the ICK pattern (3/3), the imposed")
print("2-16/9-21/15-25 connectivity is what Boltz predicts on its own — assumption validated.")
