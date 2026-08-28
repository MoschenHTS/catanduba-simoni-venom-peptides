#!/usr/bin/env python3
"""
Positive-control BENCHMARK: did Boltz reproduce the SOLVED ProTx-II / Nav1.7 VSD2
interface (PDB 6N4I)?  Run under boltz_env python (numpy + biopython).

This is the only ground-truth test in the pipeline. We score the predicted
ProTx2_Nav17VSD2 models against ref/native_protx2_vsd2.pdb using the DockQ
component metrics (computed directly, no external dep):

  fnat   = fraction of native receptor-toxin residue contacts recovered  (>0.5 good)
  L-RMSD = toxin backbone RMSD after superposing on the receptor          (<5 A good)
  i-RMSD = interface backbone RMSD                                         (<4 A good)

If the best model AND a good fraction of the 50-sample ensemble recover the
interface, the pipeline can place an ICK paddle toxin -> trust the CsTx-Kv4 poses
more. If it cannot reproduce a SOLVED complex, treat all CsTx-Kv4 poses as
unreliable regardless of ipTM.
"""
import glob, os, re, json, sys
import numpy as np
import warnings; warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

NAME = "ProTx2_Nav17VSD2"
ROOT = os.environ.get("BOLTZ_PRED_ROOT", "out/boltz_results_inputs/predictions")
PRED = f"{ROOT}/{NAME}"
# Default "ref/" preserves the original working-tree layout unchanged. Set to
# <repo>/inputs/reference_structures to run against the reorganized/portable repo instead
# (see CSIMONI_TIER1_ROOT in 08_aggregate_cofolding_statistics.py for the same pattern).
REF_ROOT = os.environ.get("BOLTZ_REF_ROOT", "ref")
NATIVE = f"{REF_ROOT}/native_protx2_vsd2.pdb"
CONTACT = 5.0   # heavy-atom contact cutoff (A)

p = PDBParser(QUIET=True)


def chains_ordered(struct):
    """return {chain_id: [residue,...]} amino acids in order."""
    m = struct[0]
    return {c.id: [r for r in c if is_aa(r)] for c in m}


def heavy(res):
    return np.array([a.coord for a in res if a.element != "H"])


def ca(res):
    return res["CA"].coord if "CA" in res else res.child_list[0].coord


def native_contacts(rec, tox):
    """set of (i,j) order-indices with any heavy-atom pair < CONTACT."""
    out = set()
    th = [heavy(t) for t in tox]
    for i, r in enumerate(rec):
        rh = heavy(r)
        for j, t in enumerate(tox):
            d = np.linalg.norm(rh[:, None, :] - th[j][None, :, :], axis=2)
            if (d < CONTACT).any():
                out.add((i, j))
    return out


def superpose(mob, ref):
    mc, rc = mob.mean(0), ref.mean(0)
    H = (mob - mc).T @ (ref - rc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, rc - R @ mc


def identify(model_chains):
    """which model chain is toxin (~30 aa) vs receptor (~125 aa)."""
    sizes = {c: len(rs) for c, rs in model_chains.items()}
    tox = min(sizes, key=sizes.get)
    rec = max(sizes, key=sizes.get)
    return rec, tox


def score(model_path, nat_rec, nat_tox, nat_cset, iface_rec, iface_tox):
    mc = chains_ordered(p.get_structure("m", model_path))
    rc_id, tx_id = identify(mc)
    m_rec, m_tox = mc[rc_id], mc[tx_id]
    if len(m_rec) != len(nat_rec) or len(m_tox) != len(nat_tox):
        return None  # sequence/length mismatch
    # fnat
    mset = native_contacts(m_rec, m_tox)
    fnat = len(mset & nat_cset) / max(1, len(nat_cset))
    # L-RMSD: superpose on receptor CA, measure toxin CA
    R, t = superpose(np.array([ca(r) for r in m_rec]),
                     np.array([ca(r) for r in nat_rec]))
    mt = (R @ np.array([ca(r) for r in m_tox]).T).T + t
    nt = np.array([ca(r) for r in nat_tox])
    lrms = float(np.sqrt(((mt - nt) ** 2).sum(1).mean()))
    # i-RMSD: superpose on interface CA only
    mi = np.array([ca(m_rec[i]) for i in iface_rec] + [ca(m_tox[j]) for j in iface_tox])
    ni = np.array([ca(nat_rec[i]) for i in iface_rec] + [ca(nat_tox[j]) for j in iface_tox])
    Ri, ti = superpose(mi, ni)
    irms = float(np.sqrt((((Ri @ mi.T).T + ti - ni) ** 2).sum(1).mean()))
    return fnat, lrms, irms


def main():
    if not os.path.isdir(PRED):
        print(f"[skip] no predictions for {NAME} yet"); return
    nat = chains_ordered(p.get_structure("n", NATIVE))
    rec_id, tox_id = identify(nat)
    nat_rec, nat_tox = nat[rec_id], nat[tox_id]
    nat_cset = native_contacts(nat_rec, nat_tox)
    iface_rec = sorted({i for i, _ in nat_cset})
    iface_tox = sorted({j for _, j in nat_cset})
    print(f"Native 6N4I interface: {len(nat_cset)} residue contacts "
          f"({len(iface_rec)} VSD2 + {len(iface_tox)} ProTx2 residues)")

    pdbs = sorted(glob.glob(os.path.join(PRED, f"{NAME}_model_*.pdb")),
                  key=lambda x: int(re.search(r"model_(\d+)", x).group(1)))
    # best by confidence
    def conf(pp):
        idx = re.search(r"model_(\d+)", pp).group(1)
        cj = os.path.join(PRED, f"confidence_{NAME}_model_{idx}.json")
        return json.load(open(cj)).get("confidence_score", 0) if os.path.exists(cj) else 0
    best = max(pdbs, key=conf)

    results = []
    for pp in pdbs:
        s = score(pp, nat_rec, nat_tox, nat_cset, iface_rec, iface_tox)
        if s:
            results.append((pp, *s))
    if not results:
        print("[error] no scorable models (chain length mismatch?)"); return

    print(f"\n{'model':28}{'fnat':>7}{'L-RMSD':>9}{'i-RMSD':>9}")
    bf = next(r for r in results if r[0] == best)
    print(f"{'BEST (by confidence)':28}{bf[1]:7.2f}{bf[2]:9.2f}{bf[3]:9.2f}")
    fn = np.array([r[1] for r in results])
    lr = np.array([r[2] for r in results])
    print(f"{'ensemble mean':28}{fn.mean():7.2f}{lr.mean():9.2f}")
    print(f"{'ensemble best fnat':28}{fn.max():7.2f}")
    print(f"\n{(fn>0.5).sum()}/{len(fn)} samples recover >50% of the native interface "
          f"(fnat>0.5); {(lr<5).sum()}/{len(lr)} have L-RMSD < 5 A.")
    ok = bf[1] > 0.5 and bf[2] < 5
    print("\nVERDICT:", "PASS — pipeline reproduces the solved paddle-toxin interface."
          if ok else "WEAK — pipeline did NOT cleanly reproduce 6N4I; "
          "interpret CsTx-Kv4 poses with strong caution.")


if __name__ == "__main__":
    main()
