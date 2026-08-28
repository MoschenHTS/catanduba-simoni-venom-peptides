#!/usr/bin/env python3
"""
Evaluate Boltz-2 predictions — ENSEMBLE-aware (run under the boltz_env python:
needs numpy). Reads out/boltz_results_inputs/predictions/<name>/.

Why ensemble-aware: with 50 diffusion samples the single best-by-confidence model
is a cherry-pick. A binding pose is only trustworthy if the ensemble CONVERGES on
it. So for every complex we report, across all samples:
  * ipTM distribution (best / mean / std)            <- is the score robust or lucky?
  * pose convergence: fraction of samples whose toxin pose (chain A) reproduces the
    top pose after receptor (chain B) superposition   <- do samples agree?
  * interface PAE of the best model (inter-chain)     <- per-pair confidence
  * ICK disulfide check (per-peptide expected pattern)
  * interface residues (peptide + receptor)
DockQ on the ProTx2 benchmark and TM-align of monomers vs 1LA4 are separate
scripts (04_validate_benchmark.py, 05_validate_fold.py).
"""
import json, glob, os, math, re
import numpy as np

PRED = os.environ.get("BOLTZ_PRED_ROOT", "out/boltz_results_inputs/predictions")
ICK = {(2, 16), (9, 21), (15, 25)}            # CsTx / ProTx2 core
ICK_SGTX1 = {(2, 16), (9, 21), (15, 28)}      # SGTx1 (3rd Cys at 28)
POSE_RMSD_CUT = 4.0                           # A, "same binding pose" threshold


def expected_ss(name):
    return ICK_SGTX1 if name.startswith("SGTx1") else ICK


def parse_pdb(path):
    """chain -> {'ca': {resi: xyz}, 'sg': {resi: xyz}}; preserves residue order."""
    ch = {}
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip(); c = line[21]; resi = int(line[22:26])
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        d = ch.setdefault(c, {"ca": {}, "sg": {}, "order": []})
        if atom == "CA":
            d["ca"][resi] = xyz
            if resi not in d["order"]:
                d["order"].append(resi)
        elif atom == "SG":
            d["sg"][resi] = xyz
    return ch


def ca_array(chain):
    return np.array([chain["ca"][r] for r in chain["order"] if r in chain["ca"]])


def superpose(mob, ref):
    """R,t so that R@mob+t best fits ref (Kabsch). mob,ref: (N,3) anchors."""
    mc, rc = mob.mean(0), ref.mean(0)
    H = (mob - mc).T @ (ref - rc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, rc - R @ mc


def ligand_rmsd(mi, mj):
    """toxin(chainA) Cα RMSD between two models after superposing on receptor(chainB)."""
    Ri, ti = superpose(ca_array(mi["B"]), ca_array(mj["B"]))
    a = (Ri @ ca_array(mi["A"]).T).T + ti
    b = ca_array(mj["A"])
    return float(np.sqrt(((a - b) ** 2).sum(1).mean()))


def convergence(models):
    """Greedy: pose with most neighbours within POSE_RMSD_CUT. Returns (frac, medRMSD)."""
    n = len(models)
    if n < 2:
        return 1.0, 0.0
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = ligand_rmsd(models[i], models[j])
    neigh = (D < POSE_RMSD_CUT).sum(1)
    hub = int(neigh.argmax())
    members = np.where(D[hub] < POSE_RMSD_CUT)[0]
    return len(members) / n, float(np.median(D[hub][members]))


def interface_pae(name, best_idx, lenA):
    """Mean inter-chain PAE (A<->B) from the pae npz of the best model."""
    p = os.path.join(PRED, name, f"pae_{name}_model_{best_idx}.npz")
    if not os.path.exists(p):
        return None
    pae = np.load(p)["pae"]
    if pae.shape[0] <= lenA:
        return None
    inter = np.concatenate([pae[:lenA, lenA:].ravel(), pae[lenA:, :lenA].ravel()])
    return float(inter.mean())


def contacts(model, ca_cut=8.0):
    """Interface vicinity by Cα-Cα distance (only CA is parsed here; the precise
    heavy-atom interface is computed in 04_validate_benchmark.py). 8 A ~ side-chain contact."""
    A = model.get("A"); B = model.get("B")
    if not A or not B:
        return None, None
    Aatoms = [(r, np.array(A["ca"][r])) for r in A["ca"]]
    Batoms = [(r, np.array(B["ca"][r])) for r in B["ca"]]
    cA, cB = set(), set()
    for ra, xa in Aatoms:
        for rb, xb in Batoms:
            if np.linalg.norm(xa - xb) < ca_cut:
                cA.add(ra); cB.add(rb)
    return sorted(cA), sorted(cB)


def analyse(name):
    d = os.path.join(PRED, name)
    pdbs = sorted(glob.glob(os.path.join(d, f"{name}_model_*.pdb")),
                  key=lambda p: int(re.search(r"model_(\d+)", p).group(1)))
    if not pdbs:
        return None
    confs = []
    for p in pdbs:
        idx = int(re.search(r"model_(\d+)", p).group(1))
        cj = os.path.join(d, f"confidence_{name}_model_{idx}.json")
        c = json.load(open(cj)) if os.path.exists(cj) else {}
        confs.append((idx, p, c))
    # best by confidence_score
    bi = max(range(len(confs)), key=lambda k: confs[k][2].get("confidence_score", 0))
    best_idx, best_pdb, best_c = confs[bi]
    models = [parse_pdb(p) for _, p, _ in confs]
    bm = models[bi]
    is_cplx = "B" in bm and bm["B"]["ca"]

    iptms = np.array([c.get("iptm", c.get("protein_iptm", 0)) for _, _, c in confs])
    out = {"name": name, "is_cplx": bool(is_cplx),
           "iptm_best": float(iptms.max()), "iptm_mean": float(iptms.mean()),
           "iptm_std": float(iptms.std()),
           "ptm": best_c.get("ptm", 0), "plddt": best_c.get("complex_plddt", 0),
           "conf": best_c.get("confidence_score", 0), "n": len(confs)}
    # SS check
    ick = expected_ss(name)
    sg = bm.get("A", {}).get("sg", {})
    formed = {(a, b) for a, b in ick
              if a in sg and b in sg and math.dist(sg[a], sg[b]) < 2.5}
    out["ss"] = "3/3" if formed == ick else f"{len(formed)}/3"
    if is_cplx:
        cplx_models = [m for m in models if "B" in m and m["B"]["ca"]]
        out["conv_frac"], out["conv_rmsd"] = convergence(cplx_models)
        out["ipae"] = interface_pae(name, best_idx, len(bm["A"]["ca"]))
        cA, cB = contacts(bm)
        out["pepC"], out["recC"] = cA, cB
    return out


print("=" * 92)
print("BOLTZ-2 ENSEMBLE EVALUATION — C. simoni (corrected VSDs, 50 samples, controls)")
print("=" * 92)
hdr = f"{'input':17}{'iptm(best/mean±sd)':>22} {'conv':>6} {'ipae':>6} {'pLDDT':>6} {'SS':>5}"
print(hdr)
print("-" * 92)

groups = {"MONOMER": [], "REAL COMPLEX": [], "POSITIVE": [], "NEGATIVE": []}
def group_of(n):
    if n in ("SGTx1_Kv42", "ProTx2_Nav17VSD2"): return "POSITIVE"
    if n in ("CsTx2_ubiquitin", "scramCsTx2_Kv42", "CsTx2_CiVSP") or n.startswith("scramCsTx2_Kv42_s"):
        return "NEGATIVE"
    if "_Kv4" in n: return "REAL COMPLEX"
    return "MONOMER"

rows = []
for name in sorted(os.listdir(PRED)) if os.path.isdir(PRED) else []:
    r = analyse(name)
    if r:
        rows.append(r); groups[group_of(name)].append(r)

for g in ("MONOMER", "REAL COMPLEX", "POSITIVE", "NEGATIVE"):
    if not groups[g]:
        continue
    print(f"\n## {g}")
    for r in groups[g]:
        if r["is_cplx"]:
            ipae = f"{r['ipae']:.1f}" if r.get("ipae") is not None else "  - "
            conv = f"{r['conv_frac']*100:.0f}%"
            it = f"{r['iptm_best']:.2f}/{r['iptm_mean']:.2f}±{r['iptm_std']:.2f}"
        else:
            ipae, conv, it = "  - ", "  - ", "      -        "
        print(f"{r['name']:17}{it:>22} {conv:>6} {ipae:>6} {r['plddt']:6.2f} {r['ss']:>5}")
        if r["is_cplx"] and r.get("recC"):
            print(f"{'':17}  pep contacts: {r['pepC']}")
            print(f"{'':17}  rec contacts: {r['recC']}  (converged pose RMSD {r['conv_rmsd']:.1f} A)")

print("\n" + "=" * 92)
print("READING IT:")
print("  iptm  : best is a cherry-pick; trust it only if conv is HIGH and std is LOW.")
print("  conv  : %% of 50 samples reproducing the top toxin pose (<%.0f A). LOW = no real" % POSE_RMSD_CUT)
print("          binding mode, just diffusion noise (this is the real signal, not best-iptm).")
print("  ipae  : mean inter-chain PAE (A). <10 confident interface, >20 essentially random.")
print("  SS    : ICK disulfides formed at <2.5 A (constrained runs ~tautological; the de-novo")
print("          test is CsTx*_free, where 3/3 means the knot forms WITHOUT being imposed).")
print("CONTROLS: ProTx2_Nav17VSD2 must score HIGH and recover the 6N4I interface (run")
print("  04_validate_benchmark.py). Negatives (ubiquitin, scram) should score LOW/non-converged.")
print("  Only trust a CsTx-Kv4 pose that beats both negatives AND converges.")
