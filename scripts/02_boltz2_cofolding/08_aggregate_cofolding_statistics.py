#!/usr/bin/env python3
"""
Boltz-2 co-folding statistics (run under boltz_env):
  (1) SCRAMBLE NULL  — convergence & ipTM for 6 independent scrambled CsTx2 on the
      Kv4.2 VSD (main seed-42 scramble + 5 new shuffles). Gives a noise baseline.
  (2) MULTI-SEED     — the 4 real complexes across seeds 15/13/1/999 -> mean ± sd
      (between-run error bars on convergence & ipTM).
  (3) VERDICT        — is each real complex's convergence an outlier vs the null?

Reproduces the "MULTI-SEED REPRODUCIBILITY AND SCRAMBLE-NULL STATISTICS" section
of results/boltz2_cofolding_summary.txt.
"""
import os, glob, re, json
import numpy as np
import evaluate_boltz as E

REAL = ["CsTx2_Kv42", "CsTx2_Kv43", "CsTx3_Kv42", "CsTx3_Kv43"]
SEED_DIRS = {15: "out", 13: "out_s13", 1: "out_s1", 999: "out_s999"}

# Root containing the seed-named directories (out/, out_s13/, out_extra/, ...) from the
# primary Boltz-2 prediction data (not shipped in this repository - see the top-level
# README for how to obtain it). Set CSIMONI_TIER1_ROOT to wherever you've placed it,
# e.g. generated/predictions_raw; default "." assumes it's the current directory.
DATA_ROOT = os.environ.get("CSIMONI_TIER1_ROOT", ".")


def resolve_predictions(seed_dirname, extra=False):
    """Locate <seed_dirname>/predictions, supporting both directory shapes:
      - the original working tree: <seed_dirname>/boltz_results_inputs[_extra]/predictions
      - this repo's flattened layout: <seed_dirname>/predictions
    Tries the original shape first (matches every environment this ever ran in before
    the repo reorganization), falls back to the flat one, else returns the original
    shape unchanged so error messages point at the expected/documented path."""
    base = os.path.join(DATA_ROOT, seed_dirname)
    wrapper = "boltz_results_inputs_extra" if extra else "boltz_results_inputs"
    wrapped = os.path.join(base, wrapper, "predictions")
    flat = os.path.join(base, "predictions")
    if os.path.isdir(wrapped):
        return wrapped
    if os.path.isdir(flat):
        return flat
    return wrapped


def metrics(pred_root, name):
    d = os.path.join(pred_root, name)
    pdbs = sorted(glob.glob(os.path.join(d, f"{name}_model_*.pdb")),
                  key=lambda x: int(re.search(r"model_(\d+)", x).group(1)))
    if not pdbs:
        return None
    iptms = []
    for pp in pdbs:
        idx = re.search(r"model_(\d+)", pp).group(1)
        cj = os.path.join(d, f"confidence_{name}_model_{idx}.json")
        c = json.load(open(cj)) if os.path.exists(cj) else {}
        iptms.append(c.get("iptm", c.get("protein_iptm", 0)))
    iptms = np.array(iptms)
    cm = [m for m in (E.parse_pdb(pp) for pp in pdbs) if "B" in m and m["B"]["ca"]]
    conv, _ = E.convergence(cm)
    return dict(best=float(iptms.max()), mean=float(iptms.mean()), conv=conv * 100)


print("=" * 80)
print("BOLTZ-2 CO-FOLDING STATISTICS — scramble null + multi-seed error bars")
print("=" * 80)

# (1) scramble null
null = []
print("\n[1] SCRAMBLE NULL (CsTx2 shuffled, knot fixed, vs Kv4.2 VSD):")
print(f"  {'scramble':22}{'conv%':>7}{'best ipTM':>11}{'mean ipTM':>11}")
sources = [(resolve_predictions("out"), "scramCsTx2_Kv42", "seed42")]
for i in range(1, 6):
    sources.append((resolve_predictions("out_extra", extra=True),
                    f"scramCsTx2_Kv42_s{i}", f"shuf{i}"))
for root, name, tag in sources:
    m = metrics(root, name)
    if m:
        null.append(m); print(f"  {tag:22}{m['conv']:7.0f}{m['best']:11.2f}{m['mean']:11.2f}")
if null:
    cv = np.array([m["conv"] for m in null]); bi = np.array([m["best"] for m in null])
    print(f"  {'NULL mean±sd':22}{cv.mean():6.0f}±{cv.std():.0f}{bi.mean():10.2f}±{bi.std():.2f}")
    print(f"  {'NULL max':22}{cv.max():7.0f}{bi.max():11.2f}")
    null_conv_mean, null_conv_sd = cv.mean(), cv.std()
else:
    null_conv_mean = null_conv_sd = None

# (2) multi-seed real complexes
print("\n[2] REAL COMPLEXES across seeds 15/13/1/999:")
print(f"  {'complex':14}{'conv% (mean±sd)':>18}{'best ipTM (mean±sd)':>22}  seeds")
real_summary = {}
for name in REAL:
    convs, bests, used = [], [], []
    for seed, od in SEED_DIRS.items():
        m = metrics(resolve_predictions(od), name)
        if m:
            convs.append(m["conv"]); bests.append(m["best"]); used.append(seed)
    if convs:
        cv = np.array(convs); bi = np.array(bests)
        real_summary[name] = (cv.mean(), cv.std())
        print(f"  {name:14}{cv.mean():9.0f} ± {cv.std():<5.0f}{bi.mean():14.2f} ± {bi.std():<5.2f}  {used}")

# (3) verdict
if null_conv_mean is not None:
    thr = null_conv_mean + 2 * null_conv_sd
    print(f"\n[3] VERDICT (real convergence vs scramble null):")
    print(f"  null convergence = {null_conv_mean:.0f} ± {null_conv_sd:.0f}%  "
          f"-> outlier threshold (mean+2sd) = {thr:.0f}%")
    for name, (mu, sd) in real_summary.items():
        verdict = "OUTLIER (specific)" if mu > thr else "within noise"
        print(f"  {name:14} {mu:.0f}% -> {verdict}")
    print("\n  A real complex above the null threshold means its pose convergence is")
    print("  specific to the cognate sequence, not achievable by a random ICK peptide.")
