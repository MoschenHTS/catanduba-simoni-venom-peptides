#!/usr/bin/env python3
"""
Single source of truth for every number that feeds the membrane-MD statistics and figures.

Reads ONLY primary analysis outputs and writes tidy CSVs. Previously each of
02_stats_test.py, 03_power_analysis.py, 04_make_figures.py and 05_make_hbond_figure.py carried
its own hand-transcribed copy of these values (MM-GBSA appeared in three independent
copies, RMSD in two) - correct at the time, but nothing prevented them drifting apart
after a partial re-run, and two of the transcribed H-bond values turned out to be
display-cutoff artifacts rather than real zeros (see hbond note below).

Sources (per-replica outputs of scripts/04_membrane_md_analysis/, not shipped in
this repository - see the top-level README for how to obtain the primary MD data):
  <sys>/prod_repN_rmsd.dat                      -> md_rmsd_paddle_contact.csv
  <sys>/mmpbsa/repN_FINAL_RESULTS_MMPBSA.dat    -> md_mmgbsa_binding_energy.csv
  <sys>/prod_repN_hbonds.csv                    -> md_hbond_occupancy.csv

Usage: 01_collect_data.py [outdir]   (default: ../../results, i.e. this repo's results/)
"""
import os
import sys
import csv
import numpy as np

BASE = os.environ.get("CSIMONI_MD_BASE", os.path.dirname(os.path.abspath(__file__)))
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")

SYSTEMS = ["cstx2_kv42", "cstx3_kv42", "cstx2_kv43", "cstx3_kv43"]
LABEL = {"cstx2_kv42": "CsTx2-Kv4.2", "cstx3_kv42": "CsTx3-Kv4.2",
         "cstx2_kv43": "CsTx2-Kv4.3", "cstx3_kv43": "CsTx3-Kv4.3"}
REPS = range(1, 9)
LAST_NS = 50.0   # averaging window at the end of each 200 ns replica


def rmsd_table():
    """Mean ligand RMSD over the last LAST_NS of each replica, plus the paddle-contact
    and R3 salt-bridge occupancies from the same file."""
    rows = []
    for s in SYSTEMS:
        for r in REPS:
            d = np.loadtxt(os.path.join(BASE, s, f"prod_rep{r}_rmsd.dat"))
            t, rmsd, paddle, salt = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
            m = t >= (t.max() - LAST_NS)
            rows.append({
                "system": s, "complex": LABEL[s], "replica": r,
                "rmsd_last50ns_A": round(float(rmsd[m].mean()), 4),
                "rmsd_full_mean_A": round(float(rmsd.mean()), 4),
                "paddle_contact_pct": round(100 * float(paddle.mean()), 4),
                "r3_saltbridge_pct": round(100 * float(salt.mean()), 4),
            })
    return rows, ["system", "complex", "replica", "rmsd_last50ns_A",
                  "rmsd_full_mean_A", "paddle_contact_pct", "r3_saltbridge_pct"]


def mmgbsa_table():
    """ΔTOTAL (and its components) from each replica's gmx_MMPBSA results file.

    NOTE: enthalpy-only MM-GBSA (igb=5, no entropy term) - a RELATIVE ranking metric,
    never an absolute binding affinity. See METHODS/RESULTS_INDEX.
    """
    wanted = {"ΔVDWAALS": "dVDWAALS", "ΔEEL": "dEEL", "ΔEGB": "dEGB",
              "ΔESURF": "dESURF", "ΔGGAS": "dGGAS", "ΔGSOLV": "dGSOLV",
              "ΔTOTAL": "dTOTAL"}
    rows = []
    for s in SYSTEMS:
        for r in REPS:
            p = os.path.join(BASE, s, "mmpbsa", f"rep{r}_FINAL_RESULTS_MMPBSA.dat")
            rec = {"system": s, "complex": LABEL[s], "replica": r}
            with open(p, encoding="utf-8") as fh:
                for line in fh:
                    key = line.split()[0] if line.split() else ""
                    if key in wanted:
                        rec[wanted[key]] = float(line.split()[1])
            missing = set(wanted.values()) - set(rec)
            if missing:
                raise ValueError(f"{p}: missing {sorted(missing)}")
            rows.append(rec)
    return rows, ["system", "complex", "replica", "dVDWAALS", "dEEL", "dEGB",
                  "dESURF", "dGGAS", "dGSOLV", "dTOTAL"]


def hbond_table():
    """Occupancy of the Trp4(toxin)-Thr22(VSD) interaction per replica, split by
    acceptor atom.

    Two variants are tracked separately and MUST NOT be summed or conflated:
      backbone_O  - Trp4 amide N-H to the Thr22 backbone carbonyl (the Kv4.3-specific
                    anchor plotted in the figure)
      sidechain_OG1 - the same donor to the Thr22 side-chain hydroxyl (a different,
                    weaker interaction that also appears on some Kv4.2 replicas)

    Reading these from the per-replica CSVs rather than the batch logs fixes two real
    defects in the previously transcribed values: (1) the logs only printed pairs at
    >=10% occupancy, so genuinely small but nonzero occupancies were recorded as 0
    (cstx2_kv43 rep4 = 5%, cstx3_kv43 rep5 = 3%); (2) a replica skipped by an idempotent
    re-run left no log record at all, indistinguishable from a measured zero
    (cstx2_kv42 rep1).
    """
    rows = []
    for s in SYSTEMS:
        for r in REPS:
            p = os.path.join(BASE, s, f"prod_rep{r}_hbonds.csv")
            occ = {"backbone_O": 0.0, "sidechain_OG1": 0.0}
            with open(p, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh):
                    if row["donor"] != "TOX:TRP4:N":
                        continue
                    if row["acceptor"] == "VSD:THR22:O":
                        occ["backbone_O"] = float(row["occupancy_pct"])
                    elif row["acceptor"] == "VSD:THR22:OG1":
                        occ["sidechain_OG1"] = float(row["occupancy_pct"])
            rows.append({"system": s, "complex": LABEL[s], "replica": r,
                         "trp4_thr22_backbone_O_pct": round(occ["backbone_O"], 4),
                         "trp4_thr22_sidechain_OG1_pct": round(occ["sidechain_OG1"], 4)})
    return rows, ["system", "complex", "replica", "trp4_thr22_backbone_O_pct",
                  "trp4_thr22_sidechain_OG1_pct"]


def write(name, rows, fields):
    os.makedirs(OUTDIR, exist_ok=True)
    p = os.path.join(OUTDIR, name)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {p}  ({len(rows)} rows)")
    return rows


if __name__ == "__main__":
    rmsd = write("md_rmsd_paddle_contact.csv", *rmsd_table())
    gbsa = write("md_mmgbsa_binding_energy.csv", *mmgbsa_table())
    hb = write("md_hbond_occupancy.csv", *hbond_table())

    print("\nper-complex summary (mean +/- SD over 8 replicas):")
    print(f"{'complex':14}{'RMSD last50 (A)':>20}{'MM-GBSA dTOTAL':>20}")
    for s in SYSTEMS:
        rv = np.array([r["rmsd_last50ns_A"] for r in rmsd if r["system"] == s])
        gv = np.array([r["dTOTAL"] for r in gbsa if r["system"] == s])
        print(f"{LABEL[s]:14}{rv.mean():11.2f} +/-{rv.std(ddof=1):5.2f}"
              f"{gv.mean():12.2f} +/-{gv.std(ddof=1):6.2f}")
