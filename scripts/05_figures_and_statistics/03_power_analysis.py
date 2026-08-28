#!/usr/bin/env python3
"""
Power analysis: given the CURRENT n=8/complex data (extended from the
original n=3 pilot that motivated the extension in the first place), how many
replicas would actually be needed to reach 80% power at alpha=0.05 for each
comparison that still isn't significant? Uses observed Cohen's d from n=8 -
IMPORTANT CAVEAT: this is still a plug-in effect-size estimate, not a true
prospective power analysis; required-N numbers are planning estimates, not
guarantees. (Replicas were extended from an n=3 pilot to n=8 for statistical
power - several n=3-based low-N predictions turned out to be
inflated pilot effect sizes that shrank sharply with more data.)
"""
import numpy as np
from statsmodels.stats.power import TTestIndPower
from itertools import combinations
from load_tables import rmsd_last50, mmgbsa_total

# Read from the derived tables (01_collect_data.py), same source as 02_stats_test.py.
RMSD = rmsd_last50()
MMGBSA = mmgbsa_total()

analysis = TTestIndPower()


def cohens_d(a, b):
    a, b = np.array(a), np.array(b)
    n1, n2 = len(a), len(b)
    pooled_sd = np.sqrt(((n1-1)*a.var(ddof=1) + (n2-1)*b.var(ddof=1)) / (n1+n2-2))
    return (a.mean() - b.mean()) / pooled_sd if pooled_sd > 0 else np.nan


def report(name, data):
    print(f"\n{'='*90}\n{name}\n{'='*90}")
    print(f"{'comparison':30}{'|d|':>7}{'power@n=3':>11}{'power@n=6':>11}{'power@n=10':>12}{'n for 80% power':>18}")
    for l1, l2 in combinations(data.keys(), 2):
        d = abs(cohens_d(data[l1], data[l2]))
        if np.isnan(d) or d == 0:
            continue
        pw3 = analysis.power(effect_size=d, nobs1=3, alpha=0.05, ratio=1.0)
        pw6 = analysis.power(effect_size=d, nobs1=6, alpha=0.05, ratio=1.0)
        pw10 = analysis.power(effect_size=d, nobs1=10, alpha=0.05, ratio=1.0)
        try:
            n_req = analysis.solve_power(effect_size=d, alpha=0.05, power=0.8, ratio=1.0)
            n_req_s = f"{n_req:.0f}"
        except Exception:
            n_req_s = ">100"
        print(f"{l1+' vs '+l2:30}{d:7.2f}{pw3:11.2f}{pw6:11.2f}{pw10:12.2f}{n_req_s:>18}")


report("RMSD", RMSD)
report("MM-GBSA", MMGBSA)

print(f"\n{'='*90}")
print("CAVEAT: SD is a plug-in estimate from the current n=8/complex data (7 degrees of")
print("freedom) - these required-N figures are planning estimates, with substantial")
print("uncertainty of their own, not guarantees that the required N will actually reach")
print("significance if run.")
