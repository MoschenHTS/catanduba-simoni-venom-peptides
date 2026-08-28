#!/usr/bin/env python3
"""
Statistical testing for the RMSD-stability and MM-GBSA comparisons.
UNIT OF ANALYSIS = replica mean (n=8 per complex, extended from an n=3 pilot
for statistical power — see 03_power_analysis.py), NOT individual frames —
frames within a 200ns trajectory are highly autocorrelated (not independent
samples), so pooling all frames into one n=200*8=1600 test would be
pseudo-replication, a common and serious error in MD statistics.
Effect sizes (Cohen's d) are still reported alongside p-values since they
remain more informative at this sample size, but several comparisons now
reach p<0.05 that did not at n=3.
"""
import numpy as np
from scipy import stats
from itertools import combinations
from load_tables import rmsd_last50, mmgbsa_total

# Read from the derived tables (01_collect_data.py) rather than hardcoded literals, so
# these can never drift out of sync with the primary analysis outputs.
RMSD = rmsd_last50()      # mean last-50ns RMSD (A), 8 replicas each
MMGBSA = mmgbsa_total()   # dTOTAL (kcal/mol), 8 replicas each


def cohens_d(a, b):
    a, b = np.array(a), np.array(b)
    n1, n2 = len(a), len(b)
    pooled_sd = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    return (a.mean() - b.mean()) / pooled_sd if pooled_sd > 0 else np.nan


def analyse(name, data):
    print(f"\n{'='*76}\n{name}  (unit of analysis: replica mean, n=8/complex)\n{'='*76}")

    groups = list(data.values())
    labels = list(data.keys())
    f, p_anova = stats.f_oneway(*groups)
    print(f"One-way ANOVA (omnibus, all 4 groups): F={f:.2f}, p={p_anova:.3f}")
    h, p_kw = stats.kruskal(*groups)
    print(f"Kruskal-Wallis (non-parametric omnibus): H={h:.2f}, p={p_kw:.3f}")

    d_hdr = "Cohen's d"
    print(f"\n{'comparison':30}{'mean diff':>11}{'Welch t':>9}{'p-value':>9}{d_hdr:>11}  interpretation")
    for l1, l2 in combinations(labels, 2):
        a, b = data[l1], data[l2]
        t, p = stats.ttest_ind(a, b, equal_var=False)  # Welch's t-test
        d = cohens_d(a, b)
        md = np.mean(a) - np.mean(b)
        mag = "negligible" if abs(d) < 0.2 else "small" if abs(d) < 0.5 else \
              "medium" if abs(d) < 0.8 else "large" if abs(d) < 1.2 else "very large"
        sig = "*" if p < 0.05 else " "
        print(f"{l1+' vs '+l2:30}{md:+11.2f}{t:9.2f}{p:9.3f}{sig}{d:+11.2f}  {mag} effect")

    print("\none-sample t-test vs 0 (is each complex's mean significantly non-zero?):")
    for l, v in data.items():
        t, p = stats.ttest_1samp(v, 0)
        print(f"  {l:16} mean={np.mean(v):+7.2f}  t={t:+.2f}  p={p:.3f}{'  *' if p<0.05 else ''}")


analyse("RMSD (Angstrom, lower = more converged pose)", RMSD)
analyse("MM-GBSA dG (kcal/mol, lower/negative = more favorable)", MMGBSA)

print(f"\n{'='*76}")
print("CAVEAT: even at n=8/group, power is moderate, not high, for medium effects")
print("(see 03_power_analysis.py). Absence of significance for a given pair does NOT mean")
print("absence of a real difference - report effect sizes alongside p-values, consistent")
print("with current best practice for small-n MD replica comparisons.")
