#!/usr/bin/env python3
"""
Boltz-2 co-folding extra inputs: a SCRAMBLE NULL DISTRIBUTION for the negative control.

The main run had one scrambled CsTx2 (shuffle seed 42). A single scramble is n=1.
Here we generate 5 independent scrambles (shuffle seeds 1-5), each keeping the six
ICK cysteines fixed (knot preserved) and permuting the rest, paired with the Kv4.2
VSD. Running all five gives a NULL distribution of pose-convergence and ipTM for
"a non-cognate ICK peptide on the Kv4.2 paddle" — so we can say whether the real
CsTx2/CsTx3 convergence is an outlier vs noise, with a proper spread (not n=1).
"""
import os
from make_inputs import PEPTIDES, VSD, CSTX_SS, CSTX2_CYS0, scramble, complex_yaml

OUT = os.path.join(os.path.dirname(__file__), "inputs_extra")
os.makedirs(OUT, exist_ok=True)


def write(name, text):
    with open(os.path.join(OUT, name + ".yaml"), "w") as f:
        f.write(text)
    print("wrote", name + ".yaml")


for s in range(1, 6):
    scr = scramble(PEPTIDES["CsTx2"], CSTX2_CYS0, seed=s)
    write(f"scramCsTx2_Kv42_s{s}", complex_yaml(scr, VSD["Kv42"], CSTX_SS))
    print(f"   scramble seed {s}: {scr}")

print("\n5 independent scrambles (Kv4.2 VSD) in", OUT)
