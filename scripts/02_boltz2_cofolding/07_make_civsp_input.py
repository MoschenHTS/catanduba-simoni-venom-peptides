#!/usr/bin/env python3
"""
Boltz-2 co-folding extra negative control: a REAL, non-target voltage-sensor domain, to test
whether CsTx2 binds any S1-S4-shaped four-helix bundle or specifically the Kv4
paddle. The existing negative controls (CsTx2_ubiquitin, scramCsTx2_Kv42*) test
"wrong receptor shape" and "wrong ligand sequence" respectively but never "right
shape, wrong (non-paddle-toxin) target" - this fills that gap.

Receptor: the isolated voltage-sensing domain of Ci-VSP (Ciona intestinalis
voltage-sensor-containing phosphatase) - a real, structurally solved (PDB 4G80,
Li et al. 2014 NSMB) S1-S4 voltage sensor with NO known gating-modifier-toxin
pharmacology (it's not an ion channel; toxins that target Kv/Nav paddles have
never been reported to act on it). UniProt F6XHE4 residues 100-238 (139 aa),
boundary verified against the 4G80 crystallized construct - not typed from
memory, pulled and checked programmatically (see conversation/session notes):
  * full UniProt sequence length 570 aa confirmed
  * 4G80's actual crystallized VSD sequence (minus the His-tag/TEV artifact)
    confirmed as an exact substring of this 100-238 slice
  * residues 217/223/226 confirmed Arg, matching the literature S4 gating-charge
    positions R1-R3 (Murata & Okamura 2005)
"""
import os
from make_inputs import PEPTIDES, CSTX_SS, complex_yaml

OUT = os.path.join(os.path.dirname(__file__), "inputs_civsp")
os.makedirs(OUT, exist_ok=True)

# Ci-VSP isolated VSD, UniProt F6XHE4 res 100-238 (see docstring for verification)
CIVSP_VSD = ("QFRVRAVIDHLGMRVFGVFLIFLDIILMIIDLSLPGKSESSQSFYDGMALALSCYFMLDLGLRIFAYGPKNFFTNPWEVAD"
             "GLIIVVTFVVTIFYTVLDEYVQETGADGLGRLVVLARLLRVVRLARIFYSHQQMKASS")
assert len(CIVSP_VSD) == 139, f"unexpected length {len(CIVSP_VSD)}"


def write(name, text):
    with open(os.path.join(OUT, name + ".yaml"), "w") as f:
        f.write(text)
    print("wrote", name + ".yaml")


write("CsTx2_CiVSP", complex_yaml(PEPTIDES["CsTx2"], CIVSP_VSD, CSTX_SS))
print(f"\nCsTx2_CiVSP: peptide {len(PEPTIDES['CsTx2'])} aa, receptor {len(CIVSP_VSD)} aa")
print("Non-target VSD negative control in", OUT)
