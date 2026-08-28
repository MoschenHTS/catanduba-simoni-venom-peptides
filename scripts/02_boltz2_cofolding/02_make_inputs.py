#!/usr/bin/env python3
"""
Generate Boltz-2 input YAMLs for the C. simoni structural work.

CORRECTED + EXTENDED design (2026-06-26 lab-PC re-run). Supersedes the earlier
make_inputs that carried two fatal errors from the Mac run:

  1. RECEPTOR CONSTRUCTS WERE WRONG. The old "Kv4.2" chain B was residues
     235-405 (mid-S2 through the S6 pore) and the old "Kv4.3" chain B was
     residues 73-243 (the T1 cytoplasmic tetramerization domain — NO voltage
     sensor at all). They were not even homologous. Both are replaced here with
     the verified, homologous S1-S4 voltage-sensor domains taken from UniProt
     transmembrane annotations:
        Kv4.2 (Q9NZV8 / KCND2)  VSD = res 185-307
        Kv4.3 (Q9UK17 / KCND3)  VSD = res 183-306
     Both contain the S4 paddle gating charges (RVFRVFRIFK). Pore excluded so
     the toxin cannot dock a non-physiological monomeric pore.

  2. DISULFIDES use the CORRECT ICK connectivity C1-C4,C2-C5,C3-C6
     (CsTx Cys 2,9,15,16,21,25 -> bonds 2-16,9-21,15-25).

Validation added in this version:
  * UNCONSTRAINED monomers (CsTx2_free, CsTx3_free) — does the ICK knot form
    de novo, i.e. is the imposed connectivity what Boltz independently predicts?
  * Positive BINDING benchmark with a SOLVED structure: ProTx-II + Nav1.7 VSD2
    (PDB 6N4I). ProTx-II is an ICK tarantula gating modifier with the SAME knot
    pattern (2-16,9-21,15-25). Scored by DockQ vs the experimental interface
    (the only ground-truth test that Boltz places a paddle toxin correctly).
  * Negative controls: decoy receptor (ubiquitin) + scrambled ligand.

Run boltz with --use_msa_server (both chains) and a fixed --seed.

Inputs (inputs/) — 13 total:
  monomers (figure)   : CsTx2, CsTx3
  monomers (free)     : CsTx2_free, CsTx3_free          (fold/connectivity check)
  fold control        : SGTx1_control                   (TM-align vs PDB 1LA4)
  complexes (S1-S4)   : CsTx{2,3}_Kv4{2,3}              (4 real targets)
  pos. binding (homo) : SGTx1_Kv42                       (paddle toxin, real VSD)
  pos. benchmark      : ProTx2_Nav17VSD2                 (solved -> DockQ)
  neg. decoy receptor : CsTx2_ubiquitin
  neg. scrambled lig. : scramCsTx2_Kv42
"""
import os
import random

OUT = os.path.join(os.path.dirname(__file__), "inputs")
os.makedirs(OUT, exist_ok=True)

# --- peptides (chain A) ------------------------------------------------------
PEPTIDES = {
    "CsTx2": "GCRWMFGACKTTADCCKALACVGTCIWDGTYG",    # Gc32.3, 32 aa
    "CsTx3": "GCRWMFGACKTTADCCKALACVGTCIWDGTYGN",   # Gc33.1, 33 aa
}
# ICK disulfides for the CsTx core: Cys at 2,9,15,16,21,25 -> C1-C4,C2-C5,C3-C6
CSTX_SS = [(2, 16), (9, 21), (15, 25)]

# --- voltage-sensor domains (chain B): VERIFIED S1-S4, UniProt TM annotation --
# Kv4.2 Q9NZV8 res 185-307 ; Kv4.3 Q9UK17 res 183-306. Homologous, pore-free,
# both carry the S4 gating charges. (S1-S2 loop has Cys -> left unconstrained.)
VSD = {
    "Kv42": ("LVFYYVTGFFIAVSVIANVVETVPCGSSPGHIKELPCGERYAVAFFCLDTACVMIFTVEYLLRL"
             "AAAPSRYRFVRSVMSIIDVVAILPYYIGLVMTDNEDVSGAFVTLRVFRVFRIFKFSRHS"),
    "Kv43": ("LVFYYVTGFFIAVSVITNVVETVPCGTVPGSKELPCGERYSVAFFCLDTACVMIFTVEYLLRL"
             "FAAPSRYRFIRSVMSIIDVVAIMPYYIGLVMTNNEDVSGAFVTLRVFRVFRIFKFSRHSQG"),
}

# --- positive control: SGTx1 (closest homolog, exp. NMR structure PDB 1LA4) --
SGTX1 = "TCRYLFGGCKTTADCCKHLACRSDGKYCAWDGTF"          # P56855, Cys at 2,9,15,16,21,28
SGTX1_SS = [(2, 16), (9, 21), (15, 28)]               # ICK C1-C4, C2-C5, C3-C6

# --- positive BENCHMARK: ProTx-II + Nav1.7 VSD2 (solved, PDB 6N4I) -----------
# ProTx-II (P83476, beta/omega-theraphotoxin-Tp2a): ICK, SAME knot as CsTx.
PROTX2 = "YCQKWMWTCDSERKCCEGMVCRLWCKKKLW"             # Cys at 2,9,15,16,21,25
PROTX2_SS = [(2, 16), (9, 21), (15, 25)]
# VSD2 receptor = modeled S1-S4 of the Nav1.7 VSD2-NavAb chimera, chain C res
# 720-844 of 6N4I (the full ProTx-II epitope, res 763-823, lies inside it).
NAV17_VSD2 = ("SHMYLRITNIVESSFFTKFIIYLIVLNTLFMAMEHHPMTEEFKNVLAIGNLVFTGIFAIEIIL"
              "RIYVHRISFFKDPWSLFDSLIVTLSLVELFLADVEGLSVLRSFRLLRVFRLVTAVPQMRKIV")

# --- decoy receptor (negative): human ubiquitin (76 aa, non-channel, 0 Cys) --
UBIQUITIN = ("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQ"
             "KESTLHLVLRLRGG")

# --- scrambled ligand (negative): shuffle CsTx2 keeping the 6 knot Cys fixed --
CSTX2_CYS0 = {1, 8, 14, 15, 20, 24}  # 0-indexed Cys (1-indexed 2,9,15,16,21,25)


def scramble(seq, fixed_idx, seed=42):
    """Permute residues except those at fixed_idx (0-indexed). Composition kept;
    linear motif/interface destroyed; knot topology preserved."""
    rng = random.Random(seed)
    movable_pos = [i for i in range(len(seq)) if i not in fixed_idx]
    movable_res = [seq[i] for i in movable_pos]
    rng.shuffle(movable_res)
    out = list(seq)
    for p, r in zip(movable_pos, movable_res):
        out[p] = r
    s = "".join(out)
    assert all(s[i] == "C" for i in fixed_idx), "scramble moved a knot Cys!"
    assert sorted(s) == sorted(seq), "scramble changed composition!"
    return s


SCRAM_CSTX2 = scramble(PEPTIDES["CsTx2"], CSTX2_CYS0)


# --- YAML builders -----------------------------------------------------------
def ss_block(pairs):
    lines = ["constraints:"]
    for a, b in pairs:
        lines += ["  - bond:", f"      atom1: [A, {a}, SG]", f"      atom2: [A, {b}, SG]"]
    return "\n".join(lines)


def monomer_yaml(seq, ss=None):
    y = ("version: 1\n"
         "sequences:\n"
         "  - protein:\n"
         "      id: A\n"
         f"      sequence: {seq}\n")
    if ss:
        y += ss_block(ss) + "\n"
    return y


def complex_yaml(pep_seq, rec_seq, ss):
    return ("version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {pep_seq}\n"
            "  - protein:\n"
            "      id: B\n"
            f"      sequence: {rec_seq}\n"
            f"{ss_block(ss)}\n")


def write(name, text):
    with open(os.path.join(OUT, name + ".yaml"), "w") as f:
        f.write(text)
    print("wrote", name + ".yaml")


# --- monomers ----------------------------------------------------------------
for pname, pseq in PEPTIDES.items():
    write(pname, monomer_yaml(pseq, CSTX_SS))            # constrained (figure)
    write(pname + "_free", monomer_yaml(pseq, ss=None))  # UNCONSTRAINED (validation)
write("SGTx1_control", monomer_yaml(SGTX1, SGTX1_SS))    # fold control vs 1LA4

# --- real complexes (corrected S1-S4 VSD) -----------------------------------
for pname, pseq in PEPTIDES.items():
    for vname, vseq in VSD.items():
        write(f"{pname}_{vname}", complex_yaml(pseq, vseq, CSTX_SS))

# --- positive controls -------------------------------------------------------
write("SGTx1_Kv42", complex_yaml(SGTX1, VSD["Kv42"], SGTX1_SS))     # homology ref
write("ProTx2_Nav17VSD2", complex_yaml(PROTX2, NAV17_VSD2, PROTX2_SS))  # benchmark

# --- negative controls -------------------------------------------------------
write("CsTx2_ubiquitin", complex_yaml(PEPTIDES["CsTx2"], UBIQUITIN, CSTX_SS))
write("scramCsTx2_Kv42", complex_yaml(SCRAM_CSTX2, VSD["Kv42"], CSTX_SS))

print("\nDone. 5 monomers + 4 real complexes + 2 positive + 2 negative = 13 inputs")
print("scram-CsTx2 =", SCRAM_CSTX2)
