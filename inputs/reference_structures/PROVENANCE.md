# Reference structures

- **`1LA4.pdb`** — solved NMR structure of SGTx1 (Swartz lab), used as the ground truth
  for the monomer fold validation (TM-align, `scripts/02_boltz2_cofolding/05_validate_fold.py`).
  SGTx1 is the closest solved homolog to CsTx2/CsTx3 used in this study (67.7% sequence
  identity to CsTx3).
- **`native_protx2_vsd2.pdb`** — extracted from PDB 6N4I (ProTx-II bound to the Nav1.7
  VSD2–NavAb chimera), used as the ground-truth interface for the positive co-folding
  benchmark (`scripts/02_boltz2_cofolding/04_validate_benchmark.py`, DockQ-style scoring).
  Chain G = ProTx-II (toxin), chain C = the VSD2 receptor — the original 6N4I chain IDs are
  preserved, unlike every other structure in this repository which uses chain A
  (toxin)/B (receptor).
