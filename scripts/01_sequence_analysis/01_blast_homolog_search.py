#!/usr/bin/env python3
"""
Step 1 — BLASTP search for the closest characterized homologs of each
C. simoni peptide (Cs32.1, CsTx1, CsTx2, CsTx3) against UniProtKB/Swiss-Prot.

Results feed into 02_msa_alignment_figure.py, which selects the top-10
distinct homologs by pairwise identity and builds the alignment figure.

Run (from the repository root):
    python scripts/01_sequence_analysis/01_blast_homolog_search.py

Output:
    generated/blast_results/<name>.xml   (one NCBI BLAST XML report per peptide)
"""
from Bio.Blast import NCBIWWW
import time
import os

SEQUENCES = {
    "Cs32.1": "GCRWMFGACKTTADCCKALACVGTCIWDGS",
    "CsTx1":  "GCRWMFGACKTTADCCKALACVGTCIWDG",
    "CsTx2":  "GCRWMFGACKTTADCCKALACVGTCIWDGTYG",
    "CsTx3":  "GCRWMFGACKTTADCCKALACVGTCIWDGTYGN",
}

RESULTS_DIR = "generated/blast_results"


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    for name, seq in SEQUENCES.items():
        out_file = os.path.join(RESULTS_DIR, f"{name}.xml")
        if os.path.exists(out_file):
            print(f"[skip] {name} — {out_file} already exists")
            continue
        print(f"[blast] {name} ({len(seq)} aa) vs swissprot ...")
        result_handle = NCBIWWW.qblast(
            "blastp",
            "swissprot",
            seq,
            matrix_name="BLOSUM62",
            word_size=2,
            expect=10,
            hitlist_size=100,
        )
        with open(out_file, "w") as f:
            f.write(result_handle.read())
        result_handle.close()
        print(f"  -> saved {out_file}")
        time.sleep(5)  # stay well under NCBI's request-rate guidance

    print("\nAll BLAST searches complete.")


if __name__ == "__main__":
    main()
