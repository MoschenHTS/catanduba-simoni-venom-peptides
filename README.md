# Catanduba simoni venom peptides — Kv4 gating-modifier computational pipeline

Code, sequence data, and result summaries supporting the paper *"Novel Peptides
from Catanduba simoni Venom Enhance Depolarization-Dependent Potassium Currents
in Cardiomyocytes."*

This repository covers the two computational parts of the study: the sequence/
homology analysis (BLAST search, multiple sequence alignment) and the structural
biology campaign (Boltz-2 co-folding + membrane molecular dynamics simulations of
the toxins against the Kv4.2/Kv4.3 voltage-sensor domain). Wet-lab methods (venom
fractionation, mass spectrometry, Edman sequencing, patch-clamp electrophysiology)
are described in the paper's Methods section and are not scripted here.

## Layout

```
peptide_sequences.csv    the four purified peptides (IDs, sequences, masses)
scripts/                 all analysis code, organized by pipeline stage (below)
inputs/                  force field, pre-equilibrated bilayer, MD parameters, reference PDBs
env/                     exported conda/micromamba environment specs + key tool versions
results/                 final summary tables (see below) — no images; the published
                         figures are in the paper itself
```

## Pipeline stages (`scripts/`)

```
01_sequence_analysis/     BLASTP search + Figure 2 (MSA vs closest homologs)
02_boltz2_cofolding/       Boltz-2 structure prediction, evaluation, benchmarks;
                           renders Figures 3-5
03_membrane_md_build/      membrane system construction (orient, embed, equilibrate)
04_membrane_md_analysis/   trajectory analysis (RMSD, H-bonds, membrane partitioning,
                           MM-GBSA)
05_figures_and_statistics/ statistics tests, power analysis, Figure 6 and
                           Supplementary Figure S3
```

Scripts are meant to be run from the repository root, e.g.:
```bash
python scripts/01_sequence_analysis/01_blast_homolog_search.py
chimerax scripts/02_boltz2_cofolding/10_render_figure3_superposition.cxc
```

## Results (`results/`)

Final, human-readable summary tables — the numbers behind Table S1 and Figures
5-6 — without the multi-GB intermediate data they were computed from:

| File | Contents |
|---|---|
| `boltz2_cofolding_summary.txt` | ipTM/convergence/interface-PAE per complex, multi-seed reproducibility, scramble-null statistics, the ProTx-II/6N4I benchmark, and TM-align/disulfide fold validation |
| `interface_pharmacophore_occupancy.txt` | Per-residue interface contact frequency (50-model ensemble) for all four CsTx complexes plus the SGTx1 and ProTx-II references |
| `md_rmsd_paddle_contact.csv` | Per-replica ligand RMSD, paddle-contact occupancy, and R3 salt-bridge occupancy (n=8 replicas x 4 complexes) |
| `md_hbond_occupancy.csv` | Per-replica Trp4-Thr22 backbone/side-chain H-bond occupancy |
| `md_mmgbsa_binding_energy.csv` | Per-replica single-trajectory MM-GBSA energy components |

## Primary data (not included)

The Boltz-2 prediction ensembles (`.pdb`/confidence files, ~600 samples), the raw
per-frame MD analysis outputs, and the 200 ns x 32-replica production trajectories
(27.7 GB) are **not** included in this repository — only the code that produced
them and the final summary tables in `results/` above, which are sufficient to see
every reported number and statistic without re-running anything. The primary data
is available from the corresponding author on request.

Scripts that operate on that primary data (everything under `scripts/02_.../04_...`
except the figure-composite steps, plus `scripts/05_figures_and_statistics/01_collect_data.py`)
read it via environment variables (`CSIMONI_BASE`, `CSIMONI_MD_BASE`,
`CSIMONI_TIER1_ROOT`, `BOLTZ_PRED_ROOT`, `BOLTZ_REF_ROOT` — set the ones each script
needs to wherever you've placed the data) and write their own output — predictions,
post-MD structures, regenerated figures — into a local `generated/` folder at the
repository root, which is git-ignored.

## Third-party data

Two categories of input data in `inputs/` originate outside this project and
retain their original terms rather than this repository's license:

- **Reference protein structures** (`inputs/reference_structures/`, from PDB
  entries 1LA4 and 6N4I) — RCSB PDB coordinate data carries no copyright
  restriction. See `PROVENANCE.md` in that folder for how each was used.
- **CHARMM36m force field files** (`inputs/forcefield/`, `inputs/bilayer/POPC.itp`)
  — the standard GROMACS-format port of the MacKerell-lab CHARMM36m force field
  (Huang et al. 2017), freely available for research use.

## License

Code and original data in this repository are released under the [MIT
License](LICENSE), except for the third-party inputs listed above.

## Citation

If you use this code or data, please cite the paper (full citation to be added
on publication). For citing this repository specifically, see [CITATION.cff](CITATION.cff).
