# psilocin-5ht2a

Structure-based modeling of psilocin (the active metabolite of psilocybin) engaging
the 5-HT2A receptor: comparative docking against native serotonin and related
tryptamines, followed by GPU-accelerated MD validation of the top binding pose.

Scoped to run end-to-end on a single consumer GPU (developed against an RTX 5060 8GB).

## Project phases

1. **Docking screen** (CPU-only, cheap)
   Dock psilocin, serotonin, DMT, 5-MeO-DMT, and LSD against 5-HT2A (PDB 6WGT)
   with AutoDock Vina. Compare poses and scores, identify likely selectivity-driving
   residues.

2. **MD validation** (GPU-accelerated, local)
   Take the best-scoring psilocin pose, embed in a POPC bilayer via CHARMM-GUI,
   run 100–200 ns in GROMACS with GPU offload. Check RMSD/RMSF for pose stability
   and track key contacts (Asp3.32 salt bridge, TM5 Ser/His network).

3. **Signaling bias (stretch)**
   Compare TM6 outward displacement in the psilocin trajectory against a
   Gq-biased agonist reference structure, as a rough proxy for active-state
   conformation.

## Repo layout

```
data/
  raw/            # downloaded PDB structures, unmodified
  processed/      # cleaned/protonated structures ready for docking or MD
  ligands/        # ligand SDF/PDB/SMILES files
docking/
  scripts/        # prep + Vina run scripts
  results/        # docking logs, scored poses
md/
  system_prep/    # CHARMM-GUI outputs, topology/coordinate files
  scripts/        # staged MD pipeline: import -> minimize -> equilibrate ->
                   # production -> analysis (docs/METHODOLOGY.md) --
                   # scaffolded, not yet executed against real GROMACS output
  analysis/       # RMSD/RMSF/contact analysis output (written by 05_analysis.py)
notebooks/        # exploratory analysis (Jupyter)
envs/             # conda environment files
docs/             # writeups, notes, references
visualization/    # ChimeraX video scripts (see visualization/README.md)
tests/            # pytest coverage for md/scripts/ gate logic
```

## Setup

```bash
conda env create -f envs/environment.yml
conda activate psilocin-5ht2a
```

See `docs/setup.md` for GROMACS GPU build notes and CHARMM-GUI walkthrough.

## Status

> **Redock gate: FAIL.** `05_validation_redock.py --strict` currently
> exits non-zero — the pipeline does not reproduce 6WGT's own
> co-crystallized ligand's pose (RMSD 5.11 Å vs. a 2.0 Å pass
> threshold, seed 42, reproducible). Every comparative
> affinity number below (single-seed and multi-seed) is **provisional**
> until this gate passes — see "The redock gate" below before reading
> anything into the ranking, especially psilocybin > psilocin.

- [x] Phase 1a — receptor/box setup (`docs/setup.md`)
  - [x] Receptor fetch + clean (`01_fetch_receptor.py`) — strips BRIL fusion,
    handles missing-loop chain breaks, verified against 6WGT
  - [x] Orthosteric box center sourced from the co-crystallized agonist
    (7LD = LSD, Kim et al. *Cell* 2020 — see correction in `docs/setup.md`)
  - [x] Receptor → PDBQT — `01b_prep_receptor.py` (meeko) hits a parsing
    edge case on this structure; `01c_prep_receptor_adfr.py` (ADFRsuite)
    works and is the current path — see `docs/setup.md`
- [x] Phase 1b — **redock gate** (`05_validation_redock.py`, ported from
  SERT's `validation_redock.py`) — extracts 6WGT's own bound ligand and
  redocks it into the prepared receptor/box. **Currently FAILS**
  (RMSD 5.11 Å) despite reproducing the ligand's own blind-screen
  affinity almost exactly. Run this — and get it to PASS — before
  trusting Phase 1c's ranking as anything more than "the scoring
  function likes this ligand's features." `sert-s438t-escitalopram`'s
  own redock gate fails too (6.63 Å) — this isn't a psilocin-5ht2a-
  specific quirk, see `docs/setup.md`.
- [x] Phase 1c — comparative screen (**provisional**, gated by 1b)
  - [x] Ligand prep (`02_prep_ligands.py`) — all 6 comparators (psilocin,
    psilocybin, DMT, 5-MeO-DMT, serotonin, LSD) verified end-to-end
  - [x] Vina docking run — LSD −10.1, psilocybin −7.9, serotonin −7.1,
    psilocin −6.9, 5-MeO-DMT −6.9, DMT −6.8 kcal/mol (single seed)
  - [x] Multi-seed confirmation (`04_docking_multiseed.py`) — 5-seed
    means match the single-seed baseline within ≤0.05 kcal/mol SD, so
    the ranking is **not** search-noise. It is, however, ungated: see
    the FAIL above before treating it as a pose-level result.
- [ ] Phase 1d — flexible/ensemble redocking to try to clear the gate
  (see `docs/setup.md`, "Flexible-residue docking: tooling gap" —
  meeko and this ADFRsuite build both currently block the natural
  approach; ensemble docking against MD-generated conformers is the
  likely path, folding this into Phase 2)
- [ ] Phase 2: MD validation — pipeline scaffolded (`md/scripts/`,
  `docs/METHODOLOGY.md`, `tests/`), not yet executed: no CHARMM-GUI
  export or working GROMACS build exists in this repo yet, see
  `docs/setup.md`
- [ ] Phase 3: signaling bias (stretch)

## Related work

- [`sert-s438t-escitalopram`](https://github.com/fabzy4L/sert-s438t-escitalopram) —
  same docking→MD methodology applied to the S438T SERT variant and
  escitalopram binding. Mechanistically adjacent: SERT (reuptake) and
  5-HT2A (postsynaptic agonism) sit on opposite sides of the serotonergic
  synapse, and serotonin itself is the shared reference ligand in both
  studies. SSRI pretreatment blunting psilocybin's subjective effects is
  an open clinical question this pairing could eventually speak to.
