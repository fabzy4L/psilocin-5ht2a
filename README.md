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

> **Redock gate: PASS.** `05_validation_redock.py --strict` reproduces
> 6WGT's own co-crystallized ligand's pose (RMSD 0.780 Å vs. a 2.0 Å pass
> threshold, seed 42, top-scoring pose) and its blind-screen affinity
> (−10.073 vs. −10.1 kcal/mol). The ergoline ring nitrogen sits 3.2–3.6 Å
> from Asp155 (Asp3.32), matching the crystal structure's 2.8–3.3 Å salt
> bridge. An earlier FAIL report (RMSD 5.107 Å) was a parsing bug, not a
> real docking failure — see `docs/REDOCK_BUG_HANDOFF.md` for the full
> writeup.

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
  redocks it into the prepared receptor/box. **PASSES** (RMSD 0.780 Å)
  and reproduces the ligand's own blind-screen affinity almost exactly.
  `sert-s438t-escitalopram`'s own redock gate passes too (1.515 Å) once
  the same parsing bug is fixed there — see `docs/setup.md` and
  `docs/REDOCK_BUG_HANDOFF.md`.
- [x] Phase 1c — comparative screen
  - [x] Ligand prep (`02_prep_ligands.py`) — all 6 comparators (psilocin,
    psilocybin, DMT, 5-MeO-DMT, serotonin, LSD) verified end-to-end
  - [x] Vina docking run — LSD −10.1, psilocybin −7.9, serotonin −7.1,
    psilocin −6.9, 5-MeO-DMT −6.9, DMT −6.8 kcal/mol (single seed)
  - [x] Multi-seed confirmation (`04_docking_multiseed.py`) — 5-seed
    means match the single-seed baseline within ≤0.05 kcal/mol SD, so
    the ranking is **not** search-noise. The redock gate now passing
    means the box/receptor/scoring protocol is validated against a known
    pose — it does **not** by itself confirm the psilocybin > psilocin
    ordering reflects real binding energetics rather than a rigid-receptor
    scoring artifact on psilocybin's phosphate group; that's still an
    open question (see Phase 1d).
- [x] Phase 1d — flexible-residue redocking: tooling gap closed and
  validated on 7LD/LSD (0.958 Å PASS, `05b_validation_redock_flex.py`),
  then the 6-ligand comparator set re-screened flexibly
  (`06_flex_screen.py`, `docking/results/flex_screen_report.txt`). **The
  psilocybin > psilocin ordering survives flexibility** at these 9
  residues — gap widens slightly (1.08 kcal/mol vs ~0.9 rigid
  multi-seed) rather than closing. This rules out "rigid-receptor
  artifact recoverable by side-chain flexibility" as the explanation; it
  does not rule out a solvation/desolvation effect (psilocybin's
  phosphate group), which needs MD (Phase 2) with explicit water to
  test, not more docking variants. Single-seed so far — not yet
  multi-seed confirmed the way the rigid screen was.
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
