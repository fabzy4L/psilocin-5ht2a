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
  scripts/        # .mdp files, run scripts
  analysis/       # RMSD/RMSF/contact analysis scripts and plots
notebooks/        # exploratory analysis (Jupyter)
envs/             # conda environment files
docs/             # writeups, notes, references
```

## Setup

```bash
conda env create -f envs/environment.yml
conda activate psilocin-5ht2a
```

See `docs/setup.md` for GROMACS GPU build notes and CHARMM-GUI walkthrough.

## Status

- [x] Phase 1: docking screen — baseline complete, n=1 seed; see
  `docs/setup.md` for the full pipeline and caveats
  - [x] Receptor fetch + clean (`01_fetch_receptor.py`) — strips BRIL fusion,
    handles missing-loop chain breaks, verified against 6WGT
  - [x] Orthosteric box center sourced from the co-crystallized agonist
    (see `docs/setup.md`)
  - [x] Ligand prep (`02_prep_ligands.py`) — all 6 comparators (psilocin,
    psilocybin, DMT, 5-MeO-DMT, serotonin, LSD) verified end-to-end
  - [x] Receptor → PDBQT — `01b_prep_receptor.py` (meeko) hits a parsing
    edge case on this structure; `01c_prep_receptor_adfr.py` (ADFRsuite)
    works and is the current path — see `docs/setup.md`
  - [x] Vina docking run — LSD −10.1, psilocybin −7.9, serotonin −7.1,
    psilocin −6.9, 5-MeO-DMT −6.9, DMT −6.8 kcal/mol (single-seed; not
    yet validated by redocking or multi-seed averaging — see
    `docs/setup.md` before reading anything into the psilocybin > psilocin
    ordering)
  - [ ] Multi-seed docking + redocking validation (pattern to port from
    SERT's `docking_multiseed.py` / `validation_redock.py`)
- [ ] Phase 2: MD validation
- [ ] Phase 3: signaling bias (stretch)

## Related work

- [`sert-s438t-escitalopram`](https://github.com/fabzy4L/sert-s438t-escitalopram) —
  same docking→MD methodology applied to the S438T SERT variant and
  escitalopram binding. Mechanistically adjacent: SERT (reuptake) and
  5-HT2A (postsynaptic agonism) sit on opposite sides of the serotonergic
  synapse, and serotonin itself is the shared reference ligand in both
  studies. SSRI pretreatment blunting psilocybin's subjective effects is
  an open clinical question this pairing could eventually speak to.
