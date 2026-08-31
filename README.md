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

- [ ] Phase 1: docking screen
- [ ] Phase 2: MD validation
- [ ] Phase 3: signaling bias (stretch)
