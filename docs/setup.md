# Setup notes

## Environment

```bash
conda env create -f envs/environment.yml
conda activate psilocin-5ht2a
```

## GROMACS with GPU offload

conda-forge's `gromacs` package is CPU-only. For GPU offload on the RTX
5060, build from source with CUDA support:

```bash
# prerequisites: CUDA toolkit matching your driver, cmake >= 3.18
wget https://ftp.gromacs.org/gromacs/gromacs-2024.3.tar.gz
tar xf gromacs-2024.3.tar.gz && cd gromacs-2024.3
mkdir build && cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON \
         -DGMX_GPU=CUDA \
         -DCMAKE_INSTALL_PREFIX=$HOME/gromacs
make -j$(nproc)
make check      # optional but worth it before a long run
make install
source $HOME/gromacs/bin/GMXRC
```

Run with GPU offload for nonbonded and PME:
```bash
gmx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu -ntmpi 1 -ntomp 8
```

Expect roughly 20–50 ns/day for a ~50-100k atom membrane-protein system on
an 8GB card — scale expectations to system size; check `nvidia-smi` during
a short test run to confirm GPU utilization before committing to a long job.

## System building with CHARMM-GUI

1. Go to charmm-gui.org → Membrane Builder.
2. Upload the cleaned receptor PDB from `data/processed/`.
3. Select POPC bilayer (or POPC/cholesterol mix if modeling a more
   physiological membrane).
4. Add the docked ligand pose (from `docking/results/`) into the binding
   site — CHARMM-GUI's ligand reader step handles this, but you'll need to
   generate ligand force field parameters first (CGenFF via the ParamChem
   server, or use GAFF2 with acpype if staying in the Amber ecosystem).
5. Download the resulting GROMACS-ready topology/coordinate files into
   `md/system_prep/`.
6. Equilibration `.mdp` files come bundled with the CHARMM-GUI output —
   copy them into `md/scripts/` and adjust production run length there.

## Notes on the receptor structure

6WGT is a cryo-EM structure and will need inspection before use:
- Check for missing loops (common in ICL3 of GPCRs) — may need modeling
  in if they're near the binding site, otherwise can often be truncated.
- Strip any fusion partners, nanobodies, or antibody fragments used for
  cryo-EM stabilization — these are not part of the biological receptor.
- Check protonation states of key residues (especially the conserved
  Asp3.32) at your simulation pH before force field assignment.

### Confirmed specifics (from actually walking the PDB)

- **BRIL fusion**: single-chain 5-HT2A/cytochrome-b562RIL construct, not a
  separate chain. BRIL replaces native ICL3 in-frame between receptor
  residues 265 and 311, renumbered 1001-1106 in the file (CA trace jumps
  265→1001, then 1106→311). `01_fetch_receptor.py` now drops this range
  by default.
- **Orthosteric pocket center** (co-crystallized agonist 7LD /
  25-CN-NBOH, chain A, resi A1201, centroid of 24 heavy atoms):
  `(25.115, 40.909, 54.225)`. Already wired into `03_run_vina.py`.
- **Missing loop density** (beyond the BRIL splice): gaps at residues
  182-183 and 216-218. `01_fetch_receptor.py` inserts `TER` at every
  numbering discontinuity so downstream tools don't try to bond across
  them.
- **Truncated side chains near the pocket**: at ~3-4 Å resolution, ~80
  residues have partially unresolved side chains. Most are safely far
  from the binding site and can be auto-repaired/deleted, but nine sit
  within 5 Å of the docking box and need real attention rather than
  silent deletion: **A:89, 135, 136, 146, 147, 348, 350, 351, 356**.
  `01b_prep_receptor.py --bad-res-radius` flags this distinction
  automatically.

### Receptor → PDBQT: meeko fails, ADFR works

`mk_prepare_receptor` (meeko, `01b_prep_receptor.py`) does strict
RDKit-template-based residue parsing. Even after repairing missing
side-chain atoms with PDBFixer, it still fails on this structure — an
`AtomValenceException` at the true C-terminus (TYR 399) that looks like
a geometry/bond-perception edge case in meeko's own parsing, not
something worth chasing further blind.

**Working path**: ADFRsuite's `prepare_receptor` (`01c_prep_receptor_adfr.py`,
wrapping the classic AutoDockTools `prepare_receptor4.py`) does
geometry-based bond building instead of strict template matching, and
handles this structure fine. ADFRsuite isn't on PyPI/conda-forge
directly — its official installer is a legacy Tcl/Tk binary that won't
run headless — but the community `hcc` conda channel has a working
build:

```bash
conda create -n adfr -c hcc -c conda-forge adfr-suite
conda activate adfr
# no native Windows build — run this from WSL2 (Ubuntu) on Windows
```

Full pipeline that produced the current `docking/results/summary.json`:

```bash
python 01_fetch_receptor.py --pdb-id 6WGT --out ../../data/raw/6WGT.pdb
python 01a_repair_sidechains.py \
    --in-pdb ../../data/processed/6WGT_chainA.pdb \
    --out-pdb ../../data/processed/6WGT_chainA_repaired.pdb \
    --true-terminus-resnum 399
conda run -n adfr python 01c_prep_receptor_adfr.py \
    --in-pdb ../../data/processed/6WGT_chainA_repaired.pdb \
    --out-pdbqt ../../data/processed/6WGT_chainA_adfr.pdbqt
python 02_prep_ligands.py
python 03_run_vina.py --receptor ../../data/processed/6WGT_chainA_adfr.pdbqt
```

`01a_repair_sidechains.py` requires `pdbfixer` + `openmm` (pip-installable;
not yet added to `envs/environment.yml` — they only need to be present
in whichever env runs that one script, not necessarily the main one).

One cosmetic warning is expected and harmless: `prepare_receptor` can't
find Gasteiger parameters for TYR 399's `OXT` (the genuine but
structurally irrelevant terminus of this truncated crystallization
construct) and assigns it zero charge. That residue isn't part of the
binding pocket, so it doesn't affect docking.

### Phase 1 baseline result (rigid single-structure docking, n=1 seed)

| Ligand | Best affinity (kcal/mol) |
|---|---:|
| LSD | −10.1 |
| psilocybin | −7.9 |
| serotonin | −7.1 |
| psilocin | −6.9 |
| 5-MeO-DMT | −6.9 |
| DMT | −6.8 |

LSD's outsized affinity matches its well-documented extreme 5-HT2A
residence time. Psilocybin scoring *better* than psilocin here is
probably a rigid-docking artifact (its phosphate group is finding
favorable polar contacts that wouldn't survive an induced-fit or MD
treatment) rather than a real result — psilocybin is a prodrug and is
not generally thought to engage 5-HT2A directly as well as psilocin
does. Don't read anything into this ranking until: (1) multi-seed
docking for a proper SD (see SERT's `docking_multiseed.py` for the
pattern), and (2) redocking validation against the co-crystallized 7LD
ligand to confirm the pipeline reproduces a known pose (again, SERT's
`validation_redock.py` is the template to port over).
