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
- **Orthosteric pocket center** (co-crystallized agonist 7LD, chain A,
  resi A1201, centroid of 24 heavy atoms): `(25.115, 40.909, 54.225)`.
  Already wired into `03_run_vina.py`. **Correction**: 7LD is LSD
  (local HETNAM: "LYSERGIC ACID DIETHYLAMIDE"), not 25-CN-NBOH as this
  note previously said. Independently confirmed against RCSB's own
  entry page for 6WGT, which names 7LD
  "(8alpha)-N,N-diethyl-6-methyl-9,10-didehydroergoline-8-carboxamide"
  (LSD's systematic name) and cites Kim, K.L., Che, T., Panova, O.,
  DiBerto, J.F., Lyu, J., Krumm, B.E., et al., "Structure of a
  Hallucinogen-Activated Gq-Coupled 5-HT2A Serotonin Receptor," *Cell*
  182:1574–1588 (2020), doi:10.1016/j.cell.2020.08.024. This means LSD
  in the Phase 1 comparator set is both the co-crystallized ligand and
  one of the six screened tryptamines, which is what makes the
  self-redocking check in `05_validation_redock.py` meaningful.
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

### Validation: self-redocking of the co-crystallized ligand (LSD/7LD)

`05_validation_redock.py` extracts 7LD's crystal pose directly from
`data/raw/6WGT.pdb`, redocks it into the prepared receptor with the same
box used for the Phase 1 screen, and compares the top-scoring redocked
pose to the crystal pose by heavy-atom RMSD (name-matched, 2.0 Å
pass/fail threshold — standard convention). Since 7LD is LSD (see
correction above), this is a true self-redock: the ligand being
validated against is also one of the six Phase 1 comparators.

**Result: PASS (seed 42, top-scoring pose).** Best redocked affinity
(−10.073 kcal/mol) reproduces the blind-screen LSD score
(−10.1 kcal/mol), and the pose itself has RMSD **0.780 Å** from the
crystal pose — well under the 2.0 Å threshold. The ergoline ring
nitrogen sits 3.2–3.6 Å from Asp155 (Asp3.32), matching the crystal
structure's own 2.8–3.3 Å salt-bridge geometry almost exactly.

An earlier run of this script reported FAIL at RMSD 5.107 Å. That number
was real, but it was computed against the wrong pose: `compute_rmsd()`
read every `MODEL` in the 10-pose docked PDBQT into a plain dict keyed
by atom name with no `MODEL`/`ENDMDL` awareness, so pose 10 (the
worst-ranked pose) silently overwrote pose 1 (the best-scoring one) by
the time RMSD was computed — the 5.107 Å figure was pose 10's RMSD, not
pose 1's. The per-atom deviation analysis that ruled out a
symmetric-group name-matching artifact (roughly uniform 3.15–8.10 Å
across all 24 heavy atoms) was real too, it was just describing pose
10's deviation, not evidence about pose 1. Fixed by stopping the parse
at the first `ENDMDL`; full writeup in `docs/REDOCK_BUG_HANDOFF.md`. No
re-docking was needed — the correct pose was already sitting in
`docking/results/7LD_redocked.pdbqt`.

This means the box/receptor/scoring protocol **is** validated against a
known pose. It does not, by itself, validate the comparative ranking
across the six screened ligands (Phase 1c) — see Phase 1d below for that
open question.

Use `--strict` to make this a real pipeline gate (non-zero exit on
FAIL, plus a machine-readable `docking/results/gate_status.json`)
rather than just a logged warning:

```bash
python 05_validation_redock.py --strict
```

**Not independent corroboration across projects — same bug, twice.**
`sert-s438t-escitalopram`'s own `validation_redock.py` (68P re-dock into
5I6Z) was also reported as failing, RMSD 6.63 Å. This was cited here as
weak independent evidence that rigid single-structure Vina docking on a
GPCR orthosteric pocket generally can't recover the true pose. That
citation no longer holds: SERT's script had the *identical* last-model-
wins bug (BioPython's `struct.get_atoms()` iterating every model instead
of `struct[0]`), and once fixed, its redock also PASSES — RMSD 1.515 Å
for pose 1. Both "independent" failures were the same one bug, fixed the
same way, in two scripts derived from a common ancestor. Neither repo's
redock gate has actually demonstrated a rigid-docking pose-recovery
failure; both pass.

### Flexible-residue docking: tooling gap

Phase 1d (see README) — Vina flexible-residue redocking, treating the 9
pocket-adjacent side chains (A:89 ILE, 135 ILE, 136 LEU, 146 LYS, 147
LEU, 348 ILE, 350 LYS, 351 GLU, 356 ASP) as rotatable — is now optional
follow-up work rather than a requirement blocking Phase 2: the redock
gate above already validates that the rigid protocol can recover a known
pose. It's still worth doing to test whether the psilocybin > psilocin
ranking (Phase 1c) survives side-chain flexibility around the pocket,
since that's a separate question the redock gate doesn't answer. Two
attempts at generating the required rigid/flexible PDBQT split both hit
walls:

- **ADFRsuite** (`prepare_receptor`, already working for the rigid
  case) needs a companion `prepare_flexreceptor` to do the split. The
  `hcc` conda build used here doesn't ship it — `$CONDA_PREFIX/bin` has
  `prepare_receptor`/`prepare_ligand` only, no `prepare_flexreceptor*`,
  and there's no `AutoDockTools/Utilities24` tree bundled either.
- **meeko** (`mk_prepare_receptor --flexres`) does support this in
  principle, but it still has to residue-template-match the *entire*
  polymer first — same failure mode as `01b_prep_receptor.py`. Tested
  directly: even with residue 399 (the C-terminus that broke 01b)
  stripped out, it fails on a **different** residue, `AtomValenceException`
  on ILE A:135 — one of the exact 9 pocket residues this step needs to
  make flexible. Not a one-off edge case; meeko's template matching is
  fragile broadly across this repaired structure.

Hand-authoring the flexible PDBQT torsion trees directly (bypassing
both tools) was considered and rejected: getting per-residue chi-angle
branching wrong silently produces a chemically invalid flexible residue
with no obvious error, which is a bad trade for a docking result headed
into a research proposal. The safer path is either (a) install full
MGLTools/AutoDockTools for a real `prepare_flexreceptor4.py`, or
(b) skip flexible side chains and go straight to ensemble docking
against a handful of receptor conformers (e.g. short restrained MD or
normal-mode-perturbed structures) — architecturally simpler and folds
naturally into the Phase 2 MD work already planned. (b) is the current
default plan.

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
does.

### Multi-seed result (n=5 seeds, `04_docking_multiseed.py`)

| Ligand | mean (kcal/mol) | SD | n |
|---|---:|---:|---:|
| LSD | −10.045 | 0.021 | 5 |
| psilocybin | −7.837 | 0.033 | 5 |
| serotonin | −7.083 | 0.019 | 5 |
| 5-MeO-DMT | −6.982 | 0.018 | 5 |
| psilocin | −6.937 | 0.021 | 5 |
| DMT | −6.872 | 0.054 | 5 |

SD is ≤0.05 kcal/mol for every ligand — an order of magnitude below
Vina's commonly-cited ~0.5 kcal/mol scoring noise floor. **This rules
out search stochasticity as the explanation for psilocybin > psilocin**:
the ranking is highly reproducible *given this box, this receptor
conformation, and this scoring function*. It says nothing about whether
that ranking reflects real binding-pose energetics — see the
self-redocking result immediately above, which suggests it may not.
