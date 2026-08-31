# Handoff: docking video rendering (ChimeraX device)

**Repo:** `psilocin-5ht2a`, latest commit `02ff978` on `master`.
**What's left:** render two ChimeraX clips on the machine that actually
has ChimeraX installed — this repo's own environment doesn't have it.

## Do this

```bash
git pull
cd visualization
chimerax --offscreen --exit --script "scripts/01_redock_overlay.py"
chimerax --offscreen --exit --script "scripts/02_ligand_parade.py"
```

Output: `visualization/output/01_redock_overlay.mp4` and
`02_ligand_parade.mp4` (gitignored, render locally — not meant to be
committed).

**Before the real run:** neither script has been tested against an
actual ChimeraX install. Open ChimeraX and paste the `run(session, "...")`
calls from `01_redock_overlay.py` a few at a time into its command line
first. The most likely trip-ups are `2dlabels` argument names and
`movie encode` flags — both have shifted slightly across ChimeraX
versions. If something errors, ChimeraX's log will point at the exact
`run()` call; fix is almost always a one-line tweak, not a rewrite. Full
details in `visualization/README.md`.

## What each clip is

1. **`01_redock_overlay.py`** — the actual finding, not stock footage.
   Receptor pocket, 6WGT's crystal LSD (ligand code 7LD) pose in blue,
   the pipeline's own redocked pose in orange, camera orbiting so the
   **5.11 Å gap** between them (seed 42, reproducible) is visible
   directly. This is meant to be the opening shot.
2. **`02_ligand_parade.py`** — the comparative screen. Fixed receptor,
   orbiting camera, all 6 ligands cycled through with name + multi-seed
   mean affinity, captioned **"provisional — redock gate currently
   FAILS"** on screen, since that's true and shouldn't get cut for the
   sake of a cleaner video.

Everything the ChimeraX scripts need is self-contained under
`visualization/` (`poses/`, `receptor/`) — no dependency on
`data/processed/` or this repo's WSL/conda docking environment, which
don't exist on the ChimeraX machine. If docking results ever change,
regenerate those inputs from this repo's own environment with
`python visualization/scripts/00_export_top_poses.py` before re-rendering.

## Context, if picking this up cold

This project docks 6 serotonergic ligands (LSD, psilocybin, serotonin,
5-MeO-DMT, psilocin, DMT) against 5-HT2A (PDB 6WGT) with AutoDock Vina.
The headline result so far: the affinity ranking is highly reproducible
across search seeds (SD ≤0.05 kcal/mol), **but** self-redocking 6WGT's
own bound ligand (7LD = LSD — corrected from an earlier mislabel as
25-CN-NBOH, independently verified against RCSB) misses its own crystal
pose by 5.11 Å against a 2.0 Å pass threshold, despite reproducing its
affinity almost exactly. So: score-reproducible, pose-unvalidated. Full
writeup in `docs/setup.md` and `docs/RESEARCH_PROPOSAL.md`
(`docs/research_proposal.html` for a styled read). The sibling project
`sert-s438t-escitalopram` shows the same pattern on its own redock
benchmark (6.63 Å) — not specific to this repo's receptor prep.

`docking/scripts/05_validation_redock.py --strict` is the enforcement
mechanism: it exits non-zero and writes `docking/results/gate_status.json`
when the redock gate fails, so "provisional" in the README is something
checkable, not just asserted.

## Not done yet / open threads

- Ensemble docking design for Aim 1 (the proposal's plan for actually
  clearing the redock gate — flexible-residue docking was tried and
  blocked by tooling gaps in both ADFRsuite and meeko, see
  `docs/setup.md`, "Flexible-residue docking: tooling gap").
- Porting the same `--strict`/`gate_status.json` pattern to
  `sert-s438t-escitalopram`'s own `validation_redock.py`.
- `docs/Redocking LSD.pdf` sits untracked in the repo — confirmed via
  `git log --all --full-history` to have no git history, i.e. a
  forgotten local export, not something orphaned from a branch. Harmless,
  left alone.
