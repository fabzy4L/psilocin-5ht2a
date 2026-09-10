# Phase 2 (MD) Methodology

**Scope:** GROMACS system building, equilibration, production, and analysis
for the psilocin/5-HT2A complex. Docking (Phase 1) methodology lives in
`docs/setup.md`; this document picks up where that one leaves off, once a
validated docked pose exists (see `docs/REDOCK_BUG_HANDOFF.md` for the
current state of that validation).

**Status: scaffolded, not yet executed.** No CHARMM-GUI export or working
GROMACS build exists in this repo yet — see `docs/setup.md`'s "GROMACS
with GPU offload" and "System building with CHARMM-GUI" sections for the
prerequisite build/download steps. Everything below describes what
`md/scripts/` does once those exist, not results that exist today.

## 1. Pipeline overview

```
Docked psilocin pose (docking/results/, Phase 1)
         |
         v
CHARMM-GUI Membrane Builder (interactive, browser -- docs/setup.md)
         |
         v
01_import_charmm_gui.py   -- validate + stage the download
         |
         v
02_minimize.py             -- energy minimization, gated on Fmax convergence
         |
         v
03_equilibrate.py --phase nvt   -- gated on temperature stability
03_equilibrate.py --phase npt   -- gated on pressure/density stability
         |
         v
04_production.py           -- 100-200 ns GPU-offloaded production
         |
         v
05_analysis.py              -- RMSD/RMSF, Asp3.32 salt bridge, gated on
                                ligand-pose RMSD plateauing
```

## 2. Design pattern: every automatable stage is a gate

This phase follows the same convention `docking/scripts/05_validation_redock.py`
established: each stage that produces a checkable result writes a
human-readable `*_report.txt`, a machine-readable `gate_status.json`
(`{"gate": "PASS"|"FAIL", ...criteria}`), and accepts `--strict` to exit
non-zero on FAIL so a stage can't silently feed a bad result to the next
one.

That convention exists *because* of a real incident, not preemptively:
`docking/scripts/05_validation_redock.py` reported a FAIL for weeks that
turned out to be a parsing bug (it graded the worst of 10 Vina poses
instead of the best one) — see `docs/REDOCK_BUG_HANDOFF.md`. Two things
changed here in direct response to that:

1. **The report/gate plumbing is factored into `md/scripts/gate_utils.py`**
   instead of being re-implemented per script (as the docking phase's
   scripts each did independently). One tested implementation is safer
   than five duplicated ones.
2. **Every gate has a real unit test** (`tests/`) against a synthetic
   fixture — a fake `gmx` log, a fake `.xvg`, a synthetic RMSD series —
   not just "the script ran without an exception." A test asserting
   `is_plateaued()` correctly distinguishes a settled trajectory tail
   from a still-drifting one, for instance, is exactly the kind of check
   that would have caught the redock bug's "grades the wrong
   frame/pose" failure mode on day one.

## 3. Stage details

### 3.1 Import (`01_import_charmm_gui.py`)

CHARMM-GUI's Membrane Builder step is interactive (browser upload of the
receptor + docked ligand, POPC bilayer selection — see `docs/setup.md`)
and can't be scripted. This stage validates what comes back from that
download: a topology (`topol.top`), an index (`index.ndx`), initial
coordinates, and at least one `.mdp` file, before anything downstream
touches it. A partial or wrong download fails here with a specific
missing-file list, not three stages later as a cryptic `gmx grompp` error.

### 3.2 Minimization (`02_minimize.py`)

Steepest-descent energy minimization via `gmx grompp` + `gmx mdrun`.
Gated on the *actual* final max force parsed out of the `mdrun` log
(default threshold: Fmax < 1000 kJ/mol/nm), not on `mdrun` merely exiting
0 — a run can hit its step limit without converging and still exit
cleanly.

### 3.3 Equilibration (`03_equilibrate.py`)

Two phases, run separately: NVT (gated on temperature stability) then NPT
(gated on pressure and density stability). Stability is checked as the
standard deviation over the last `--tail-ps` of the run (default 100 ps)
against a per-observable tolerance (2 K temperature, 50 bar pressure, 5
kg/m³ density by default) — a run that's still drifting when it ends
shouldn't be treated as equilibrated input for production.

### 3.4 Production (`04_production.py`)

Long GPU-offloaded `gmx mdrun` (100–200 ns per `docs/setup.md`'s scoping
for an RTX 5060 8GB). Not pass/fail gated the way the earlier stages are
— there's no single number that means "production succeeded" at this
point, that's what the analysis stage's plateau gate is for. Supports
`--smoke-test` (a few ps) so GPU utilization can be checked via
`nvidia-smi` before committing to a multi-day run, per `docs/setup.md`'s
own advice.

### 3.5 Analysis (`05_analysis.py`)

Loads the trajectory with MDAnalysis (already in `envs/environment.yml`)
rather than hand-parsing GRO/XTC frames — multi-frame binary trajectory
parsing is exactly the class of problem a maintained library should own,
and hand-rolling it would repeat the redock bug's mistake at a different
layer. Computes:

- **Ligand pose RMSD** over the trajectory, gated on whether the tail has
  *plateaued* (std dev within tolerance over the last `--tail-frac` of
  frames), not just "RMSD was computed."
- **Asp155 (Asp3.32) salt-bridge distance** — the ligand's basic ring
  nitrogen to the nearer of Asp155's two carboxylate oxygens. The
  geometry function (`salt_bridge_distance()`) generalizes the same
  check done by hand against the crystal structure during the redock-bug
  investigation (2.8–3.3 Å in the 6WGT crystal pose, reproduced at
  3.2–3.6 Å by the corrected redocked pose — see
  `docs/REDOCK_BUG_HANDOFF.md`).

TM5 Ser/His contact tracking (per `README.md`'s Phase 2 scope) is not yet
implemented — `05_analysis.py` currently covers ligand RMSD and the
Asp3.32 salt bridge only; extending `run_analysis()` with an additional
`MDAnalysis` distance selection is the natural next step once a real
trajectory exists to validate the selection against.

## 4. What's still a scaffold vs. what's real

Everything in `md/scripts/` is real, runnable code with real parsing/gate
logic (see `tests/`), but none of it has been run against real GROMACS
output — there isn't any yet. Specifically not done in this pass:

- No `gmx` build exists in this repo's tracked state (`docs/setup.md`'s
  CUDA source build is a local, un-tracked prerequisite).
- No CHARMM-GUI export has been generated or imported.
- TM5 Ser/His contacts (analysis stage).
- Phase 3 (signaling-bias / TM6 displacement proxy) has no scaffold yet
  — it's a stretch goal per `README.md` and depends on a reference
  Gq-biased agonist structure not yet identified.
