# Redock gate bug — handoff notes (2026-09-09)

## What was found

The "redock gate FAIL" in both `psilocin-5ht2a` and `sert-s438t-escitalopram`
was a **parsing bug, not a real docking failure**. Both validation scripts
read every pose in a multi-`MODEL` docked PDBQT into a plain `dict` keyed by
atom name, with no `MODEL`/`ENDMDL` boundary awareness. Since every pose
reuses the same atom names, later (worse-ranked) poses silently overwrote
pose 1 — so the reported RMSD was actually the *worst* pose's RMSD, not the
best-scoring pose Vina reported the affinity for.

| Repo | Old reported RMSD (bug, worst pose) | Real RMSD (pose 1, fixed) | Verdict |
|---|---:|---:|---|
| `psilocin-5ht2a` (7LD/LSD self-redock) | 5.107 A (pose 10/10) | **0.780 A** | PASS |
| `sert-s438t-escitalopram` (68P/escitalopram redock) | 6.632 A (pose 5/5) | **1.515 A** | PASS |

`psilocin-5ht2a`'s pose 1 also reproduces the crystal salt bridge almost
exactly: the ergoline ring N sits 3.2-3.6 A from Asp155 (Asp3.32), vs. 2.8-3.3
A in the crystal structure.

## What was changed (uncommitted, local only, in both repos)

- **`psilocin-5ht2a/docking/scripts/05_validation_redock.py`** —
  `get_heavy_atom_coords_by_name()` now breaks at the first `ENDMDL` line.
- **`psilocin-5ht2a/docking/results/validation_report.txt`** and
  **`gate_status.json`** (gitignored) — regenerated from the *existing*
  docked poses (no re-docking needed) → gate now reads **PASS, 0.780 A**.
- **`sert-s438t-escitalopram/scripts/validation_redock.py`** —
  `get_ligand_coords_by_name()` now restricts to `struct[0]` (first model)
  instead of `struct.get_atoms()` (all models).
- **`sert-s438t-escitalopram/output/validation_report.txt`** — regenerated
  the same way → gate now reads **PASS, 1.515 A**.

Neither repo's changes are committed yet — both are sitting as unstaged
working-tree edits (`git status` in each repo shows only these two modified
files).

## Still open — narrative cleanup not yet done

These still describe the old, incorrect "rigid docking hit a fundamental
limitation" conclusion and need revisiting once you're comfortable with the
fix:

- `psilocin-5ht2a/README.md` — status checklist, Phase 1b/1c framing,
  "psilocybin > psilocin" caveat language.
- `psilocin-5ht2a/docs/setup.md` — the "Validation: self-redocking" section
  (result writeup) and the "not a psilocin-5ht2a-specific artifact" claim
  (built on the *same bug* firing in both repos, not independent evidence).
- `psilocin-5ht2a/docs/RESEARCH_PROPOSAL.md` + `docs/research_proposal.html`
  — the whole proposal is framed around the score/pose dissociation that no
  longer holds.
- The Phase 1d flexible-docking pivot (README, `docs/setup.md`) was a
  reaction to the false failure — worth deciding whether it's still needed
  now that the redock gate actually passes.
- Commit messages already pushed to GitHub in both repos state the FAIL
  result as fact (e.g. psilocin-5ht2a's `d850aec`/now `a98c302` "Turn the
  redock check into a real gate; confirm citation; log the flex-docking dead
  end") — not urgent to rewrite, but worth knowing the history says FAIL
  where the corrected code now says PASS.

## Also worth a look, not yet investigated

- Whether `sert-s438t-escitalopram`'s other output files
  (`docking_results_report.txt`, `multiseed_docking_report.txt`,
  `structural_audit_report.txt`) or scripts have similar multi-model
  parsing patterns worth auditing.
- Independent lead surfaced during this debugging (now probably moot given
  the pose-1 pass, but recorded in case the picture changes again): 6WGT
  chain A is missing 3 residues (Gln216-Asp217-Asp218, confirmed via
  `REMARK 465`) right in the ECL2 loop between the helix ending at Phe213 and
  the beta-hairpin at Phe222/Lys223-Ser226/Cys227 — structurally the same
  region reported to form a ligand-clamping "lid" in the homologous
  5-HT2B/LSD structure (residues 207-214, Wacker et al. 2017). Not needed to
  explain the redock result now that pose 1 passes, but could matter later
  for MD stability/RMSF if that loop is ever added back in via homology/loop
  modeling.
