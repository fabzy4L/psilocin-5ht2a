# Docking video

Two short ChimeraX-rendered clips, meant to lead with the actual finding
rather than generic docking footage:

1. **`01_redock_overlay.py`** — the redock-gate result. Receptor pocket,
   6WGT's crystal 7LD (LSD) pose in blue, the pipeline's own redocked pose
   in orange, orbiting so the 5.11 Å gap between them is visible directly.
   This is the hook: it's what "score reproduced, pose missed" (see
   `docs/setup.md`) actually looks like.
2. **`02_ligand_parade.py`** — the comparative screen. Same receptor,
   fixed camera orbit, each of the six ligands' top pose swapped in with
   its name and multi-seed mean affinity, captioned as **provisional**
   since the redock gate above is currently failing.

## Why ChimeraX instead of PyMOL

PyMOL (`pymol-open-source`) is already in `envs/environment.yml` and would
run in this repo's own environment, but ChimeraX generally renders nicer
out of the box for this kind of outreach clip and is what's available —
just on a different machine than the rest of this pipeline. That's why
these are two separate steps:

- **Data prep** (this repo's environment, already run once and committed):
  `scripts/00_export_top_poses.py` converts each ligand's top Vina pose
  from multi-model PDBQT to a plain single-model PDB via OpenBabel, so the
  ChimeraX scripts have no dependency on this repo's WSL/conda docking
  environment. Only re-run it if the docking results change:
  ```bash
  cd visualization/scripts
  python 00_export_top_poses.py
  ```
- **Rendering** (on the machine with ChimeraX):
  ```bash
  cd visualization
  chimerax --offscreen --exit --script "scripts/01_redock_overlay.py"
  chimerax --offscreen --exit --script "scripts/02_ligand_parade.py"
  ```
  Output lands in `output/01_redock_overlay.mp4` and
  `output/02_ligand_parade.mp4` (gitignored — pull the repo, render
  locally, the mp4s aren't meant to be committed).

Splicing the two clips together (title cards, music, trimming) is left to
whatever editor you'd normally use — ChimeraX's `movie encode` just
produces the two raw clips.

## Important: untested

Both scripts were authored without a local ChimeraX install — it isn't in
this repo's environment, by design (see above). The commands used
(`open`, `cartoon`, `style`, `color`, `2dlabels`, `movie record/stop/
encode`, `turn`/`wait`, `view`) are all standard, well-documented
ChimeraX commands, and the scripts use ChimeraX's Python API to capture
each opened model's real `atomspec` rather than guessing numeric model
IDs (`#1`, `#2`, ...), which is the usual way this kind of script breaks.
Still, **do a short interactive dry run first**: open ChimeraX, and paste
the `run(session, "...")` calls from one script a few at a time into its
command line before trusting a full unattended `--offscreen` render — the
most likely trip-ups are `2dlabels` argument names and `movie encode`
flags, which have shifted slightly across ChimeraX versions.

If something doesn't match your version, the fix is almost always a
one-line tweak to the specific `run()` call ChimeraX's command line
reports the error on, not a rewrite.
