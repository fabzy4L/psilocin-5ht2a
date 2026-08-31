"""
ChimeraX video: the redocking-gate finding, visualized. Loads the 5-HT2A
receptor, the 7LD (LSD) crystal pose extracted from 6WGT, and the
pipeline's own redocked top pose for the same ligand -- then orbits the
camera so the 5.11 A gap between them (seed 42; see docs/setup.md,
"Validation: self-redocking") is visible directly, not just quoted as a
number. This is meant to be the opening shot of the video: the actual
finding, not a generic "ligand drops into pocket" clip.

UNTESTED locally: authored without a ChimeraX install in this repo's own
environment (see visualization/README.md for why -- ChimeraX runs on a
different machine than the rest of this pipeline). Before a full
unattended render, open ChimeraX interactively and run the `run(...)`
calls below a few at a time from its command line to confirm each one
behaves as expected on your ChimeraX version, especially the `2dlabels`
and `movie` commands, which vary slightly across releases.

Run headless from the visualization/ directory:
    chimerax --offscreen --exit --script "scripts/01_redock_overlay.py"

Or from inside ChimeraX:
    open scripts/01_redock_overlay.py
"""
from chimerax.core.commands import run

RECEPTOR_PDB = "receptor/6WGT_chainA_repaired.pdb"
CRYSTAL_PDB = "receptor/7LD_crystal_reference.pdb"
REDOCKED_PDB = "poses/7LD_redocked_top_pose.pdb"
OUT_MP4 = "output/01_redock_overlay.mp4"

ORBIT_DEGREES_PER_FRAME = 2
ORBIT_FRAMES = 180  # 2 * 180 = 360 degrees, 6s at 30fps


def open_one(session, path):
    models = run(session, f"open {path}")
    return models[0]


def main(session) -> None:
    run(session, "set bgColor black")
    run(session, "lighting soft")

    receptor = open_one(session, RECEPTOR_PDB)
    run(session, f"hide {receptor.atomspec} atoms")
    run(session, f"cartoon {receptor.atomspec}")
    run(session, f"color {receptor.atomspec} gray70")

    crystal = open_one(session, CRYSTAL_PDB)
    run(session, f"style {crystal.atomspec} stick")
    run(session, f"color {crystal.atomspec} dodgerblue")

    redocked = open_one(session, REDOCKED_PDB)
    run(session, f"style {redocked.atomspec} stick")
    run(session, f"color {redocked.atomspec} orange")

    run(session, f"view {crystal.atomspec},{redocked.atomspec}")

    run(session, 'lighting soft')
    run(session, '2dlabels create title text "Redocking LSD: score reproduced, pose missed" '
                 'xpos 0.04 ypos 0.92 size 28 color white')
    run(session, '2dlabels create legend1 text "crystal pose, PDB 6WGT" '
                 'xpos 0.04 ypos 0.10 size 18 color dodgerblue')
    run(session, '2dlabels create legend2 text "pipeline redock" '
                 'xpos 0.04 ypos 0.06 size 18 color orange')
    run(session, '2dlabels create rmsd text "RMSD 5.11 A vs. 2.0 A pass threshold (seed 42)" '
                 'xpos 0.04 ypos 0.02 size 18 color white')

    run(session, "movie record")
    run(session, f"turn y {ORBIT_DEGREES_PER_FRAME} {ORBIT_FRAMES}")
    run(session, f"wait {ORBIT_FRAMES}")
    run(session, "movie stop")
    run(session, f"movie encode {OUT_MP4} framerate 30 quality high")

    run(session, "2dlabels delete all")
    print(f"Wrote {OUT_MP4}")


main(session)  # noqa: F821 -- `session` is injected by ChimeraX when running a script
