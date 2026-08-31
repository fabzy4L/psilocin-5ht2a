"""
ChimeraX video: comparative 6-ligand "parade" through the 5-HT2A
orthosteric pocket. Receptor stays fixed under one continuous camera
orbit; each ligand's top Vina pose is swapped in in turn, labeled with
its name and multi-seed mean affinity (see docs/setup.md /
docking/results/multiseed_report.txt). Intended as the second shot of
the video, after 01_redock_overlay.py's finding -- the caveat label is
deliberate: this ranking is provisional until the redock gate passes.

UNTESTED locally: see the note at the top of 01_redock_overlay.py --
same caveat applies here (no ChimeraX in this repo's own environment).
Sanity-check the `run(...)` calls interactively on your ChimeraX
install before a full unattended render, particularly the model-open/
close loop below.

Run headless from the visualization/ directory:
    chimerax --offscreen --exit --script "scripts/02_ligand_parade.py"

Or from inside ChimeraX:
    open scripts/02_ligand_parade.py
"""
from chimerax.core.commands import run

RECEPTOR_PDB = "receptor/6WGT_chainA_repaired.pdb"
OUT_MP4 = "output/02_ligand_parade.mp4"

# name, multi-seed mean affinity (kcal/mol), top-pose PDB -- ranked as in
# docking/results/multiseed_report.txt
LIGANDS = [
    ("LSD", -10.045, "poses/LSD_top_pose.pdb"),
    ("psilocybin", -7.837, "poses/psilocybin_top_pose.pdb"),
    ("serotonin", -7.083, "poses/serotonin_top_pose.pdb"),
    ("5-MeO-DMT", -6.982, "poses/5-MeO-DMT_top_pose.pdb"),
    ("psilocin", -6.937, "poses/psilocin_top_pose.pdb"),
    ("DMT", -6.872, "poses/DMT_top_pose.pdb"),
]

FRAMES_PER_LIGAND = 90  # 3s at 30fps


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
    run(session, f"view {receptor.atomspec}")

    run(session, '2dlabels create title text "Comparative screen: 5-HT2A, PDB 6WGT" '
                 'xpos 0.04 ypos 0.92 size 26 color white')
    run(session, '2dlabels create ligand text "" xpos 0.04 ypos 0.10 size 22 color white')
    run(session, '2dlabels create score text "" xpos 0.04 ypos 0.06 size 18 color orange')
    run(session, '2dlabels create caveat text "provisional -- redock gate currently FAILS, see docs/setup.md" '
                 'xpos 0.04 ypos 0.02 size 14 color gray')

    run(session, "movie record")

    for name, affinity, pdb_path in LIGANDS:
        ligand = open_one(session, pdb_path)
        run(session, f"style {ligand.atomspec} stick")
        run(session, f"color {ligand.atomspec} byhetero")
        run(session, f'2dlabels change ligand text "{name}"')
        run(session, f'2dlabels change score text "{affinity:.2f} kcal/mol (5-seed mean)"')
        run(session, f"view {receptor.atomspec},{ligand.atomspec}")
        run(session, f"turn y 1 {FRAMES_PER_LIGAND}")
        run(session, f"wait {FRAMES_PER_LIGAND}")
        run(session, f"close {ligand.atomspec}")

    run(session, "movie stop")
    run(session, f"movie encode {OUT_MP4} framerate 30 quality high")

    run(session, "2dlabels delete all")
    print(f"Wrote {OUT_MP4}")


main(session)  # noqa: F821 -- `session` is injected by ChimeraX when running a script
