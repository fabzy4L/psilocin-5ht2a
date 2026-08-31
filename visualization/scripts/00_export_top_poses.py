"""
Regenerate everything the ChimeraX visualization scripts in this
directory need, so `visualization/` is self-contained and portable to a
machine that only has this repo's git history, not its WSL/conda docking
environment:

1. Extract the top (mode 1) pose from each multi-model Vina PDBQT output
   and convert it to a plain PDB via OpenBabel (already a project
   dependency). `-f 1 -l 1` takes only the first model (the top-ranked
   pose) rather than all 9-10.
2. Copy the receptor and crystal-reference PDBs used by the redock gate
   into visualization/receptor/ — data/processed/* is gitignored (only
   .gitkeep is tracked there), so the ChimeraX scripts can't reach it
   directly on a fresh clone; docking/results/7LD_crystal_reference.pdb
   *is* tracked, but is duplicated here too so both scripts only ever
   read from within visualization/.

Usage:
    python 00_export_top_poses.py
"""
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "docking" / "results"
POSES_DIR = ROOT / "visualization" / "poses"
RECEPTOR_DIR = ROOT / "visualization" / "receptor"

POSE_SOURCES = {
    "LSD": RESULTS / "LSD_docked.pdbqt",
    "psilocybin": RESULTS / "psilocybin_docked.pdbqt",
    "serotonin": RESULTS / "serotonin_docked.pdbqt",
    "psilocin": RESULTS / "psilocin_docked.pdbqt",
    "5-MeO-DMT": RESULTS / "5-MeO-DMT_docked.pdbqt",
    "DMT": RESULTS / "DMT_docked.pdbqt",
    "7LD_redocked": RESULTS / "7LD_redocked.pdbqt",
}

RECEPTOR_SOURCES = {
    "6WGT_chainA_repaired.pdb": ROOT / "data" / "processed" / "6WGT_chainA_repaired.pdb",
    "7LD_crystal_reference.pdb": RESULTS / "7LD_crystal_reference.pdb",
}


def main() -> None:
    POSES_DIR.mkdir(parents=True, exist_ok=True)
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)

    missing = [name for name, path in POSE_SOURCES.items() if not path.exists()]
    missing += [name for name, path in RECEPTOR_SOURCES.items() if not path.exists()]
    if missing:
        raise SystemExit(
            f"Missing source file(s), run the Phase 1 pipeline first "
            f"(see docs/setup.md): {missing}"
        )

    for name, src in POSE_SOURCES.items():
        out_pdb = POSES_DIR / f"{name}_top_pose.pdb"
        cmd = ["obabel", str(src), "-O", str(out_pdb), "-f", "1", "-l", "1"]
        subprocess.run(cmd, check=True)
        print(f"Wrote {out_pdb}")

    for name, src in RECEPTOR_SOURCES.items():
        out_path = RECEPTOR_DIR / name
        shutil.copyfile(src, out_path)
        print(f"Copied {out_path}")

    print(f"\n{len(POSE_SOURCES)} top-pose PDBs in {POSES_DIR}, "
          f"{len(RECEPTOR_SOURCES)} receptor file(s) in {RECEPTOR_DIR}")


if __name__ == "__main__":
    main()
