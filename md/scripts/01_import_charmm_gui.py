"""
Import and validate a CHARMM-GUI Membrane Builder output archive into
md/system_prep/.

CHARMM-GUI's membrane builder download is a directory (or a .tar.gz/.zip
of one) containing GROMACS-ready topology, coordinates, and staged .mdp
files. This script does NOT talk to charmm-gui.org -- that step (upload
the cleaned receptor PDB, place the docked ligand, pick a POPC bilayer,
download) is interactive and has to happen in a browser first (see
docs/setup.md, "System building with CHARMM-GUI"). This script only
validates and stages what comes back from that download, so a missing or
partial export fails loudly here instead of silently breaking
02_minimize.py three steps later.

Usage:
    python 01_import_charmm_gui.py --input path/to/charmm-gui-download --strict

Requires: nothing beyond stdlib.
"""
import argparse
import shutil
import sys
import tarfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict

# CHARMM-GUI's gromacs/ output directory uses these canonical names.
REQUIRED_TOPOLOGY = ["topol.top"]
REQUIRED_INDEX = ["index.ndx"]
# The initial coordinate file name varies by CHARMM-GUI version/step
# (step5_input.gro is typical for Membrane Builder); accept common variants
# rather than hardcoding one.
COORDINATE_CANDIDATES = ["step5_input.gro", "step5_input.pdb", "system.gro"]
# CHARMM-GUI bundles a staged sequence of equilibration + production .mdp
# files; require at least one exists rather than hardcoding exact step
# numbers/counts, since that varies with membrane complexity.
MDP_GLOB = "*.mdp"


def find_archive_root(path: Path) -> Path:
    """If path is a CHARMM-GUI tar.gz/zip, extract it and return the
    directory to validate; if it's already a directory, return it as-is."""
    if path.is_dir():
        return path
    if not path.exists():
        sys.exit(f"ERROR: {path} does not exist")
    name = path.name
    if name.endswith((".tar.gz", ".tgz")):
        extract_to = path.parent / name.split(".tar")[0].removesuffix(".tgz")
        with tarfile.open(path) as tf:
            tf.extractall(extract_to)
        return extract_to
    if path.suffix == ".zip":
        extract_to = path.parent / path.stem
        with zipfile.ZipFile(path) as zf:
            zf.extractall(extract_to)
        return extract_to
    sys.exit(f"ERROR: {path} is neither a directory nor a .tar.gz/.zip archive")


def locate(root: Path, names) -> Path | None:
    for name in names:
        hits = list(root.rglob(name))
        if hits:
            return hits[0]
    return None


def validate_charmm_gui_dir(root: Path) -> dict:
    """Returns a dict describing what was found/missing. Pure function, no
    side effects beyond reading the filesystem -- kept separate from
    main() so it's directly unit-testable against a synthetic fixture
    directory (see tests/test_01_import_charmm_gui.py) without needing a
    real CHARMM-GUI download."""
    topol = locate(root, REQUIRED_TOPOLOGY)
    index = locate(root, REQUIRED_INDEX)
    coords = locate(root, COORDINATE_CANDIDATES)
    mdp_files = sorted(root.rglob(MDP_GLOB))

    missing = []
    if topol is None:
        missing.append(f"topology ({'/'.join(REQUIRED_TOPOLOGY)})")
    if index is None:
        missing.append(f"index ({'/'.join(REQUIRED_INDEX)})")
    if coords is None:
        missing.append(f"initial coordinates ({'/'.join(COORDINATE_CANDIDATES)})")
    if not mdp_files:
        missing.append("no .mdp files found (need at least equilibration + production)")

    return {
        "topol": str(topol) if topol else None,
        "index": str(index) if index else None,
        "coordinates": str(coords) if coords else None,
        "mdp_files": [str(p) for p in mdp_files],
        "missing": missing,
    }


def stage_into(out_dir: Path, found: dict) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for key in ("topol", "index", "coordinates"):
        src = found.get(key)
        if src:
            shutil.copy2(src, out_dir / Path(src).name)
    for mdp in found["mdp_files"]:
        shutil.copy2(mdp, out_dir / Path(mdp).name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True,
                         help="CHARMM-GUI Membrane Builder download: a "
                              "directory, or a .tar.gz/.zip archive of one")
    parser.add_argument("--out-dir", default="../system_prep")
    parser.add_argument("--strict", action="store_true",
                         help="exit 1 if required files are missing, so "
                              "02_minimize.py can't run against an "
                              "incomplete/miscopied CHARMM-GUI export")
    args = parser.parse_args()

    gl = GateLogger()
    gl.log("=" * 60)
    gl.log("CHARMM-GUI IMPORT VALIDATION")
    gl.log("=" * 60)

    root = find_archive_root(Path(args.input))
    gl.log(f"\n[1] Validating {root} ...")
    found = validate_charmm_gui_dir(root)

    gl.log(f"  Topology:     {found['topol'] or 'MISSING'}")
    gl.log(f"  Index:        {found['index'] or 'MISSING'}")
    gl.log(f"  Coordinates:  {found['coordinates'] or 'MISSING'}")
    gl.log(f"  .mdp files:   {len(found['mdp_files'])} found")
    for mdp in found["mdp_files"]:
        gl.log(f"    - {mdp}")

    verdict = "PASS" if not found["missing"] else "FAIL"
    if found["missing"]:
        gl.log("\n  Missing:")
        for m in found["missing"]:
            gl.log(f"    - {m}")

    out_dir = Path(args.out_dir)
    if verdict == "PASS":
        gl.log(f"\n[2] Staging into {out_dir} ...")
        stage_into(out_dir, found)
        gl.log("  Done.")
    else:
        gl.log("\n[2] Skipped staging -- required files missing.")

    gl.log(f"\n  GATE: {verdict}")

    out_dir.mkdir(parents=True, exist_ok=True)
    gl.write_report(out_dir / "import_report.txt")
    write_gate_status(out_dir / "gate_status.json", verdict, {
        "topol": found["topol"],
        "index": found["index"],
        "coordinates": found["coordinates"],
        "n_mdp_files": len(found["mdp_files"]),
        "missing": found["missing"],
        "source": str(root),
    })

    enforce_strict(args.strict, verdict,
        f"GATE FAILED: {root} is missing required CHARMM-GUI output "
        f"({', '.join(found['missing'])}). Re-download from charmm-gui.org "
        f"Membrane Builder (docs/setup.md) before running 02_minimize.py.")


if __name__ == "__main__":
    main()
