"""
Docking box / protocol validation: extract the co-crystallized ligand from
6WGT itself and re-dock it into the prepared receptor, then check the
top-scoring redocked pose against the experimental pose by heavy-atom RMSD.

6WGT's HET record "7LD" is LSD (HETNAM: "LYSERGIC ACID DIETHYLAMIDE";
independently confirmed against RCSB's entry, which names it
"(8alpha)-N,N-diethyl-6-methyl-9,10-didehydroergoline-8-carboxamide" and
cites Kim et al., "Structure of a Hallucinogen-Activated Gq-Coupled
5-HT2A Serotonin Receptor," Cell 2020, doi:10.1016/j.cell.2020.08.024) —
correcting an earlier note in docs/setup.md that mislabeled it as
25-CN-NBOH. That means this validation is a genuine self-redock test:
LSD is both the co-crystallized agonist AND one of the six comparators in
the Phase 1 screen, so a pass here is direct evidence the box/protocol
that produced the Phase 1 numbers can reproduce a known real pose, not
just an internally consistent one.

Result as of the current receptor/box (see docs/setup.md for the full
writeup): FAIL, RMSD 5.2 A. Run with --strict to make that failure
block downstream scripts rather than just get logged.

Unlike sert-s438t-escitalopram/scripts/validation_redock.py, no
cross-structure alignment is needed here (SERT's script aligns two
different PDB entries, 5I71 and 5I6Z, into a shared frame first) — 6WGT
already has the receptor and the ligand in one coordinate frame, so the
reference pose is extracted directly.

The reference pose is converted to PDBQT with OpenBabel rather than
meeko/RDKit (02_prep_ligands.py's route): meeko's route re-embeds 3D
coordinates from SMILES, which would throw away the experimental pose
this script exists to check against. OpenBabel's `-h` (no `--gen3d`)
keeps the input heavy-atom coordinates and just adds/merges hydrogens.

Usage:
    python 05_validation_redock.py \
        --raw-pdb ../../data/raw/6WGT.pdb \
        --receptor ../../data/processed/6WGT_chainA_adfr.pdbqt \
        --out-dir ../results

Requires: numpy, openbabel (obabel on PATH)
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from vina import Vina

BOX_CENTER = (25.115, 40.909, 54.225)
BOX_SIZE = (20.0, 20.0, 20.0)
EXHAUSTIVENESS = 16
N_POSES = 10
RMSD_THRESHOLD = 2.0  # Angstrom, standard redocking pass/fail cutoff

LIGAND_RESNAME = "7LD"
LIGAND_CHAIN = "A"


def extract_reference_pose(raw_pdb: Path, out_pdb: Path) -> None:
    kept = []
    for line in raw_pdb.read_text().splitlines(keepends=True):
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        chain = line[21]
        if resname == LIGAND_RESNAME and chain == LIGAND_CHAIN:
            kept.append(line)
    if not kept:
        sys.exit(f"ERROR: no {LIGAND_RESNAME} chain {LIGAND_CHAIN} HETATM "
                  f"records found in {raw_pdb}")
    kept.append("END\n")
    out_pdb.write_text("".join(kept))
    print(f"  Extracted {len(kept) - 1} atoms -> {out_pdb}")


def get_heavy_atom_coords_by_name(pdb_or_pdbqt: Path) -> dict:
    """
    Coords keyed by atom name, heavy atoms only. A multi-MODEL PDBQT (Vina's
    n_poses output) has the same atom names repeated once per pose — without
    stopping at the first ENDMDL, later poses silently overwrite earlier ones
    in the dict, so this would grade the *last* (worst-ranked) pose instead of
    the top one. Files with no MODEL record (plain reference PDBs) are read
    in full as before.
    """
    coords = {}
    for line in pdb_or_pdbqt.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        name = line[12:16].strip()
        element = (line[76:78].strip() if len(line) >= 78 else "") or name[0]
        if element.upper().startswith("H"):
            continue
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        coords[name] = np.array([x, y, z])
    return coords


def compute_rmsd(ref_path: Path, mob_path: Path) -> float | None:
    ref = get_heavy_atom_coords_by_name(ref_path)
    mob = get_heavy_atom_coords_by_name(mob_path)
    common = sorted(set(ref) & set(mob))
    if not common:
        return None
    diff = np.array([ref[k] for k in common]) - np.array([mob[k] for k in common])
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()))


def prepare_ligand_pdbqt(ref_pdb: Path, out_pdbqt: Path) -> bool:
    if not shutil.which("obabel"):
        return False
    cmd = ["obabel", str(ref_pdb), "-O", str(out_pdbqt),
           "--partialcharge", "gasteiger", "-h"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if not out_pdbqt.exists():
        print(f"  OpenBabel conversion failed: {result.stderr.strip()}")
        return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-pdb", default="../../data/raw/6WGT.pdb")
    parser.add_argument("--receptor",
                         default="../../data/processed/6WGT_chainA_adfr.pdbqt")
    parser.add_argument("--out-dir", default="../results")
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Vina search seed (default 42) — pinned so the RMSD/verdict "
             "is reproducible run to run, matching 04_docking_multiseed.py's "
             "convention")
    parser.add_argument(
        "--strict", action="store_true",
        help="exit 1 if the redock gate fails, so this can block a CI/"
             "pipeline step (e.g. 03_run_vina.py / 04_docking_multiseed.py) "
             "from running against an unvalidated receptor/box")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_pose = out_dir / "7LD_crystal_reference.pdb"
    lig_pdbqt = out_dir / "7LD_redock_input.pdbqt"
    docked_pdbqt = out_dir / "7LD_redocked.pdbqt"
    report_path = out_dir / "validation_report.txt"

    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log("=" * 60)
    log("DOCKING VALIDATION - 7LD (LSD) SELF-REDOCK")
    log("=" * 60)

    log("\n[1] Extracting co-crystallized 7LD pose from raw 6WGT...")
    extract_reference_pose(Path(args.raw_pdb), ref_pose)

    log("\n[2] Box proximity check...")
    ref_coords = np.array(list(get_heavy_atom_coords_by_name(ref_pose).values()))
    centroid = ref_coords.mean(axis=0)
    half = np.array(BOX_SIZE) / 2
    dist = float(np.linalg.norm(centroid - np.array(BOX_CENTER)))
    inside = bool(np.all(np.abs(centroid - np.array(BOX_CENTER)) <= half))
    log(f"  Crystal pose centroid: {centroid.round(3).tolist()}")
    log(f"  Box centre:            {BOX_CENTER}")
    log(f"  Distance:               {dist:.3f} A")
    log(f"  Inside grid:            {'YES' if inside else 'NO'}")
    log("  (expected to trivially pass: the box center IS this ligand's "
        "own crystal centroid, per docs/setup.md — this just confirms "
        "that wiring is still correct, not that docking works)")

    log("\n[3] Preparing ligand PDBQT (OpenBabel, coordinates preserved)...")
    if not prepare_ligand_pdbqt(ref_pose, lig_pdbqt):
        log("  SKIPPED — OpenBabel not available or conversion failed.")
        log(f"\nReport: {report_path}")
        report_path.write_text("\n".join(lines))
        return

    log("\n[4] Redocking into the prepared receptor...")
    v = Vina(sf_name="vina", seed=args.seed)
    v.set_receptor(args.receptor)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    v.set_ligand_from_file(str(lig_pdbqt))
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    v.write_poses(str(docked_pdbqt), n_poses=N_POSES, overwrite=True)
    best_affinity = float(v.energies(n_poses=N_POSES)[0][0])
    log(f"  Best redocked affinity: {best_affinity:.3f} kcal/mol")

    log("\n[5] Heavy-atom RMSD vs crystal pose (top pose only)...")
    rmsd = compute_rmsd(ref_pose, docked_pdbqt)
    if rmsd is None:
        log("  Could not compute RMSD (no matching atom names).")
        verdict = "FAIL"
    else:
        verdict = "PASS" if rmsd < RMSD_THRESHOLD else "FAIL"
        log(f"  RMSD: {rmsd:.3f} A  [{verdict} (<{RMSD_THRESHOLD} A)]")

    log("\n" + "=" * 60)
    log("SUMMARY")
    log("=" * 60)
    log(f"  Redock RMSD < {RMSD_THRESHOLD} A: {verdict}"
        + (f"  ({rmsd:.3f} A)" if rmsd is not None else ""))
    log(f"  Best redock affinity:    {best_affinity:.3f} kcal/mol")
    log(f"  Reference pose:          {ref_pose}")
    log(f"  Redocked poses:          {docked_pdbqt}")
    log(f"\n  GATE: {verdict}")

    report_path.write_text("\n".join(lines))
    (out_dir / "gate_status.json").write_text(json.dumps({
        "gate": verdict,
        "rmsd_angstrom": rmsd,
        "rmsd_threshold_angstrom": RMSD_THRESHOLD,
        "redock_affinity_kcal_mol": best_affinity,
        "ligand": LIGAND_RESNAME,
        "receptor": args.receptor,
        "seed": args.seed,
    }, indent=2))
    print(f"\nFull report saved: {report_path}")

    if args.strict and verdict != "PASS":
        sys.exit(
            f"GATE FAILED (RMSD {rmsd:.2f} A >= {RMSD_THRESHOLD} A threshold): "
            f"refusing to treat {args.receptor} + this box as validated. "
            f"03_run_vina.py / 04_docking_multiseed.py output against this "
            f"receptor should be reported as provisional until this gate "
            f"passes (see docs/setup.md)."
        )


if __name__ == "__main__":
    main()
