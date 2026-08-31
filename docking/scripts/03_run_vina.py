"""
Run AutoDock Vina for each prepared ligand against the prepared receptor.

Box center/size should be set to the orthosteric binding pocket of 5-HT2A —
find this by inspecting the co-crystallized ligand position in the raw PDB
(6WGT's bound agonist coordinates) before running. Placeholder values below
must be replaced.

Requires: receptor PDBQT (prepare with e.g. ADFR/prepare_receptor or
OpenBabel from the cleaned chain PDB — not scripted here since it typically
needs manual inspection for missing loops/waters).

Usage:
    python 03_run_vina.py --receptor ../../data/processed/6WGT_receptor.pdbqt \
        --ligand-dir ../../data/ligands/pdbqt \
        --out-dir ../results
"""
import argparse
import json
from pathlib import Path

from vina import Vina

# TODO: replace with the actual orthosteric pocket center from 6WGT
BOX_CENTER = (0.0, 0.0, 0.0)
BOX_SIZE = (20.0, 20.0, 20.0)


def dock_ligand(v: Vina, ligand_pdbqt: Path, out_dir: Path) -> dict:
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.dock(exhaustiveness=16, n_poses=10)

    name = ligand_pdbqt.stem
    out_pdbqt = out_dir / f"{name}_docked.pdbqt"
    v.write_poses(str(out_pdbqt), n_poses=10, overwrite=True)

    energies = v.energies(n_poses=10)
    return {"name": name, "best_score_kcal_mol": float(energies[0][0])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand-dir", default="../../data/ligands/pdbqt")
    parser.add_argument("--out-dir", default="../results")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    v = Vina(sf_name="vina")
    v.set_receptor(args.receptor)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)

    results = []
    for ligand_pdbqt in sorted(Path(args.ligand_dir).glob("*.pdbqt")):
        print(f"Docking {ligand_pdbqt.stem} ...")
        results.append(dock_ligand(v, ligand_pdbqt, out_dir))

    results.sort(key=lambda r: r["best_score_kcal_mol"])
    summary_path = out_dir / "summary.json"
    summary_path.write_text(json.dumps(results, indent=2))
    print(f"\nSummary written to {summary_path}")
    for r in results:
        print(f"  {r['name']:15s} {r['best_score_kcal_mol']:.2f} kcal/mol")


if __name__ == "__main__":
    main()
