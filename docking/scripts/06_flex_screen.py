"""
Flexible-residue comparator screen: redock the 5 remaining Phase 1c
ligands (psilocybin, psilocin, serotonin, 5-MeO-DMT, DMT -- LSD is already
validated in 05b_validation_redock_flex.py) with the same 9 pocket-adjacent
side chains flexible, using the corrected 44x44x36 A box (see docs/setup.md,
"Flexible-residue docking: tooling gap closed" -- the rigid gate's
20x20x20 A box is sized to the ligand only and cuts off 8 of the 9 flex
residues' actual reach).

This is Aim 1's real test (README Phase 1d): does the psilocybin >
psilocin ordering from the rigid blind screen (03_run_vina.py) survive
once the receptor can flex around the pocket, or was it a rigid-receptor
scoring artifact?

Requires the rigid/flex receptor split from 01d_prep_flexreceptor.py.
That split is a gitignored build artifact (not committed) -- regenerate
it first if data/processed/6WGT_chainA_adfr_{rigid,flex}.pdbqt don't
exist:

    ~/ADFRsuite-1.1dev/bin/pythonsh docking/scripts/01d_prep_flexreceptor.py \
        --receptor data/processed/6WGT_chainA_adfr.pdbqt \
        --residues 89,135,136,146,147,348,350,351,356 \
        --out-dir data/processed

Usage (run in the vina micromamba env, which has the vina python package):
    python 06_flex_screen.py \
        --rigid-receptor ../../data/processed/6WGT_chainA_adfr_rigid.pdbqt \
        --flex-receptor ../../data/processed/6WGT_chainA_adfr_flex.pdbqt \
        --ligand-dir ../../data/ligands/pdbqt \
        --out-dir ../results
"""
import argparse
import json
from pathlib import Path

from vina import Vina

BOX_CENTER = (25.115, 40.909, 54.225)
BOX_SIZE = (44.0, 44.0, 36.0)  # sized to the 9 flex residues' reach, not just the ligand
EXHAUSTIVENESS = 16
N_POSES = 10
SEED = 42  # matches 05b_validation_redock_flex.py's convention

SKIP_LIGANDS = {"LSD"}  # already validated in 05b_validation_redock_flex.py (0.958 A PASS)


def dock_ligand(v: Vina, ligand_pdbqt: Path, out_dir: Path) -> dict:
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)

    name = ligand_pdbqt.stem
    out_pdbqt = out_dir / f"{name}_flex_docked.pdbqt"
    v.write_poses(str(out_pdbqt), n_poses=N_POSES, overwrite=True)

    energies = v.energies(n_poses=N_POSES)
    return {"name": name, "best_score_kcal_mol": float(energies[0][0])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rigid-receptor",
                         default="../../data/processed/6WGT_chainA_adfr_rigid.pdbqt")
    parser.add_argument("--flex-receptor",
                         default="../../data/processed/6WGT_chainA_adfr_flex.pdbqt")
    parser.add_argument("--ligand-dir", default="../../data/ligands/pdbqt")
    parser.add_argument("--out-dir", default="../results")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    v = Vina(sf_name="vina", seed=SEED)
    v.set_receptor(rigid_pdbqt_filename=args.rigid_receptor,
                    flex_pdbqt_filename=args.flex_receptor)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)

    results = []
    for ligand_pdbqt in sorted(Path(args.ligand_dir).glob("*.pdbqt")):
        if ligand_pdbqt.stem in SKIP_LIGANDS:
            print(f"Skipping {ligand_pdbqt.stem} (already validated in "
                  f"05b_validation_redock_flex.py)")
            continue
        print(f"Flexibly docking {ligand_pdbqt.stem} ...")
        results.append(dock_ligand(v, ligand_pdbqt, out_dir))

    results.sort(key=lambda r: r["best_score_kcal_mol"])
    summary_path = out_dir / "flex_screen_summary.json"
    summary_path.write_text(json.dumps(results, indent=2))
    print(f"\nSummary written to {summary_path}")
    for r in results:
        print(f"  {r['name']:15s} {r['best_score_kcal_mol']:.2f} kcal/mol")


if __name__ == "__main__":
    main()
