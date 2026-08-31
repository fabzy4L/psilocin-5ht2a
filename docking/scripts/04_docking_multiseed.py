"""
Multi-seed docking for the Phase 1 comparator screen.

Vina's search is stochastic; a single-seed run (03_run_vina.py) gives one
best-pose score per ligand with no sense of how much that score would move
on a different search trajectory. This reruns every ligand across several
seeds and reports mean +/- SD, so the psilocybin > psilocin ordering (and
anything else close to the ~0.5 kcal/mol noise floor typical of Vina) can
be judged against actual variance instead of a single point estimate.

Ported from sert-s438t-escitalopram/scripts/docking_multiseed.py (same
seed set, same idea) but adapted to this repo's Vina Python-API style
(03_run_vina.py) rather than the CLI-binary + config-file style the SERT
scripts use. The Vina Python bindings fix `seed` at object construction
(no per-dock seed argument), so a fresh Vina instance is built per seed
and maps are recomputed per seed — but only once per seed, not once per
ligand, since the receptor/box don't change across ligands.

Usage:
    python 04_docking_multiseed.py \
        --receptor ../../data/processed/6WGT_chainA_adfr.pdbqt \
        --ligand-dir ../../data/ligands/pdbqt \
        --out-dir ../results
"""
import argparse
import json
from pathlib import Path

import numpy as np
from vina import Vina

# Same orthosteric box as 03_run_vina.py — centroid of the co-crystallized
# agonist 7LD in 6WGT chain A, residue A1201. See docs/setup.md: 7LD is
# LSD (HETNAM confirms "LYSERGIC ACID DIETHYLAMIDE"), not 25-CN-NBOH as an
# earlier note here mistakenly said.
BOX_CENTER = (25.115, 40.909, 54.225)
BOX_SIZE = (20.0, 20.0, 20.0)

EXHAUSTIVENESS = 16
N_POSES = 10
SEEDS = [42, 123, 456, 789, 1001]


def dock_all_ligands(v: Vina, ligand_dir: Path) -> dict:
    scores = {}
    for ligand_pdbqt in sorted(ligand_dir.glob("*.pdbqt")):
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        scores[ligand_pdbqt.stem] = float(energies[0][0])
    return scores


def run_seed(receptor_pdbqt: str, ligand_dir: Path, seed: int) -> dict:
    print(f"  seed={seed} ...")
    v = Vina(sf_name="vina", seed=seed)
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    scores = dock_all_ligands(v, ligand_dir)
    for name, score in sorted(scores.items(), key=lambda kv: kv[1]):
        print(f"    {name:15s} {score:.2f} kcal/mol")
    return scores


def summarize(per_seed: dict) -> dict:
    """per_seed: {seed: {ligand: score}} -> {ligand: {mean, std, n, scores}}"""
    ligands = sorted({name for scores in per_seed.values() for name in scores})
    summary = {}
    for name in ligands:
        vals = [scores[name] for scores in per_seed.values() if name in scores]
        mean = float(np.mean(vals))
        std = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
        summary[name] = {"mean_kcal_mol": mean, "std_kcal_mol": std,
                          "n_seeds": len(vals), "scores": vals}
    return summary


def write_report(summary: dict, out_dir: Path) -> str:
    ranked = sorted(summary.items(), key=lambda kv: kv[1]["mean_kcal_mol"])
    lines = []
    lines.append("Psilocin/5-HT2A — Multi-Seed Docking Results")
    lines.append(f"AutoDock Vina | exhaustiveness={EXHAUSTIVENESS} | seeds={SEEDS}")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"{'Ligand':<15} {'mean (kcal/mol)':>16} {'SD':>8} {'n':>4}")
    lines.append("-" * 45)
    for name, s in ranked:
        lines.append(f"{name:<15} {s['mean_kcal_mol']:>16.3f} "
                      f"{s['std_kcal_mol']:>8.3f} {s['n_seeds']:>4}")
    lines.append("")
    lines.append("-" * 70)
    lines.append("NOTE: this is still rigid single-structure docking. Multi-seed")
    lines.append("averaging bounds Vina's search-stochasticity noise; it does not")
    lines.append("address receptor flexibility, solvation, or scoring-function bias")
    lines.append("(e.g. psilocybin's phosphate group finding polar contacts a")
    lines.append("flexible or MD treatment likely wouldn't preserve). See")
    lines.append("05_validation_redock.py for a check on whether this box/protocol")
    lines.append("reproduces the known co-crystallized pose at all.")

    report_text = "\n".join(lines)
    (out_dir / "multiseed_report.txt").write_text(report_text)
    return report_text


def main() -> None:
    global SEEDS
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand-dir", default="../../data/ligands/pdbqt")
    parser.add_argument("--out-dir", default="../results")
    parser.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    args = parser.parse_args()

    SEEDS = args.seeds

    ligand_dir = Path(args.ligand_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    per_seed = {}
    for seed in SEEDS:
        per_seed[seed] = run_seed(args.receptor, ligand_dir, seed)

    summary = summarize(per_seed)
    (out_dir / "multiseed_summary.json").write_text(
        json.dumps({"seeds": SEEDS, "exhaustiveness": EXHAUSTIVENESS,
                     "per_ligand": summary}, indent=2))

    print("\n" + "=" * 60)
    report = write_report(summary, out_dir)
    print(report)
    print(f"\nSummary: {out_dir / 'multiseed_summary.json'}")
    print(f"Report:  {out_dir / 'multiseed_report.txt'}")


if __name__ == "__main__":
    main()
