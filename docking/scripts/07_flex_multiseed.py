"""
Multi-seed confirmation for the flexible-residue comparator screen
(06_flex_screen.py) -- the caveat flagged in that script's report: it was
single-seed (42), so a close score gap (like psilocybin vs psilocin, or
5-MeO-DMT vs psilocin) can't yet be told apart from search-stochasticity
noise. This reruns every ligand across the same seed set
04_docking_multiseed.py used for the rigid screen ([42, 123, 456, 789,
1001]), with the flexible receptor split and the 44x44x36 A box
(docs/setup.md), and reports mean +/- SD per ligand the same way.

Includes all 6 ligands (unlike 06_flex_screen.py, which skips LSD since
that script's purpose was the LSD pose-accuracy validation already done
in 05b_validation_redock_flex.py at seed 42). This script only concerns
affinity/ranking stability across seeds, not pose accuracy, so including
LSD here gives a complete, directly-comparable table against
04_docking_multiseed.py's rigid results.

Requires the rigid/flex receptor split from 01d_prep_flexreceptor.py
(gitignored build artifact, regenerate if missing -- see
06_flex_screen.py's docstring).

Usage (run in the vina micromamba env):
    python 07_flex_multiseed.py \
        --rigid-receptor ../../data/processed/6WGT_chainA_adfr_rigid.pdbqt \
        --flex-receptor ../../data/processed/6WGT_chainA_adfr_flex.pdbqt \
        --ligand-dir ../../data/ligands/pdbqt \
        --out-dir ../results
"""
import argparse
import json
from pathlib import Path

import numpy as np
from vina import Vina

BOX_CENTER = (25.115, 40.909, 54.225)
BOX_SIZE = (44.0, 44.0, 36.0)  # sized to the 9 flex residues' reach, not just the ligand
EXHAUSTIVENESS = 16
N_POSES = 10
SEEDS = [42, 123, 456, 789, 1001]  # same set as 04_docking_multiseed.py


def dock_all_ligands(v: Vina, ligand_dir: Path) -> dict:
    scores = {}
    for ligand_pdbqt in sorted(ligand_dir.glob("*.pdbqt")):
        v.set_ligand_from_file(str(ligand_pdbqt))
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        scores[ligand_pdbqt.stem] = float(energies[0][0])
    return scores


def run_seed(rigid_receptor: str, flex_receptor: str, ligand_dir: Path, seed: int) -> dict:
    print(f"  seed={seed} ...")
    v = Vina(sf_name="vina", seed=seed)
    v.set_receptor(rigid_pdbqt_filename=rigid_receptor, flex_pdbqt_filename=flex_receptor)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    scores = dock_all_ligands(v, ligand_dir)
    for name, score in sorted(scores.items(), key=lambda kv: kv[1]):
        print(f"    {name:15s} {score:.2f} kcal/mol")
    return scores


def summarize(per_seed: dict) -> dict:
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
    lines.append("Psilocin/5-HT2A -- Flexible-Residue Multi-Seed Confirmation")
    lines.append(f"AutoDock Vina | exhaustiveness={EXHAUSTIVENESS} | seeds={SEEDS} | "
                  f"box={BOX_SIZE} A (flex-residue-sized, not ligand-sized)")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"{'Ligand':<15} {'mean (kcal/mol)':>16} {'SD':>8} {'n':>4}")
    lines.append("-" * 45)
    for name, s in ranked:
        lines.append(f"{name:<15} {s['mean_kcal_mol']:>16.3f} "
                      f"{s['std_kcal_mol']:>8.3f} {s['n_seeds']:>4}")
    lines.append("")
    lines.append("-" * 70)
    lines.append("Compare against docking/results/multiseed_report.txt (rigid).")
    lines.append("This confirms whether the flexible screen's single-seed result")
    lines.append("(docking/results/flex_screen_report.txt) holds up against search")
    lines.append("stochasticity the same way the rigid multi-seed run confirmed the")
    lines.append("original blind screen -- not a claim about pose accuracy or")
    lines.append("solvation effects, which flexible-side-chain Vina still doesn't")
    lines.append("model (see docs/setup.md).")

    report_text = "\n".join(lines)
    (out_dir / "flex_multiseed_report.txt").write_text(report_text)
    return report_text


def main() -> None:
    global SEEDS
    parser = argparse.ArgumentParser()
    parser.add_argument("--rigid-receptor",
                         default="../../data/processed/6WGT_chainA_adfr_rigid.pdbqt")
    parser.add_argument("--flex-receptor",
                         default="../../data/processed/6WGT_chainA_adfr_flex.pdbqt")
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
        per_seed[seed] = run_seed(args.rigid_receptor, args.flex_receptor, ligand_dir, seed)

    summary = summarize(per_seed)
    (out_dir / "flex_multiseed_summary.json").write_text(
        json.dumps({"seeds": SEEDS, "exhaustiveness": EXHAUSTIVENESS,
                     "box_size_angstrom": BOX_SIZE, "per_ligand": summary}, indent=2))

    print("\n" + "=" * 60)
    report = write_report(summary, out_dir)
    print(report)
    print(f"\nSummary: {out_dir / 'flex_multiseed_summary.json'}")
    print(f"Report:  {out_dir / 'flex_multiseed_report.txt'}")


if __name__ == "__main__":
    main()
