"""
Convert the cleaned receptor PDB to PDBQT for AutoDock Vina, via Meeko's
receptor prep CLI (`mk_prepare_receptor.py`, installed alongside the
`meeko` package already listed in envs/environment.yml).

Run AFTER 01_fetch_receptor.py, BEFORE 03_run_vina.py.

6WGT is a ~3-4 A cryo-EM structure; many solvent-exposed side chains are
only partially resolved (deposited truncated at the last atom with clear
density) and will fail meeko's residue template matching. This script
auto-deletes such residues if they're safely far from the docking box
(--bad-res-radius, default 5 A) and raises an error otherwise, since a
truncated side chain near the pocket would need manual repair, not
silent deletion.

Usage:
    python 01b_prep_receptor.py \
        --in-pdb ../../data/processed/6WGT_chainA.pdb \
        --out-prefix ../../data/processed/6WGT_chainA \
        --box-center 25.115 40.909 54.225 \
        --box-size 20 20 20
"""
import argparse
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-pdb", required=True)
    parser.add_argument("--out-prefix", required=True,
                         help="output path without extension; "
                              "writes <prefix>.pdbqt plus box files")
    parser.add_argument("--box-center", nargs=3, type=float, required=True,
                         metavar=("X", "Y", "Z"))
    parser.add_argument("--box-size", nargs=3, type=float,
                         default=(20.0, 20.0, 20.0), metavar=("X", "Y", "Z"))
    parser.add_argument(
        "--bad-res-radius", type=float, default=5.0,
        help="6WGT is a moderate-resolution cryo-EM structure with many "
             "solvent-exposed side chains truncated at the Cbeta (no "
             "resolved density past that point) — these fail meeko's "
             "residue template matching. Auto-delete such residues if "
             "their center-of-mass is more than this many Angstroms "
             "outside the docking box (safe: far from the pocket, doesn't "
             "affect Vina scoring); anything closer raises an error and "
             "needs manual side-chain repair (e.g. PDBFixer/SCWRL4) "
             "instead, since it could sit in or near the binding site.")
    args = parser.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mk_prepare_receptor",
        "--read_pdb", args.in_pdb,
        "-o", str(out_prefix),
        "--box_center", *map(str, args.box_center),
        "--box_size", *map(str, args.box_size),
        "--delete_bad_res_from_box_radius", str(args.bad_res_radius),
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"Wrote receptor PDBQT + box files -> {out_prefix}.*")


if __name__ == "__main__":
    main()
