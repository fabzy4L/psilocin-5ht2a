"""
Flexible-residue variant of 05_validation_redock.py: redock 7LD (LSD) into
the receptor with the 9 pocket-adjacent side chains (A:89, 135, 136, 146,
147, 348, 350, 351, 356) treated as rotatable, using the rigid/flex PDBQT
split produced by 01d_prep_flexreceptor.py. Same (bug-fixed, first-MODEL-
only) RMSD-vs-crystal comparison as the rigid run, so the two are directly
comparable.

Run in the `vina` micromamba env (has the vina python package + obabel),
not the ADFRsuite pythonsh used for the receptor split:

    ~/micromamba/envs/vina/bin/python3 docking/scripts/05b_validation_redock_flex.py

IMPORTANT -- box size: the rigid gate's 20x20x20 A box is sized to the
LIGAND only. Once the 9 pocket residues are flexible, most of their side
chain atoms sit 10-20 A from the box center even in the static crystal
pose (e.g. LYS146/LYS350 terminal atoms, LEU147, ASP356) -- well outside
that box, and outside the Vina grid maps computed over it. Docking with
the 20A box gave a nonsensical +4.058 kcal/mol affinity and 3.724 A RMSD
(see docking/results/gate_status_flex_smallbox.json) purely from flexible
atoms falling off the map, not from anything about the pose itself. Sized
to actually cover their reach (44x44x36 A below), the same run gives
-10.001 kcal/mol and 0.958 A RMSD -- PASS, and consistent with the rigid
gate's -10.073 kcal/mol / 0.780 A. Anyone reusing this script for a
different flexible-residue set needs to re-check this box against the new
residues' reach, not just assume the ligand-sized box is fine.
"""
import json
from pathlib import Path

import numpy as np
from vina import Vina

BOX_CENTER = (25.115, 40.909, 54.225)
BOX_SIZE = (44.0, 44.0, 36.0)  # sized to the 9 flex residues' reach -- see box-size note above, NOT just the ligand
EXHAUSTIVENESS = 16
N_POSES = 10
RMSD_THRESHOLD = 2.0
SEED = 42


def get_heavy_atom_coords_by_name(path: Path) -> dict:
    """Same fix as 05_validation_redock.py: stop at first ENDMDL so a
    multi-MODEL PDBQT is graded on pose 1, not the last pose in the file."""
    coords = {}
    for line in path.read_text().splitlines():
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


def compute_rmsd(ref_path: Path, mob_path: Path):
    ref = get_heavy_atom_coords_by_name(ref_path)
    mob = get_heavy_atom_coords_by_name(mob_path)
    common = sorted(set(ref) & set(mob))
    if not common:
        return None
    diff = np.array([ref[k] for k in common]) - np.array([mob[k] for k in common])
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()))


def main():
    out_dir = Path("docking/results")
    ref_pose = out_dir / "7LD_crystal_reference.pdb"
    lig_pdbqt = out_dir / "7LD_redock_input.pdbqt"
    docked_pdbqt = out_dir / "7LD_redocked_flex_bigbox.pdbqt"

    rigid = "data/processed/6WGT_chainA_adfr_rigid.pdbqt"
    flex = "data/processed/6WGT_chainA_adfr_flex.pdbqt"

    print("Redocking 7LD with 9 flexible pocket side chains...")
    v = Vina(sf_name="vina", seed=SEED)
    v.set_receptor(rigid_pdbqt_filename=rigid, flex_pdbqt_filename=flex)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    v.set_ligand_from_file(str(lig_pdbqt))
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
    v.write_poses(str(docked_pdbqt), n_poses=N_POSES, overwrite=True)
    best_affinity = float(v.energies(n_poses=N_POSES)[0][0])
    print(f"Best flexible-redock affinity: {best_affinity:.3f} kcal/mol")

    rmsd = compute_rmsd(ref_pose, docked_pdbqt)
    verdict = "PASS" if (rmsd is not None and rmsd < RMSD_THRESHOLD) else "FAIL"
    print(f"Ligand-only RMSD vs crystal pose (top pose): {rmsd:.3f} A [{verdict}]")

    (out_dir / "gate_status_flex_bigbox.json").write_text(json.dumps({
        "gate": verdict,
        "rmsd_angstrom": rmsd,
        "rmsd_threshold_angstrom": RMSD_THRESHOLD,
        "redock_affinity_kcal_mol": best_affinity,
        "ligand": "7LD",
        "receptor_rigid": rigid,
        "receptor_flex": flex,
        "flexible_residues": [89, 135, 136, 146, 147, 348, 350, 351, 356],
        "seed": SEED,
        "note": "flexible-residue redock, for comparison against the rigid "
                "gate result (0.780 A). Ligand-only RMSD -- does not include "
                "flexible side-chain atom displacement.",
    }, indent=2))
    print("wrote docking/results/gate_status_flex_bigbox.json")


if __name__ == "__main__":
    main()
