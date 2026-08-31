"""
Repair missing side-chain atoms in the stripped receptor PDB.

6WGT is a ~3-4 A cryo-EM structure; many solvent-exposed side chains are
only partially resolved and deposited truncated (e.g. Ile/Leu missing
their delta carbons, Lys missing everything past Cbeta). This adds the
missing heavy atoms with PDBFixer (standard geometry, not rotamer
optimization) WITHOUT filling missing *residues* (real loop gaps like
182-183, 216-218 are left alone — that needs actual loop modeling, out
of scope here).

PDBFixer also treats every residue-numbering discontinuity as a chain
terminus and adds a spurious OXT (terminal carboxylate O) there. Since
the only genuine terminus in this construct is the last resolved residue
(399 — itself just where the crystallization construct was truncated,
not the receptor's real C-terminus), this script strips every other OXT
it introduces.

Run AFTER 01_fetch_receptor.py, BEFORE receptor->PDBQT conversion
(01b_prep_receptor.py or 01c_prep_receptor_adfr.py).

Requires: pdbfixer, openmm (both pip-installable; not yet in
envs/environment.yml — add them, or install ADFRsuite instead if going
the 01c route exclusively, since MGLTools' own repair pass has looser
requirements).

Usage:
    python 01a_repair_sidechains.py \
        --in-pdb ../../data/processed/6WGT_chainA.pdb \
        --out-pdb ../../data/processed/6WGT_chainA_repaired.pdb \
        --true-terminus-resnum 399
"""
import argparse

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-pdb", required=True)
    parser.add_argument("--out-pdb", required=True)
    parser.add_argument(
        "--true-terminus-resnum", type=int, required=True,
        help="residue number of the construct's real C-terminus; OXT is "
             "kept there and stripped everywhere else")
    args = parser.parse_args()

    fixer = PDBFixer(filename=args.in_pdb)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # don't fill real loop gaps, atoms only
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    print(f"Residues with missing atoms: {len(fixer.missingAtoms)}")
    fixer.addMissingAtoms()

    tmp_out = args.out_pdb + ".tmp"
    with open(tmp_out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    removed = 0
    with open(tmp_out) as f:
        lines = f.readlines()
    kept = []
    for line in lines:
        if (line.startswith("ATOM") and line[12:16].strip() == "OXT"
                and int(line[22:26]) != args.true_terminus_resnum):
            removed += 1
            continue
        kept.append(line)
    with open(args.out_pdb, "w") as f:
        f.writelines(kept)

    import os
    os.remove(tmp_out)
    print(f"Removed {removed} spurious OXT atom(s)")
    print(f"Wrote repaired receptor -> {args.out_pdb}")


if __name__ == "__main__":
    main()
