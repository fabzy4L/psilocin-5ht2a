"""
Split the ADFR-prepared receptor into rigid + flexible PDBQT for the 9
pocket-adjacent side chains (A:89, 135, 136, 146, 147, 348, 350, 351, 356 --
see docs/setup.md). Flexible side chains, not backbone (allowed_bonds
defaults to 'backbone' meaning backbone bonds stay rigid; side-chain
torsions become rotatable).

docs/setup.md documents this as blocked: the `hcc` conda ADFRsuite build
has no `prepare_flexreceptor` CLI, and meeko's --flexres mode fails on
these exact residues. Both are true for their respective binaries, but the
underlying library class (AD4FlexibleReceptorPreparation) is bundled
inside this machine's local ADFRsuite-1.1dev install
(~/ADFRsuite-1.1dev/CCSBpckgs/AutoDockTools/MoleculePreparation.py) and
imports fine -- confirmed by hand. This script calls that class directly,
bypassing the missing CLI wrapper. Must run under ADFRsuite's bundled
Python 2 interpreter (pythonsh), not the repo's normal python3 env:

    ~/ADFRsuite-1.1dev/bin/pythonsh docking/scripts/01d_prep_flexreceptor.py \
        --receptor data/processed/6WGT_chainA_adfr.pdbqt \
        --residues 89,135,136,146,147,348,350,351,356 \
        --out-dir data/processed
"""
import getopt
import os
import sys

from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4FlexibleReceptorPreparation


def main():
    opts, _ = getopt.getopt(sys.argv[1:], "", ["receptor=", "residues=", "out-dir="])
    opt = dict(opts)
    receptor_filename = opt["--receptor"]
    residue_numbers = [s.strip() for s in opt["--residues"].split(",")]
    out_dir = opt.get("--out-dir", "data/processed")

    rec = Read(receptor_filename)[0]
    rec.buildBondsByDistance()

    all_res = []
    for num in residue_numbers:
        matched = [r for r in rec.chains.residues if r.number == num]
        if not matched:
            print "WARNING: no residue numbered", num, "found"
            continue
        all_res.extend(matched)
    print "matched", len(all_res), "of", len(residue_numbers), "requested residues:"
    for r in all_res:
        print "  ", r.full_name()

    base = os.path.splitext(os.path.basename(receptor_filename))[0]
    rigid_out = os.path.join(out_dir, base + "_rigid.pdbqt")
    flex_out = os.path.join(out_dir, base + "_flex.pdbqt")

    AD4FlexibleReceptorPreparation(
        rec,
        mode="automatic",
        rigid_filename=rigid_out,
        flexres_filename=flex_out,
        residues=all_res,
    )
    print "wrote", rigid_out
    print "wrote", flex_out


if __name__ == "__main__":
    main()
