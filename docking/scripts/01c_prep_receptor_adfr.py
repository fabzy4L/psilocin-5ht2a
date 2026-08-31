"""
Convert the repaired receptor PDB to PDBQT via ADFRsuite's
prepare_receptor (a thin wrapper around the classic AutoDockTools
prepare_receptor4.py). This is the path that actually worked end-to-end
on 6WGT — meeko's mk_prepare_receptor (01b_prep_receptor.py) does strict
RDKit-template residue matching and chokes on this structure's terminal
geometry; prepare_receptor4.py just builds bonds/adds hydrogens
geometrically and is far more forgiving of exactly this kind of
cryo-EM messiness.

Run AFTER 01a_repair_sidechains.py.

ADFRsuite is not on PyPI/conda-forge and its official Linux installer is
a legacy Tcl/Tk (InstallJammer) binary that doesn't support headless/CLI
installation cleanly. The reliable path is the community conda-forge
channel build, installed into its own env (needs Python 2-era deps that
would fight with the rest of envs/environment.yml, hence a separate env
rather than folding it in):

    conda create -n adfr -c hcc -c conda-forge adfr-suite
    conda activate adfr
    # or on Windows: run this from WSL2 (Ubuntu) - ADFRsuite has no
    # native Windows build.

Usage:
    conda run -n adfr python 01c_prep_receptor_adfr.py \
        --in-pdb ../../data/processed/6WGT_chainA_repaired.pdb \
        --out-pdbqt ../../data/processed/6WGT_chainA_adfr.pdbqt
"""
import argparse
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-pdb", required=True)
    parser.add_argument("--out-pdbqt", required=True)
    args = parser.parse_args()

    out_path = Path(args.out_pdbqt)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "prepare_receptor",
        "-r", args.in_pdb,
        "-o", str(out_path),
        "-A", "bonds_hydrogens",
        "-v",
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"Wrote receptor PDBQT -> {out_path}")


if __name__ == "__main__":
    main()
