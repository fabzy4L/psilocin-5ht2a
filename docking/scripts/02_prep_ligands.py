"""
Convert the ligand SMILES set to 3D structures, protonate at physiological
pH, and export PDBQT files for AutoDock Vina.

Requires: rdkit, meeko

Usage:
    python 02_prep_ligands.py --smi ../../data/ligands/ligands.smi \
        --out-dir ../../data/ligands/pdbqt
"""
import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy


def prep_one(smiles: str, name: str, out_dir: Path) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[skip] could not parse SMILES for {name}")
        return
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)

    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(setups[0])

    out_path = out_dir / f"{name}.pdbqt"
    out_path.write_text(pdbqt_string)
    print(f"Wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smi", default="../../data/ligands/ligands.smi")
    parser.add_argument("--out-dir", default="../../data/ligands/pdbqt")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.smi) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smiles, name = line.split("\t")
            prep_one(smiles, name, out_dir)


if __name__ == "__main__":
    main()
