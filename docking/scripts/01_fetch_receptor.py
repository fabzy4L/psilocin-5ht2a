"""
Fetch and lightly clean the 5-HT2A receptor structure for docking.

PDB 6WGT: cryo-EM structure of 5-HT2A bound to 25-CN-NBOH (agonist-bound,
active-state conformation) — a reasonable starting point for agonist docking.

Usage:
    python 01_fetch_receptor.py --pdb-id 6WGT --out ../../data/raw/6WGT.pdb
"""
import argparse
import urllib.request
from pathlib import Path


def fetch_pdb(pdb_id: str, out_path: Path) -> None:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, out_path)
    print(f"Saved {pdb_id} -> {out_path}")


def strip_to_receptor_chain(in_path: Path, out_path: Path, chain: str = "A") -> None:
    """
    Very crude first pass: keep ATOM/HETATM records for the given chain only,
    drop the bound ligand/lipids/waters. Inspect the PDB header first and
    adjust `chain` and any HETATM exclusions (nanobody, fusion partners,
    cholesterol hemisuccinate, etc. are common in GPCR cryo-EM structures
    and need to be stripped manually before use).
    """
    keep = []
    with open(in_path) as f:
        for line in f:
            if line.startswith(("ATOM", "TER")) and line[21] == chain:
                keep.append(line)
    with open(out_path, "w") as f:
        f.writelines(keep)
        f.write("END\n")
    print(f"Wrote stripped receptor -> {out_path} ({len(keep)} lines)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", default="6WGT")
    parser.add_argument("--out", default="../../data/raw/6WGT.pdb")
    parser.add_argument("--chain", default="A")
    args = parser.parse_args()

    raw_path = Path(args.out)
    fetch_pdb(args.pdb_id, raw_path)

    processed_path = Path("../../data/processed") / f"{args.pdb_id}_chain{args.chain}.pdb"
    strip_to_receptor_chain(raw_path, processed_path, args.chain)
