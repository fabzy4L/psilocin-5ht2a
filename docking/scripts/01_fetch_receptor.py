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


BRIL_RESNUM_RANGE = (1001, 1106)  # inclusive; see NOTE below


def strip_to_receptor_chain(
    in_path: Path, out_path: Path, chain: str = "A", drop_bril: bool = True
) -> None:
    """
    Very crude first pass: keep ATOM records for the given chain only, drop
    all HETATM (bound ligand 7LD, cholesterol CLR, oleic acid OLA, PEG/1PE
    crystallization additives, phosphate PO4 — confirmed via `HET` records
    in 6WGT).

    NOTE on the BRIL fusion: 6WGT's construct is a single-chain
    5-HT2A/cytochrome-b562RIL (BRIL) fusion, not a separate chain. BRIL
    replaces native ICL3 in-frame between receptor residues 265 and 311,
    and is itself renumbered 1001-1106 in the PDB (confirmed by walking
    CA residue numbers in chain A: the sequence jumps 265->1001 and then
    1106->311). Chain-only filtering does NOT remove it, so this function
    drops BRIL_RESNUM_RANGE by default (drop_bril=True) — it is not part
    of the biological receptor (docs/setup.md), and leaving it in causes
    real downstream problems: BRIL folds back spatially close to parts of
    the TM bundle despite being sequence-distant, which trips up
    distance-based bond perception in receptor-prep tools (meeko) with
    spurious inter-residue bonds. Pass drop_bril=False to keep it (e.g.
    to inspect the raw construct).

    Also inserts a TER record at every remaining residue-numbering
    discontinuity (missing loop density: 182-183, 216-218), not just
    between chains. Without this, downstream tools try to bond across the
    gap and fail residue template matching for everything after it.
    """
    keep = []
    prev_resnum = None
    with open(in_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "TER")) or line[21] != chain:
                continue
            if line.startswith("ATOM"):
                resnum = int(line[22:26])
                if drop_bril and BRIL_RESNUM_RANGE[0] <= resnum <= BRIL_RESNUM_RANGE[1]:
                    continue
                if prev_resnum is not None and resnum != prev_resnum and resnum != prev_resnum + 1:
                    keep.append("TER\n")
                prev_resnum = resnum
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
