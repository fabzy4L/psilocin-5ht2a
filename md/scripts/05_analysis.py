"""
Post-production analysis: ligand pose RMSD, the Asp155 (Asp3.32)
salt-bridge distance to the ligand's basic ring nitrogen, and a gate on
whether the ligand RMSD has actually plateaued over the trajectory tail
rather than still drifting when the run ended.

Atom-name-keyed geometry follows the same convention as
docking/scripts/05_validation_redock.py (heavy atoms only, matched by
name) -- see docs/REDOCK_BUG_HANDOFF.md for why "never silently read past
a frame/MODEL boundary" is enforced as a first-class, tested gate here
rather than left to each script to get right ad hoc.

This script uses MDAnalysis (already in envs/environment.yml) to load the
trajectory rather than hand-parsing GRO/XTC frames -- multi-frame binary
XTC parsing is exactly the kind of format MDAnalysis exists to get right,
and hand-rolling it would repeat the mistake this whole file exists to
avoid. The MDAnalysis-dependent code is isolated in run_analysis() so the
pure gate/geometry functions below it stay importable and testable
without MDAnalysis installed (see tests/test_05_analysis.py).

Usage:
    python 05_analysis.py --topology ../system_prep/md.tpr \
        --trajectory ../system_prep/md.xtc \
        --ligand-resname UNK --strict

Requires: MDAnalysis, numpy (both in envs/environment.yml).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict

ASP332_RESID_DEFAULT = 155        # see docs/setup.md, confirmed against 6WGT
RMSD_PLATEAU_TOLERANCE_A = 0.5    # A, std dev over the tail counts as "plateaued"


def is_plateaued(rmsd_series: np.ndarray, tail_frac: float = 0.5,
                  tolerance: float = RMSD_PLATEAU_TOLERANCE_A) -> tuple[bool, float]:
    """Given a 1D per-frame RMSD array, checks whether the std dev over
    the last tail_frac of frames is within tolerance -- i.e. the pose has
    settled rather than still drifting when the run ended. Pure function,
    unit-tested without a real trajectory."""
    n = len(rmsd_series)
    if n == 0:
        return False, float("nan")
    tail = rmsd_series[int(n * (1 - tail_frac)):]
    std = float(np.std(tail))
    return std <= tolerance, std


def salt_bridge_distance(ligand_basic_n: np.ndarray, asp_od1: np.ndarray,
                          asp_od2: np.ndarray) -> float:
    """Distance from the ligand's basic ring nitrogen to the nearer of
    Asp3.32's two carboxylate oxygens -- the same geometry checked ad hoc
    during the redock-bug investigation (docs/REDOCK_BUG_HANDOFF.md),
    generalized into a reusable, tested function here."""
    return float(min(
        np.linalg.norm(ligand_basic_n - asp_od1),
        np.linalg.norm(ligand_basic_n - asp_od2),
    ))


def run_analysis(topology: Path, trajectory: Path, ligand_resname: str,
                  asp_resid: int, ligand_amine_name: str):
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms

    u = mda.Universe(str(topology), str(trajectory))
    ligand = u.select_atoms(f"resname {ligand_resname} and not name H*")
    if len(ligand) == 0:
        sys.exit(f"ERROR: no atoms matched 'resname {ligand_resname}' in {topology}")

    r = rms.RMSD(ligand, ligand).run()
    rmsd_series = r.results.rmsd[:, 2]  # column 2 = RMSD in angstrom

    asp = u.select_atoms(f"resid {asp_resid} and name OD1 OD2")
    amine = u.select_atoms(f"resname {ligand_resname} and name {ligand_amine_name}")
    salt_bridge_series = []
    if len(asp) == 2 and len(amine) == 1:
        for _ts in u.trajectory:
            salt_bridge_series.append(salt_bridge_distance(
                amine.positions[0], asp.positions[0], asp.positions[1]))

    return rmsd_series, np.array(salt_bridge_series)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topology", required=True)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--ligand-resname", required=True,
                         help="residue name of the ligand in the GROMACS "
                              "topology (CHARMM-GUI/CGenFF assigns this "
                              "when the ligand is parameterized -- check "
                              "topol.top)")
    parser.add_argument("--ligand-amine-name", default="N2",
                         help="atom name of the ligand's basic ring "
                              "nitrogen, for the Asp3.32 salt-bridge check "
                              "-- default N2 matches 7LD/LSD's PDB "
                              "chemical-component naming used elsewhere in "
                              "this repo (docking/results/"
                              "7LD_crystal_reference.pdb); override for "
                              "other ligands")
    parser.add_argument("--asp-resid", type=int, default=ASP332_RESID_DEFAULT)
    parser.add_argument("--tail-frac", type=float, default=0.5,
                         help="fraction of the trajectory (from the end) "
                              "checked for RMSD stability")
    parser.add_argument("--plateau-tolerance", type=float,
                         default=RMSD_PLATEAU_TOLERANCE_A)
    parser.add_argument("--out-dir", default="../analysis")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gl = GateLogger()
    gl.log("=" * 60)
    gl.log("POST-PRODUCTION ANALYSIS GATE")
    gl.log("=" * 60)

    gl.log(f"\n[1] Loading {args.trajectory} ...")
    rmsd_series, salt_bridge_series = run_analysis(
        Path(args.topology), Path(args.trajectory), args.ligand_resname,
        args.asp_resid, args.ligand_amine_name)

    gl.log(f"\n[2] Ligand RMSD over trajectory ({len(rmsd_series)} frames)...")
    plateaued, std = is_plateaued(rmsd_series, args.tail_frac, args.plateau_tolerance)
    gl.log(f"  Tail std dev ({args.tail_frac:.0%} of frames): {std:.3f} A "
           f"(tolerance {args.plateau_tolerance} A) -> "
           f"{'PLATEAUED' if plateaued else 'STILL DRIFTING'}")

    if len(salt_bridge_series):
        gl.log(f"\n[3] Asp{args.asp_resid} (Asp3.32) salt-bridge distance...")
        gl.log(f"  Mean: {salt_bridge_series.mean():.2f} A, "
               f"range: {salt_bridge_series.min():.2f}"
               f"-{salt_bridge_series.max():.2f} A")
    else:
        gl.log(f"\n[3] Asp{args.asp_resid} salt-bridge: SKIPPED "
               f"(couldn't find both OD1/OD2 and the named amine atom)")

    verdict = "PASS" if plateaued else "FAIL"
    gl.log(f"\n  GATE: {verdict}")
    gl.write_report(out_dir / "analysis_report.txt")
    write_gate_status(out_dir / "gate_status.json", verdict, {
        "rmsd_tail_std_angstrom": std,
        "plateau_tolerance_angstrom": args.plateau_tolerance,
        "n_frames": len(rmsd_series),
        "salt_bridge_mean_angstrom": (
            float(salt_bridge_series.mean()) if len(salt_bridge_series) else None
        ),
    })
    np.savetxt(out_dir / "ligand_rmsd_per_frame.csv", rmsd_series, delimiter=",",
               header="rmsd_angstrom", comments="")

    enforce_strict(args.strict, verdict,
        f"GATE FAILED: ligand RMSD has not plateaued (tail std "
        f"{std:.2f} A > {args.plateau_tolerance} A) -- the trajectory may "
        f"be too short, or the pose may genuinely be unstable. Don't treat "
        f"contacts/salt-bridge numbers from this run as settled-state "
        f"values yet.")


if __name__ == "__main__":
    main()
