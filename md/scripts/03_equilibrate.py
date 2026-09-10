"""
NVT/NPT equilibration gate: after each phase, extract the relevant
observable via `gmx energy` (temperature for NVT, pressure + density for
NPT) and check it's actually stable over the tail of the run, not just
"the run finished" -- a run that finishes without exception can still be
mid-drift if given too few equilibration steps.

Usage:
    python 03_equilibrate.py --phase nvt --gro ../system_prep/em.gro \
        --top ../system_prep/topol.top --mdp ../system_prep/nvt.mdp \
        --out-dir ../system_prep --strict

    python 03_equilibrate.py --phase npt --gro ../system_prep/nvt.gro \
        --top ../system_prep/topol.top --mdp ../system_prep/npt.mdp \
        --out-dir ../system_prep --strict

Requires: gmx on PATH.
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict

PHASE_OBSERVABLES = {
    "nvt": ["Temperature"],
    "npt": ["Pressure", "Density"],
}
DEFAULT_TOLERANCE = {
    "Temperature": 2.0,   # K
    "Pressure": 50.0,     # bar -- pressure is inherently noisy at these timescales
    "Density": 5.0,       # kg/m^3
}


def parse_xvg(path: Path) -> list[tuple[float, float]]:
    """Returns [(time_ps, value), ...] from a gmx energy .xvg, skipping
    xmgrace header/comment lines ('#'/'@'). Pure function -- unit-tested
    against a synthetic .xvg fixture in tests/test_03_equilibrate.py."""
    points = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        points.append((float(parts[0]), float(parts[1])))
    return points


def tail_stats(points: list[tuple[float, float]], tail_ps: float) -> tuple[float, float, int]:
    """Mean, population std dev, and sample count over the last tail_ps
    of the series (by time, not by index count, since sampling interval
    can vary)."""
    if not points:
        return float("nan"), float("nan"), 0
    t_end = points[-1][0]
    tail = [v for t, v in points if t >= t_end - tail_ps]
    if not tail:
        tail = [points[-1][1]]
    n = len(tail)
    mean = sum(tail) / n
    var = sum((v - mean) ** 2 for v in tail) / n if n > 1 else 0.0
    return mean, var ** 0.5, n


def extract_xvg(edr: Path, observable: str, out_xvg: Path) -> None:
    if not shutil.which("gmx"):
        sys.exit("ERROR: gmx not on PATH -- see docs/setup.md for the GPU build")
    proc = subprocess.run(
        ["gmx", "energy", "-f", str(edr), "-o", str(out_xvg)],
        input=f"{observable}\n0\n", capture_output=True, text=True,
    )
    if proc.returncode != 0 or not out_xvg.exists():
        sys.exit(f"ERROR: gmx energy failed for {observable}: {proc.stderr}")


def run_gmx_phase(gro: Path, top: Path, mdp: Path, out_dir: Path, deffnm: str) -> Path:
    if not shutil.which("gmx"):
        sys.exit("ERROR: gmx not on PATH -- see docs/setup.md for the GPU build")
    tpr = out_dir / f"{deffnm}.tpr"
    subprocess.run(
        ["gmx", "grompp", "-f", str(mdp), "-c", str(gro), "-p", str(top),
         "-o", str(tpr), "-po", str(out_dir / f"{deffnm}_mdout.mdp")],
        cwd=out_dir, check=True,
    )
    subprocess.run(
        ["gmx", "mdrun", "-deffnm", deffnm, "-nb", "gpu", "-pme", "gpu"],
        cwd=out_dir, check=True,
    )
    return out_dir / f"{deffnm}.edr"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["nvt", "npt"], required=True)
    parser.add_argument("--gro", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--mdp", required=True)
    parser.add_argument("--out-dir", default="../system_prep")
    parser.add_argument("--tail-ps", type=float, default=100.0,
                         help="stability window: check the last N ps of "
                              "the run, not the whole trace (early "
                              "equilibration transients are expected)")
    parser.add_argument("--tolerance", type=float, default=None,
                         help="override the default per-observable std-dev "
                              "tolerance for every observable in this phase")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--skip-run", action="store_true",
                         help="grade an existing --out-dir/<phase>.edr "
                              "instead of invoking gmx")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    deffnm = args.phase
    gl = GateLogger()
    gl.log("=" * 60)
    gl.log(f"{args.phase.upper()} EQUILIBRATION GATE")
    gl.log("=" * 60)

    edr_path = out_dir / f"{deffnm}.edr"
    if args.skip_run:
        gl.log(f"\n[1] Skipping run, grading existing {edr_path}")
        if not edr_path.exists():
            sys.exit(f"ERROR: --skip-run but {edr_path} doesn't exist")
    else:
        gl.log(f"\n[1] Running gmx grompp + mdrun (deffnm={deffnm})...")
        edr_path = run_gmx_phase(Path(args.gro), Path(args.top), Path(args.mdp),
                                  out_dir, deffnm)

    observables = PHASE_OBSERVABLES[args.phase]
    results = {}
    all_pass = True
    for obs in observables:
        xvg_path = out_dir / f"{deffnm}_{obs.lower()}.xvg"
        gl.log(f"\n[2] Extracting {obs}...")
        extract_xvg(edr_path, obs, xvg_path)
        points = parse_xvg(xvg_path)
        mean, std, n = tail_stats(points, args.tail_ps)
        tol = args.tolerance if args.tolerance is not None else DEFAULT_TOLERANCE[obs]
        obs_pass = std <= tol
        all_pass = all_pass and obs_pass
        results[obs] = {"mean": mean, "std": std, "n_samples": n, "tolerance": tol,
                         "pass": obs_pass}
        gl.log(f"  Last {args.tail_ps} ps: mean={mean:.3f}, std={std:.3f} "
               f"(n={n}) -- tolerance {tol} -> {'PASS' if obs_pass else 'FAIL'}")

    verdict = "PASS" if all_pass else "FAIL"
    gl.log(f"\n  GATE: {verdict}")
    gl.write_report(out_dir / f"{deffnm}_report.txt")
    write_gate_status(out_dir / f"{deffnm}_gate_status.json", verdict, {
        "phase": args.phase,
        "tail_ps": args.tail_ps,
        "observables": results,
    })

    enforce_strict(args.strict, verdict,
        f"GATE FAILED: {args.phase.upper()} did not stabilize within "
        f"tolerance over the last {args.tail_ps} ps. Refusing to treat "
        f"{edr_path} as equilibrated input for the next phase.")


if __name__ == "__main__":
    main()
