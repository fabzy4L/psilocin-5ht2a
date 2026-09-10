"""
Energy minimization stage: gmx grompp + gmx mdrun, gated on whether the
run actually converged to the target Fmax rather than just completing.

gmx's steepest-descent EM prints a line like:
    "Steepest Descents converged to Fmax < 1000 in 4213 steps"
or, if it did NOT converge before hitting nsteps:
    "Steepest Descents did not converge to Fmax < 1000 in 5000 steps"
followed later by lines including:
    "Potential Energy  = -1.2345678e+06"
    "Maximum force     =  8.234e+02 on atom 1234"
Parsing the log rather than trusting the .mdp's emtol at face value matters
because a run can hit -nsteps (or get restarted with a softer emtol)
without actually reaching the target force, and mdrun still exits 0.

Usage:
    python 02_minimize.py --gro ../system_prep/step5_input.gro \
        --top ../system_prep/topol.top --mdp ../system_prep/minimization.mdp \
        --out-dir ../system_prep --strict

Requires: gmx on PATH (see docs/setup.md for the GPU build).
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from gate_utils import GateLogger, write_gate_status, enforce_strict

FMAX_LINE = re.compile(r"Maximum force\s*=\s*([\d.eE+-]+)")
CONVERGED_LINE = re.compile(r"converged to Fmax < (\d+(?:\.\d+)?)")


def parse_final_max_force(log_text: str) -> float | None:
    """Returns the last 'Maximum force = ...' value in an EM log, or None
    if the log doesn't have one (run didn't get that far / wrong log
    passed in). Pure function -- unit-tested against a synthetic log
    fixture in tests/test_02_minimize.py, no gmx binary required."""
    matches = FMAX_LINE.findall(log_text)
    return float(matches[-1]) if matches else None


def parse_converged_flag(log_text: str) -> bool | None:
    """True if gmx's own convergence line says it converged, False if it
    explicitly says it did NOT, None if neither line is present."""
    if "did not converge" in log_text:
        return False
    if CONVERGED_LINE.search(log_text):
        return True
    return None


def run_gmx_em(gro: Path, top: Path, mdp: Path, out_dir: Path, deffnm: str) -> Path:
    if not shutil.which("gmx"):
        sys.exit("ERROR: gmx not on PATH -- see docs/setup.md for the GPU build")
    tpr = out_dir / f"{deffnm}.tpr"
    subprocess.run(
        ["gmx", "grompp", "-f", str(mdp), "-c", str(gro), "-p", str(top),
         "-o", str(tpr), "-po", str(out_dir / f"{deffnm}_mdout.mdp")],
        cwd=out_dir, check=True,
    )
    subprocess.run(
        ["gmx", "mdrun", "-deffnm", deffnm, "-nb", "gpu"],
        cwd=out_dir, check=True,
    )
    return out_dir / f"{deffnm}.log"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--mdp", required=True)
    parser.add_argument("--out-dir", default="../system_prep")
    parser.add_argument("--deffnm", default="em")
    parser.add_argument("--fmax-threshold", type=float, default=1000.0,
                         help="kJ/mol/nm, gate threshold checked "
                              "independently of whatever emtol the .mdp "
                              "used -- catches a run that 'completed' via "
                              "-nsteps without actually converging")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--skip-run", action="store_true",
                         help="grade an existing --out-dir/<deffnm>.log "
                              "instead of invoking gmx (for re-grading, or "
                              "testing the gate logic without gmx installed)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gl = GateLogger()
    gl.log("=" * 60)
    gl.log("ENERGY MINIMIZATION GATE")
    gl.log("=" * 60)

    log_path = out_dir / f"{args.deffnm}.log"
    if args.skip_run:
        gl.log(f"\n[1] Skipping run, grading existing {log_path}")
        if not log_path.exists():
            sys.exit(f"ERROR: --skip-run but {log_path} doesn't exist")
    else:
        gl.log(f"\n[1] Running gmx grompp + mdrun (deffnm={args.deffnm})...")
        log_path = run_gmx_em(Path(args.gro), Path(args.top), Path(args.mdp),
                               out_dir, args.deffnm)

    log_text = log_path.read_text(errors="replace")
    fmax = parse_final_max_force(log_text)
    converged = parse_converged_flag(log_text)

    gl.log(f"\n[2] Parsing {log_path.name}...")
    gl.log(f"  gmx-reported converged: {converged}")
    gl.log(f"  Final max force:        {fmax} kJ/mol/nm" if fmax is not None
           else "  Final max force:        could not parse")

    if fmax is None:
        verdict = "FAIL"
        gl.log("  Could not extract a max-force value from the log.")
    else:
        verdict = "PASS" if fmax < args.fmax_threshold else "FAIL"
        gl.log(f"  Threshold: <{args.fmax_threshold} kJ/mol/nm -> {verdict}")

    gl.log(f"\n  GATE: {verdict}")
    gl.write_report(out_dir / f"{args.deffnm}_report.txt")
    write_gate_status(out_dir / f"{args.deffnm}_gate_status.json", verdict, {
        "fmax_kj_per_mol_nm": fmax,
        "fmax_threshold": args.fmax_threshold,
        "gmx_reported_converged": converged,
        "log_file": str(log_path),
    })

    enforce_strict(args.strict, verdict,
        f"GATE FAILED: {log_path} did not converge below "
        f"{args.fmax_threshold} kJ/mol/nm (fmax={fmax}). Refusing to treat "
        f"this minimized structure as safe input for equilibration.")


if __name__ == "__main__":
    main()
