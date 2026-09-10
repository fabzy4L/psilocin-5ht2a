"""
Production MD run driver: long GPU-offloaded gmx mdrun, using the
GPU-offload flags already scoped in docs/setup.md. Not gated pass/fail
like the earlier stages -- there's no single "did production work"
numeric criterion at this point, that's what 05_analysis.py's RMSD-plateau
gate is for. This script just runs it (or a short --smoke-test slice) and
reminds you to check GPU utilization before committing to the full
100-200 ns, per docs/setup.md's own advice.

Usage:
    # short GPU-utilization sanity check first (per docs/setup.md):
    python 04_production.py --gro ../system_prep/npt.gro \
        --top ../system_prep/topol.top --mdp ../system_prep/production.mdp \
        --out-dir ../system_prep --smoke-test

    # then the real run:
    python 04_production.py --gro ../system_prep/npt.gro \
        --top ../system_prep/topol.top --mdp ../system_prep/production.mdp \
        --out-dir ../system_prep

Requires: gmx on PATH, built with CUDA support (docs/setup.md).
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

SMOKE_TEST_STEPS = 5000  # a few ps at a typical 2 fs timestep


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", required=True)
    parser.add_argument("--top", required=True)
    parser.add_argument("--mdp", required=True)
    parser.add_argument("--out-dir", default="../system_prep")
    parser.add_argument("--deffnm", default="md")
    parser.add_argument(
        "--smoke-test", action="store_true",
        help=f"override the .mdp's nsteps down to {SMOKE_TEST_STEPS} (a "
             f"few ps) and run that instead of full production -- use "
             f"this first with nvidia-smi open in another terminal to "
             f"confirm GPU utilization before committing to a multi-day "
             f"run, per docs/setup.md")
    parser.add_argument("--ntomp", type=int, default=8)
    args = parser.parse_args()

    if not shutil.which("gmx"):
        sys.exit("ERROR: gmx not on PATH -- see docs/setup.md for the GPU build")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tpr = out_dir / f"{args.deffnm}.tpr"

    grompp_cmd = ["gmx", "grompp", "-f", str(args.mdp), "-c", str(args.gro),
                  "-p", str(args.top), "-o", str(tpr),
                  "-po", str(out_dir / f"{args.deffnm}_mdout.mdp")]
    if args.smoke_test:
        print(f"[SMOKE TEST] overriding nsteps to {SMOKE_TEST_STEPS} -- "
              f"watch `nvidia-smi` in another terminal now")
    subprocess.run(grompp_cmd, cwd=out_dir, check=True)

    mdrun_cmd = ["gmx", "mdrun", "-deffnm", args.deffnm,
                 "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu",
                 "-ntmpi", "1", "-ntomp", str(args.ntomp)]
    if args.smoke_test:
        mdrun_cmd += ["-nsteps", str(SMOKE_TEST_STEPS)]
    subprocess.run(mdrun_cmd, cwd=out_dir, check=True)

    print(f"\nDone. {'Smoke test' if args.smoke_test else 'Production'} "
          f"run output: {out_dir / (args.deffnm + '.xtc')}")
    if not args.smoke_test:
        print("Reminder (docs/setup.md): expect ~20-50 ns/day for a "
              "50-100k atom system on an 8GB card -- if this run's rate "
              "looks far off that, check nvidia-smi before assuming it's "
              "just slow.")


if __name__ == "__main__":
    main()
