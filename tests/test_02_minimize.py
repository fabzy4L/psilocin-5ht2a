"""
This gate exists specifically to catch a minimization that "completed"
(gmx mdrun exit 0) without actually converging -- these fixtures cover
both the converged and did-not-converge log shapes gmx actually emits.
"""
from conftest import load_script

CONVERGED_LOG = """\
Steepest Descents:
   Tolerance (Fmax)   =  1.00000e+03
   Number of steps    =         5000

Steepest Descents converged to Fmax < 1000 in 4213 steps
Potential Energy  = -1.23456780e+06
Maximum force     =  8.23400e+02 on atom 1234
Norm of force     =  1.10000e+01
"""

NOT_CONVERGED_LOG = """\
Steepest Descents:
   Tolerance (Fmax)   =  1.00000e+03
   Number of steps    =         5000

Steepest Descents did not converge to Fmax < 1000 in 5000 steps.
Potential Energy  = -1.10000000e+06
Maximum force     =  3.45000e+03 on atom 9876
Norm of force     =  4.20000e+01
"""

NO_FORCE_LINE_LOG = "some unrelated gmx log output with no force line\n"


def test_parse_final_max_force_converged():
    mod = load_script("02_minimize.py")
    assert mod.parse_final_max_force(CONVERGED_LOG) == 823.4


def test_parse_final_max_force_not_converged():
    mod = load_script("02_minimize.py")
    assert mod.parse_final_max_force(NOT_CONVERGED_LOG) == 3450.0


def test_parse_final_max_force_missing_returns_none():
    mod = load_script("02_minimize.py")
    assert mod.parse_final_max_force(NO_FORCE_LINE_LOG) is None


def test_parse_converged_flag_true():
    mod = load_script("02_minimize.py")
    assert mod.parse_converged_flag(CONVERGED_LOG) is True


def test_parse_converged_flag_false():
    mod = load_script("02_minimize.py")
    assert mod.parse_converged_flag(NOT_CONVERGED_LOG) is False


def test_parse_converged_flag_none_when_absent():
    mod = load_script("02_minimize.py")
    assert mod.parse_converged_flag(NO_FORCE_LINE_LOG) is None


def test_skip_run_gate_passes_below_threshold(tmp_path):
    """End-to-end through main()'s gating logic via --skip-run, no gmx
    binary required: a converged log under the default 1000 threshold
    should write a PASS gate_status.json."""
    import json
    import sys

    mod = load_script("02_minimize.py")
    (tmp_path / "em.log").write_text(CONVERGED_LOG)

    argv = ["02_minimize.py", "--gro", "x.gro", "--top", "x.top", "--mdp", "x.mdp",
            "--out-dir", str(tmp_path), "--skip-run"]
    old_argv = sys.argv
    sys.argv = argv
    try:
        mod.main()
    finally:
        sys.argv = old_argv

    status = json.loads((tmp_path / "em_gate_status.json").read_text())
    assert status["gate"] == "PASS"
    assert status["fmax_kj_per_mol_nm"] == 823.4


def test_skip_run_gate_fails_above_threshold(tmp_path):
    import json
    import sys

    mod = load_script("02_minimize.py")
    (tmp_path / "em.log").write_text(NOT_CONVERGED_LOG)

    argv = ["02_minimize.py", "--gro", "x.gro", "--top", "x.top", "--mdp", "x.mdp",
            "--out-dir", str(tmp_path), "--skip-run"]
    old_argv = sys.argv
    sys.argv = argv
    try:
        mod.main()
    finally:
        sys.argv = old_argv

    status = json.loads((tmp_path / "em_gate_status.json").read_text())
    assert status["gate"] == "FAIL"
    assert status["fmax_kj_per_mol_nm"] == 3450.0


def test_skip_run_strict_exits_on_fail(tmp_path):
    import sys

    import pytest

    mod = load_script("02_minimize.py")
    (tmp_path / "em.log").write_text(NOT_CONVERGED_LOG)

    argv = ["02_minimize.py", "--gro", "x.gro", "--top", "x.top", "--mdp", "x.mdp",
            "--out-dir", str(tmp_path), "--skip-run", "--strict"]
    old_argv = sys.argv
    sys.argv = argv
    try:
        with pytest.raises(SystemExit):
            mod.main()
    finally:
        sys.argv = old_argv
