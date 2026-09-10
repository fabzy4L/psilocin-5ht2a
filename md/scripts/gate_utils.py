"""
Shared helpers for md/scripts/ pipeline stages: consistent gate output
(report.txt + gate_status.json) and --strict exit behavior.

Factored out after docking/scripts/05_validation_redock.py's RMSD bug (see
docs/REDOCK_BUG_HANDOFF.md) showed the cost of a gate re-implementing its
own report/parsing logic ad hoc with no test coverage -- one tested
implementation of the report/gate plumbing here is safer than five
duplicated copies of it across the MD stages.
"""
import json
import sys
from pathlib import Path


class GateLogger:
    """Accumulates log lines for a *_report.txt while also printing them."""

    def __init__(self):
        self.lines = []

    def log(self, msg: str = "") -> None:
        print(msg)
        self.lines.append(msg)

    def write_report(self, path: Path) -> None:
        path.write_text("\n".join(self.lines))


def write_gate_status(path: Path, gate: str, criteria: dict) -> None:
    """criteria: dict of metric_name -> value (numbers/bools/strings/None).
    `gate` should always be "PASS" or "FAIL"."""
    payload = {"gate": gate}
    payload.update(criteria)
    path.write_text(json.dumps(payload, indent=2, default=str))


def enforce_strict(strict: bool, gate: str, message: str) -> None:
    """Exit non-zero if --strict was passed and the gate isn't PASS."""
    if strict and gate != "PASS":
        sys.exit(message)
